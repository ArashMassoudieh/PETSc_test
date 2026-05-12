// script2_executor.cpp
//
// Executes an S2Program by driving Grid2D objects directly.
// No dependency on sim_runner.
// -----------------------------------------------------------------------

#include "script2_executor.h"
#include "script2_parser.h"
#include "script2_types.h"

#include "grid.h"
#include "sim_helpers.h"    // createDirectory, joinPath

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <map>
#include <memory>
#include <algorithm>
#include <functional>
#include "copula_analysis.h"
#include "TimeSeriesSet.h"

// -----------------------------------------------------------------------
// Runtime state: named grids and accumulators
// -----------------------------------------------------------------------
namespace {

struct RuntimeState
{
    S2Context                                                  ctx;
    std::map<std::string, std::unique_ptr<Grid2D>>             grids;
    std::map<std::string, std::unique_ptr<CopulaAccumulator>>  accums;
    std::map<std::string, std::unique_ptr<BTCAccumulator>>     btc_accums;
    // Mean CDF accumulator: sum of per-realization CDFs (vector of (u,v) pairs)
    // keyed by a label the user assigns with g.extract_cdf output=NAME
    std::map<std::string, std::vector<TimeSeries<double>>>     cdf_sets;

    // ---- Loaded objects (from disk, e.g. precomputed mean copulas/CDFs) ----
    std::map<std::string, CMatrix>            copulas;   // <-- added
    std::map<std::string, TimeSeries<double>> cdfs;      // <-- added

    int rank = 0;

    Grid2D& grid(const std::string& name) {
        auto it = grids.find(name);
        if (it == grids.end())
            throw std::runtime_error("grid '" + name + "' not found — "
                                                       "declare it with: grid " + name + " = nx:...");
        return *it->second;
    }
    CopulaAccumulator& accum(const std::string& name) {
        auto it = accums.find(name);
        if (it == accums.end())
            throw std::runtime_error("accumulator '" + name + "' not found — "
                                                              "declare it with: accumulator " + name + " = bins:...");
        return *it->second;
    }
    BTCAccumulator& btc(const std::string& name) {
        auto it = btc_accums.find(name);
        if (it == btc_accums.end())
            throw std::runtime_error("btc accumulator '" + name + "' not found");
        return *it->second;
    }

    // ---- Accessors for loaded objects ----
    CMatrix& copula(const std::string& name) {                              // <-- added
        auto it = copulas.find(name);
        if (it == copulas.end())
            throw std::runtime_error("copula '" + name + "' not loaded — "
                                                         "load it with: g.load_copula name=" + name + ", file=...");
        return it->second;
    }
    TimeSeries<double>& cdf(const std::string& name) {                      // <-- added
        auto it = cdfs.find(name);
        if (it == cdfs.end())
            throw std::runtime_error("cdf '" + name + "' not loaded — "
                                                      "load it with: g.load_cdf name=" + name + ", file=...");
        return it->second;
    }
};
// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

// Resolve ArrayKind from field name: qx/qy → Fx/Fy, everything else → Cell
Grid2D::ArrayKind kindOf(const std::string& name, const Grid2D& g)
{
    // Check fluxes first
    if (g.hasFlux(name)) {
        // qx-family → Fx, qy-family → Fy
        const std::string low = [&]{ std::string s=name;
            std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c){return std::tolower(c);});
            return s; }();
        if (low.find("qy") != std::string::npos ||
            low.find("fy") != std::string::npos)
            return Grid2D::ArrayKind::Fy;
        return Grid2D::ArrayKind::Fx;
    }
    return Grid2D::ArrayKind::Cell;
}

// Parse a "min:max:n" or "min:max" range spec
// Returns {min, max, n_points}
struct RangeSpec { double lo; double hi; int n; };
RangeSpec parseRange(const std::string& s, int lineno)
{
    std::vector<std::string> tok;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ':')) {
        t.erase(0, t.find_first_not_of(" \t"));
        t.erase(t.find_last_not_of(" \t") + 1);
        if (!t.empty()) tok.push_back(t);
    }
    if (tok.size() < 2)
        throw std::runtime_error("line " + std::to_string(lineno)
            + ": range needs at least min:max, got '" + s + "'");
    RangeSpec r;
    r.lo = std::stod(tok[0]);
    r.hi = std::stod(tok[1]);
    r.n  = (tok.size() >= 3) ? std::stoi(tok[2]) : 10;
    if (r.n < 2) r.n = 2;
    return r;
}

// Ensure directory exists (rank 0 only)
void ensureDir(const std::string& path, int rank)
{
    if (rank == 0) createDirectory(path);
    MPI_Barrier(PETSC_COMM_WORLD);
}

// -----------------------------------------------------------------------
// Build a binned copula matrix from a set of (u1, u2) rank pairs.
// pairs: TimeSeries where .t = u1, .c = u2, both in [0,1].
// Returns flat row-major n_bins x n_bins matrix, values = count/total/n_bins
// (so each row sums to 1/n_bins when well-sampled).
// -----------------------------------------------------------------------
std::vector<double> binnedCopula(const TimeSeries<double>& pairs, int n_bins)
{
    std::vector<double> mat(n_bins * n_bins, 0.0);
    int total = 0;
    for (int k = 0; k < (int)pairs.size(); ++k) {
        double u1 = pairs.getTime(k);
        double u2 = pairs.getValue(k);
        // Clamp to [0,1)
        u1 = std::max(0.0, std::min(1.0 - 1e-12, u1));
        u2 = std::max(0.0, std::min(1.0 - 1e-12, u2));
        int r = (int)(u1 * n_bins);
        int c = (int)(u2 * n_bins);
        mat[r * n_bins + c] += 1.0;
        ++total;
    }
    if (total > 0) {
        const double du = 1.0 / n_bins;
        for (auto& v : mat) v /= (total * du);  // normalise so rows sum to du
    }
    return mat;
}

// Write binned copula to CSV  (n_bins rows, n_bins comma-separated values)
void writeCopulaCSV(const std::vector<double>& mat, int n_bins,
                    const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("writeCopulaCSV: cannot open '" + path + "'");
    f << std::scientific << std::setprecision(6);
    for (int r = 0; r < n_bins; ++r) {
        for (int c = 0; c < n_bins; ++c) {
            f << mat[r * n_bins + c];
            if (c < n_bins - 1) f << ", ";
        }
        f << "\n";
    }
}


// -----------------------------------------------------------------------
// Forward declaration
// -----------------------------------------------------------------------
bool execStmts(const std::vector<S2Stmt>& stmts, RuntimeState& st);

// -----------------------------------------------------------------------
// Dispatch a single method call:  obj.method  args...
// -----------------------------------------------------------------------
bool execMethod(const S2Stmt& stmt, RuntimeState& st)
{
    S2Context&         ctx    = st.ctx;
    const std::string  obj    = ctx.expand(stmt.obj);
    const std::string& method = stmt.method;
    const auto&        args   = stmt.args;
    const int          rank   = st.rank;
    const int          lineno = stmt.lineno;

    // ---- Accumulator methods ----
    if (st.accums.count(obj)) {
        CopulaAccumulator& acc = st.accum(obj);

        // accum.save_mean  file=path
        if (method == "save_mean") {
            if (rank == 0) {
                const std::string path = ctx.expand(ctx.argVal(args, "file", obj + "_mean.csv"));
                ensureDir(path.substr(0, path.find_last_of("/\\")), rank);
                acc.saveMean(path);
            }
            MPI_Barrier(PETSC_COMM_WORLD);
            return true;
        }

        // accum.save_mean_vti  file=path
        if (method == "save_mean_vti") {
            if (rank == 0) {
                const std::string path = ctx.expand(ctx.argVal(args, "file", obj + "_mean.vti"));
                ensureDir(path.substr(0, path.find_last_of("/\\")), rank);
                acc.saveMeanVTI(path);
            }
            MPI_Barrier(PETSC_COMM_WORLD);
            return true;
        }

        throw std::runtime_error("line " + std::to_string(lineno)
            + ": unknown accumulator method '" + method + "'");
    }

    // ---- BTC accumulator methods ----
    if (st.btc_accums.count(obj)) {
        BTCAccumulator& bacc = st.btc(obj);

        if (method == "save_mean") {
            if (rank == 0) {
                const std::string path = ctx.expand(ctx.argVal(args, "file", obj + "_mean.csv"));
                ensureDir(path.substr(0, path.find_last_of("/\\")), rank);
                bacc.saveMean(path);
                std::cout << "  " << obj << ".save_mean ("
                          << bacc.count() << " realizations) -> " << path << "\n";
            }
            MPI_Barrier(PETSC_COMM_WORLD);
            return true;
        }

        throw std::runtime_error("line " + std::to_string(lineno)
                                 + ": unknown btc-accumulator method '" + method + "'");
    }


    // ---- Grid methods ----
    Grid2D& g = st.grid(obj);

    // ------------------------------------------------------------------
    // g.create_field  name=K_ns, type=gaussian_sgs,
    //                 lx=$corr_x, ly=$corr_y, stdev=1.0, seed=$seed
    // ------------------------------------------------------------------
    if (method == "create_field") {
        const std::string name  = ctx.argVal(args, "name", "K_ns");
        const std::string type  = ctx.argVal(args, "type", "gaussian_sgs");
        const double lx         = ctx.argDouble(args, "lx", 1.0);
        const double ly         = ctx.argDouble(args, "ly", 0.1);
        const double stdev      = ctx.argDouble(args, "stdev", 1.0);
        const unsigned long seed= static_cast<unsigned long>(
                                    ctx.argInt(args, "seed", 42));
        const double nugget     = ctx.argDouble(args, "nugget", 0.0);
        const int    nneigh     = ctx.argInt(args, "nneigh", 10);

        if (type == "gaussian_sgs" || type == "sgs") {
            g.makeGaussianFieldSGS(name, lx, ly, nneigh, seed, nugget);
        } else if (type == "gaussian_fft" || type == "fft") {
            g.generateGaussianRandomField(name, lx, ly, 0.0, stdev * stdev, seed);
        } else {
            throw std::runtime_error("line " + std::to_string(lineno)
                + ": unknown field type '" + type + "'");
        }
        if (rank == 0)
            std::cout << "  create_field " << name << " (" << type
                      << ", lx=" << lx << ", ly=" << ly << ")\n";
        return true;
    }

    // ------------------------------------------------------------------
    // g.to_lognormal  input=K_ns, stdev=$stdev, mean=0.0, output=K
    // ------------------------------------------------------------------
    if (method == "to_lognormal") {
        const std::string input  = ctx.argVal(args, "input",  "K_ns");
        const std::string output = ctx.argVal(args, "output", "K");
        const double stdev       = ctx.argDouble(args, "stdev", 1.0);
        const double mean        = ctx.argDouble(args, "mean",  0.0);
        g.normalizeField(input, Grid2D::ArrayKind::Cell);
        g.createExponentialField(input, stdev, mean, output);
        if (rank == 0)
            std::cout << "  to_lognormal " << input << " -> " << output
                      << " (stdev=" << stdev << ")\n";
        return true;
    }

    // ------------------------------------------------------------------
    // g.solve_darcy  Kx=K, Ky=K, H_west=1.0, H_east=0.0
    // ------------------------------------------------------------------
    if (method == "solve_darcy") {
        const std::string kx    = ctx.argVal(args, "Kx", "K");
        const std::string ky    = ctx.argVal(args, "Ky", "K");
        const double h_west     = ctx.argDouble(args, "H_west", 1.0);
        const double h_east     = ctx.argDouble(args, "H_east", 0.0);
        g.DarcySolve(h_west, h_east, kx, ky);
        if (rank == 0) std::cout << "  solve_darcy done\n";
        return true;
    }

    // ------------------------------------------------------------------
    // g.to_ranks         input=qx, output=qx_ranks
    // g.to_normal_score  input=qx, output=qx_ns
    // ------------------------------------------------------------------
    if (method == "to_ranks" || method == "to_normal_score") {
        const std::string input  = ctx.argVal(args, "input",  "qx");
        const std::string output = ctx.argVal(args, "output",
                                              method == "to_ranks" ? "qx_ranks" : "qx_ns");
        const auto kind = kindOf(input, g);
        TimeSeries<double> values = g.exportFieldToTimeSeries(input, kind);
        TimeSeries<double> result = (method == "to_ranks")
                                    ? values.ConvertToRanks()
                                    : values.ConvertToNormalScore();
        g.assignFromTimeSeries(result, output, kind);
        if (rank == 0)
            std::cout << "  " << method << " " << input << " -> " << output << "\n";
        return true;
    }

    // ------------------------------------------------------------------
    // g.assign_constant  field=C, kind=cell, value=0.0
    //                    (kind ∈ {cell, fx, fy}; default cell)
    // Initializes or resets a field to a constant value.
    // ------------------------------------------------------------------
    if (method == "assign_constant") {
        const std::string fname = ctx.argVal(args, "field", "");
        if (fname.empty())
            throw std::runtime_error("line " + std::to_string(lineno)
                                     + ": assign_constant requires field=NAME");

        const std::string kindstr = ctx.argVal(args, "kind", "cell");
        Grid2D::ArrayKind kind =
            (kindstr == "fx" || kindstr == "Fx") ? Grid2D::ArrayKind::Fx :
                (kindstr == "fy" || kindstr == "Fy") ? Grid2D::ArrayKind::Fy :
                Grid2D::ArrayKind::Cell;

        const double value = ctx.argDouble(args, "value", 0.0);

        g.assignConstant(fname, kind, value);
        if (rank == 0)
            std::cout << "  assign_constant " << fname
                      << " (kind=" << kindstr << ") = " << value << "\n";
        return true;
    }

    // ------------------------------------------------------------------
    // g.extract_copula  field=qx_ranks, direction=x,
    //                   range=0.001:1.0:20, samples=5000,
    //                   accumulate=adv_copula,
    //                   save=output/r$i/adv_copula_x.csv
    //
    // Samples pairs of velocity ranks at log-spaced separations,
    // bins them into an empirical copula, and optionally accumulates
    // the mean across realizations.
    //
    // The copula at the representative delta (geometric mean of range)
    // is what gets saved and accumulated.  Correlation vs. delta is
    // written to save_corr= if provided.
    // ------------------------------------------------------------------
    if (method == "extract_copula") {
        const std::string field  = ctx.argVal(args, "field", "qx_ranks");
        const std::string dirstr = ctx.argVal(args, "direction", "x");
        const std::string rspec  = ctx.argVal(args, "range", "0.001:1.0:20");
        const int n_samples      = ctx.argInt(args, "samples", 5000);
        const int bins           = ctx.argInt(args, "bins", 20);
        const std::string accname= ctx.argVal(args, "accumulate", "");
        const std::string savef  = ctx.argVal(args, "save", "");
        const std::string savecorr=ctx.argVal(args, "save_corr", "");

        const auto kind  = kindOf(field, g);
        const auto rng   = parseRange(rspec, lineno);

        PerturbDir dir =
            (dirstr == "y" || dirstr == "Y") ? PerturbDir::YOnly  :
                (dirstr == "r" || dirstr == "R") ? PerturbDir::Radial :
                PerturbDir::XOnly;

        // Representative delta = geometric mean of range; this is the
        // distance at which the saved/accumulated copula is built.
        // (When lo == hi, this is just that single distance.)
        const double delta_rep = std::sqrt(rng.lo * rng.hi);

        // Sweep log-spaced deltas across the range; record correlation
        // vs delta for save_corr, and capture the copula at the delta
        // closest to delta_rep in log space.
        TimeSeries<double> corr_ts;
        std::vector<double> copula_mat;
        double best_log_dist = 1e30;

        for (int k = 0; k < rng.n; ++k) {
            const double exponent = (rng.n == 1) ? 0.0
                                                 : static_cast<double>(k) / (rng.n - 1);
            const double delta    = rng.lo * std::pow(rng.hi / rng.lo, exponent);

            try {
                TimeSeries<double> pairs = g.sampleGaussianPerturbation(
                    field, kind, n_samples, delta, 0, dir);

                corr_ts.append(delta, pairs.correlation_tc());

                const double ldist = std::abs(std::log(delta) - std::log(delta_rep));
                if (ldist < best_log_dist) {
                    best_log_dist = ldist;
                    copula_mat    = binnedCopula(pairs, bins);
                }
            } catch (...) {
                if (rank == 0)
                    std::cerr << "  [warn] extract_copula: failed at delta=" << delta << "\n";
            }
        }

        if (rank == 0) {
            std::cout << "  extract_copula " << field << " dir=" << dirstr
                      << " range=[" << rng.lo << ":" << rng.hi << ":" << rng.n << "]"
                      << " bins=" << bins
                      << "  delta_rep=" << delta_rep << "\n";

            if (!savef.empty()) {
                ensureDir(savef.substr(0, savef.find_last_of("/\\")), rank);
                if (!copula_mat.empty())
                    writeCopulaCSV(copula_mat, bins, savef);
            }

            if (!savecorr.empty()) {
                ensureDir(savecorr.substr(0, savecorr.find_last_of("/\\")), rank);
                corr_ts.writefile(savecorr);
            }
        }

        // Accumulate into the named accumulator, auto-creating it on
        // first use.  This makes per-distance accumulator names inside
        // a 'foreach' loop work without requiring a separate
        // 'accumulator NAME = bins:N' declaration for each one.
        if (!accname.empty() && !copula_mat.empty()) {
            if (st.accums.find(accname) == st.accums.end()) {
                st.accums[accname] = std::make_unique<CopulaAccumulator>(bins);
                if (rank == 0)
                    std::cout << "  [auto] created accumulator '" << accname
                              << "' (bins=" << bins << ")\n";
            }
            st.accum(accname).add(copula_mat);
        }

        return true;
    }

    // ------------------------------------------------------------------
    // g.extract_cdf  field=qx, output=qx_cdf, save=path/qx_cdf.txt
    //                bins=100, threshold=1e-6
    // ------------------------------------------------------------------
    if (method == "extract_cdf") {
        const std::string field  = ctx.argVal(args, "field", "qx");
        const std::string outvar = ctx.argVal(args, "output", "");
        const std::string savef  = ctx.argVal(args, "save", "");
        const int bins           = ctx.argInt(args, "bins", 100);
        const double thresh      = ctx.argDouble(args, "threshold", 1e-6);
        const auto kind          = kindOf(field, g);

        TimeSeries<double> cdf = g.extractFieldCDF(field, kind, bins, thresh);

        if (rank == 0) {
            std::cout << "  extract_cdf " << field << " -> ";
            if (!savef.empty()) {
                ensureDir(savef.substr(0, savef.find_last_of("/\\")), rank);
                cdf.writefile(savef);
                std::cout << savef;
            }
            std::cout << "\n";
        }

        // Accumulate into cdf_sets for later mean computation
        if (!outvar.empty())
            st.cdf_sets[outvar].push_back(cdf);

        return true;
    }

    // ------------------------------------------------------------------
    // g.save_mean_cdf  cdf=qx_cdf, file=path/mean_qx_cdf.txt
    // ------------------------------------------------------------------
    if (method == "save_mean_cdf") {
        const std::string cdfname = ctx.argVal(args, "cdf", "qx_cdf");
        const std::string savef   = ctx.argVal(args, "file", "mean_qx_cdf.txt");
        if (rank == 0) {
            auto it = st.cdf_sets.find(cdfname);
            if (it == st.cdf_sets.end() || it->second.empty()) {
                std::cerr << "[warn] save_mean_cdf: no CDFs accumulated for '"
                          << cdfname << "'\n";
                return true;
            }
            const auto& cdfs = it->second;

            // Pack into a TimeSeriesSet, treating u as "time"
            TimeSeriesSet<double> tss;
            for (const auto& ts : cdfs) tss.append(ts, "");

            // Reference-grid interpolation (handles non-uniform u-grids,
            // skips realizations whose u-support doesn't cover a given u)
            TimeSeries<double> mean_cdf =
                tss.mean_ts_longest(/*start_item=*/0, /*time_eps=*/1e-12);

            ensureDir(savef.substr(0, savef.find_last_of("/\\")), rank);
            mean_cdf.writefile(savef);
            std::cout << "  save_mean_cdf " << cdfname << " ("
                      << cdfs.size() << " realizations) -> "
                      << savef << "  (" << mean_cdf.size() << " points)\n";
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        return true;
    }

    // load_copula
    if (method == "load_copula") {
        const std::string name = ctx.argVal(args, "name", "");
        const std::string file = ctx.argVal(args, "file", "");
        if (name.empty() || file.empty()) {
            std::cerr << "[err] load_copula: need name and file\n";
            return false;
        }
        if (rank == 0) {
            try {
                CMatrix M = CMatrix::readCSV(file);
                Grid2D::checkNormalisation(M, name);   // warns on noisy rows
                st.copulas[name] = std::move(M);
                std::cout << "  load_copula " << name << " ("
                          << st.copulas[name].getnumrows() << "x"
                          << st.copulas[name].getnumcols()
                          << ") <- " << file << "\n";
            } catch (const std::exception& e) {
                std::cerr << "[err] load_copula: " << e.what() << "\n";
                return false;
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        return true;
    }

    // load_cdf
    if (method == "load_cdf") {
        const std::string name  = ctx.argVal(args, "name", "");
        const std::string file  = ctx.argVal(args, "file", "");
        const double      u_cap = std::stod(ctx.argVal(args, "u_cap", "0.95"));
        if (name.empty() || file.empty()) {
            std::cerr << "[err] load_cdf: need name and file\n";
            return false;
        }
        if (rank == 0) {
            TimeSeries<double> cdf;
            if (!cdf.readfile(file)) {
                std::cerr << "[err] load_cdf: failed to read " << file << "\n";
                return false;
            }
            // Apply u_cap by truncating at the last point with u <= u_cap
            TimeSeries<double> capped;
            double v_cap = 0.0;
            for (int i = 0; i < (int)cdf.size(); ++i) {
                const double u = cdf.getTime(i);
                const double v = cdf.getValue(i);
                if (u <= u_cap) { capped.append(u, v); v_cap = v; }
                else break;
            }
            capped.append(u_cap, v_cap);  // sentinel: v(u >= u_cap) = v_cap
            st.cdfs[name] = std::move(capped);
            std::cout << "  load_cdf " << name << " ("
                      << st.cdfs[name].size() << " pts, u_cap=" << u_cap
                      << ") <- " << file << "\n";
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        return true;
    }
    if (method == "solve_copula_transport") {
        const std::string adv_name   = ctx.argVal(args, "theta_adv",  "");
        const std::string diff_name  = ctx.argVal(args, "theta_diff", "");
        const std::string cdf_name   = ctx.argVal(args, "v_of_u",     "");
        const double lambda_a        = std::stod(ctx.argVal(args, "lambda_a",  "1.0"));
        const double lambda_d        = std::stod(ctx.argVal(args, "lambda_d",  "0.1"));
        const double diffusion       = std::stod(ctx.argVal(args, "diffusion", "0.0"));
        const double t_end           = std::stod(ctx.argVal(args, "t_end",     "10.0"));
        const double dt              = std::stod(ctx.argVal(args, "dt",        "0.01"));
        const double c_left          = std::stod(ctx.argVal(args, "c_left",    "1.0"));
        const int    output_interval = std::stoi(ctx.argVal(args, "output_interval", "0"));
        const std::string output_dir = ctx.argVal(args, "output_dir", "");
        const std::string save_btc   = ctx.argVal(args, "save_btc",   "");
        const std::string btc_x_str  = ctx.argVal(args, "btc_x",      "");

        // Look up the named objects
        auto it_adv  = st.copulas.find(adv_name);
        auto it_diff = st.copulas.find(diff_name);
        auto it_cdf  = st.cdfs.find(cdf_name);
        if (it_adv == st.copulas.end()) {
            std::cerr << "[err] solve_copula_transport: theta_adv '"
                      << adv_name << "' not loaded\n";
            return false;
        }
        if (it_diff == st.copulas.end()) {
            std::cerr << "[err] solve_copula_transport: theta_diff '"
                      << diff_name << "' not loaded\n";
            return false;
        }
        if (it_cdf == st.cdfs.end()) {
            std::cerr << "[err] solve_copula_transport: v_of_u '"
                      << cdf_name << "' not loaded\n";
            return false;
        }

        const CMatrix& theta_adv  = it_adv->second;
        const CMatrix& theta_diff = it_diff->second;
        const TimeSeries<double>& cdf = it_cdf->second;

        // Build v(u) by linear interpolation of the loaded CDF, with tail cap
        auto v_of_u = [&cdf](double u) -> double {
            const int n = (int)cdf.size();
            if (n == 0) return 0.0;
            if (u <= cdf.getTime(0))     return cdf.getValue(0);
            if (u >= cdf.getTime(n - 1)) return cdf.getValue(n - 1);
            // Binary search for the bracketing interval
            int lo = 0, hi = n - 1;
            while (hi - lo > 1) {
                int mid = (lo + hi) / 2;
                if (cdf.getTime(mid) <= u) lo = mid; else hi = mid;
            }
            const double u0 = cdf.getTime(lo),  u1 = cdf.getTime(hi);
            const double v0 = cdf.getValue(lo), v1 = cdf.getValue(hi);
            const double a = (u - u0) / (u1 - u0);
            return (1.0 - a) * v0 + a * v1;
        };

        // Set diffusion coefficient on the grid (used by SolveCopulaTransport
        // for the x-direction diffusion AND for computing r_diff = 2D/lambda_d^2)
        Grid2D& g = st.grid(obj);
        g.SetVal("diffusion",diffusion);
        g.SetVal("c_left",c_left);

        // Parse BTC locations
        if (!btc_x_str.empty()) {
            std::vector<double> locs;
            std::stringstream ss(btc_x_str);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                const auto first = tok.find_first_not_of(" \t");
                const auto last  = tok.find_last_not_of(" \t");
                if (first == std::string::npos) continue;
                tok = tok.substr(first, last - first + 1);
                if (!tok.empty()) locs.push_back(std::stod(tok));
            }
            g.setBTCLocations(locs);
        }

        // Set up output
        if (!output_dir.empty()) ensureDir(output_dir, rank);

        TimeSeriesSet<double> btc_data;
        g.SolveCopulaTransport(
            t_end, dt,
            theta_adv, theta_diff,
            lambda_a, lambda_d,
            v_of_u,
            "upscaled_",                          // ksp_prefix
            save_btc.empty() ? nullptr : &btc_data,
            output_interval,
            output_dir);

        if (!save_btc.empty() && rank == 0) {
            ensureDir(save_btc.substr(0, save_btc.find_last_of("/\\")), rank);
            btc_data.write(save_btc);
            std::cout << "  saved upscaled BTC -> " << save_btc << "\n";
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        return true;
    }

    // ------------------------------------------------------------------
    // g.save  field=NAME, file=path/file.vti
    //      or field=NAME, file=path/file.csv  (writes CSV matrix)
    //      or field=NAME, file=path/file.txt
    // ------------------------------------------------------------------
    if (method == "save") {
        const std::string field = ctx.argVal(args, "field", "");
        const std::string file  = ctx.expand(ctx.argVal(args, "file", field + ".vti"));

        if (field.empty())
            throw std::runtime_error("line " + std::to_string(lineno)
                + ": g.save requires field=NAME");

        if (rank == 0) {
            ensureDir(file.substr(0, file.find_last_of("/\\")), rank);

            // Detect kind from field name
            const auto kind = kindOf(field, g);
            const bool field_ok = g.hasField(field);
            const bool flux_ok  = g.hasFlux(field);

            if (!field_ok && !flux_ok)
                throw std::runtime_error("line " + std::to_string(lineno)
                    + ": field/flux '" + field + "' does not exist in grid");

            const std::string ext = [&]{
                auto p = file.rfind('.');
                return (p == std::string::npos) ? "" : file.substr(p);
            }();

            if (ext == ".vti" || ext == ".vtk") {
                g.writeNamedVTI(field, kind, file);
            } else {
                // CSV / TXT
                g.writeNamedMatrix(field, kind, file);
            }
            std::cout << "  save " << field << " -> " << file << "\n";
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        return true;
    }

    // ------------------------------------------------------------------
    // g.set_diffusion  value=0.01
    // g.set_porosity   value=1.0
    // g.set_left_bc    value=1.0
    // ------------------------------------------------------------------
    if (method == "set_diffusion") {
        g.SetVal("diffusion", ctx.argDouble(args, "value", 0.0));
        return true;
    }
    if (method == "set_porosity") {
        g.SetVal("porosity", ctx.argDouble(args, "value", 1.0));
        return true;
    }
    if (method == "set_left_bc") {
        g.SetVal("c_left", ctx.argDouble(args, "value", 1.0));
        return true;
    }

    // ------------------------------------------------------------------
    // g.set_transport  diffusion=1e-4, porosity=0.3, c_left=1.0
    //                  (any subset)
    // ------------------------------------------------------------------
    if (method == "set_transport") {
        if (ctx.argHas(args, "diffusion") || ctx.argHas(args, "D")) {
            const double D = ctx.argDouble(args,
                                           ctx.argHas(args, "diffusion") ? "diffusion" : "D",
                                           0.0);
            g.SetVal("diffusion", D);
        }
        if (ctx.argHas(args, "porosity")) {
            g.SetVal("porosity", ctx.argDouble(args, "porosity", 1.0));
        }
        if (ctx.argHas(args, "c_left") || ctx.argHas(args, "left_bc")) {
            const double c = ctx.argDouble(args,
                                           ctx.argHas(args, "c_left") ? "c_left" : "left_bc",
                                           1.0);
            g.SetVal("c_left", c);
        }
        if (rank == 0)
            std::cout << "  set_transport applied\n";
        return true;
    }

    // ------------------------------------------------------------------
    // g.set_btc_locations  x=0.5,1.0,1.5,2.0,2.5
    // ------------------------------------------------------------------
    if (method == "set_btc_locations") {
        const std::string xstr = ctx.argVal(args, "x", "");
        if (xstr.empty())
            throw std::runtime_error("line " + std::to_string(lineno)
                                     + ": set_btc_locations requires x=val,val,...");

        std::vector<double> locs;
        std::stringstream ss(xstr);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            const auto first = tok.find_first_not_of(" \t");
            const auto last  = tok.find_last_not_of(" \t");
            if (first == std::string::npos) continue;
            tok = tok.substr(first, last - first + 1);
            if (!tok.empty()) locs.push_back(std::stod(tok));
        }
        g.setBTCLocations(locs);
        if (rank == 0) {
            std::cout << "  set_btc_locations: " << locs.size() << " locations: ";
            for (size_t k = 0; k < locs.size(); ++k)
                std::cout << locs[k] << (k + 1 < locs.size() ? ", " : "\n");
        }
        return true;
    }

    // ------------------------------------------------------------------
    // g.solve_transport  t_end=10.0, dt=0.01,
    //                    output_interval=10,
    //                    output_dir=$out/r$i,
    //                    output_base=transport_C,
    //                    realization=$i,
    //                    save_btc=$out/r$i/btc.csv,
    //                    accumulate=btc_mean
    // ------------------------------------------------------------------
    if (method == "solve_transport") {
        const double      t_end       = ctx.argDouble(args, "t_end", 1.0);
        const double      dt_user     = ctx.argDouble(args, "dt",          -1.0);
        const double      dt_factor   = ctx.argDouble(args, "dt_factor",    0.5);
        const double      dt_max      = ctx.argDouble(args, "dt_max",      -1.0);
        const int         out_interval= ctx.argInt   (args, "output_interval", 0);
        const std::string out_dir     = ctx.argVal   (args, "output_dir",  "");
        const std::string out_base    = ctx.argVal   (args, "output_base", "transport_C");
        const int         realization = ctx.argInt   (args, "realization", -1);
        const std::string savebtc     = ctx.argVal   (args, "save_btc",    "");
        const std::string accname     = ctx.argVal   (args, "accumulate",  "");

        // Compute CFL-stable dt from the actual qx field of this realization.
        // dt_cfl = dt_factor * dx / qx_max, matching run_fine_loop_collect.
        // If the user supplied an explicit dt, use it; otherwise use dt_cfl,
        // optionally capped by dt_max.
        const auto qx_minmax = g.fieldMinMax("qx", Grid2D::ArrayKind::Fx);
        const double qx_max  = qx_minmax.second;
        const double dt_cfl  = (qx_max > 0.0) ? dt_factor * g.dx() / qx_max
                                             : t_end;   // diffusion-only fallback

        double dt = (dt_user > 0.0) ? dt_user : dt_cfl;
        if (dt_max > 0.0 && dt > dt_max) dt = dt_max;

        if (!out_dir.empty()) ensureDir(out_dir, rank);

        TimeSeriesSet<double>  btc_data;
        TimeSeriesSet<double>* btc_ptr =
            (!savebtc.empty() || !accname.empty() || !g.getBTCLocations().empty())
                ? &btc_data : nullptr;

        if (rank == 0) {
            std::cout << "  solve_transport t_end=" << t_end
                      << " dt=" << dt
                      << " (cfl=" << dt_cfl << ", qx_max=" << qx_max << ")"
                      << " realization=" << realization << "\n";
        }

        const std::string ksp_pref = ctx.argVal(args, "ksp_prefix", "transport_");
        g.SolveTransport(t_end, dt,
                         ksp_pref.empty() ? nullptr : ksp_pref.c_str(),
                         out_interval, out_dir, out_base,
                         btc_ptr, realization);

        if (!savebtc.empty() && btc_ptr) {
            if (rank == 0) {
                ensureDir(savebtc.substr(0, savebtc.find_last_of("/\\")), rank);
                btc_data.write(savebtc);
                std::cout << "  saved per-realization BTC -> " << savebtc << "\n";
            }
            MPI_Barrier(PETSC_COMM_WORLD);
        }

        if (!accname.empty() && btc_ptr) {
            if (st.btc_accums.find(accname) == st.btc_accums.end()) {
                st.btc_accums[accname] = std::make_unique<BTCAccumulator>();
                if (rank == 0)
                    std::cout << "  [auto] created BTC accumulator '" << accname << "'\n";
            }
            st.btc(accname).add(btc_data);
        }

        return true;
    }
    throw std::runtime_error("line " + std::to_string(lineno)
        + ": unknown method '" + obj + "." + method + "'");
}

// -----------------------------------------------------------------------
// Execute a list of statements
// -----------------------------------------------------------------------
bool execStmts(const std::vector<S2Stmt>& stmts, RuntimeState& st)
{
    for (const auto& stmt : stmts) {
        switch (stmt.kind) {

        case S2Kind::Blank:
            break;

        case S2Kind::Set: {
            const std::string val = st.ctx.expand(stmt.val);
            st.ctx.set(stmt.var, val);
            if (st.rank == 0)
                std::cout << "set " << stmt.var << " = " << val << "\n";
            break;
        }

        case S2Kind::GridCreate: {
            const int    nx = st.ctx.argInt(stmt.args, "nx", 60);
            const int    ny = st.ctx.argInt(stmt.args, "ny", 20);
            const double Lx = [&]{ double v=3.0; for(auto& a:stmt.args) if(a.key=="Lx") v=std::stod(st.ctx.expand(a.value)); return v; }();
            const double Ly = [&]{ double v=1.0; for(auto& a:stmt.args) if(a.key=="Ly") v=std::stod(st.ctx.expand(a.value)); return v; }();
            st.grids[stmt.obj] = std::make_unique<Grid2D>(nx, ny, Lx, Ly);
            if (st.rank == 0)
                std::cout << "grid " << stmt.obj << " = "
                          << nx << "x" << ny << "  Lx=" << Lx << " Ly=" << Ly << "\n";
            break;
        }

        case S2Kind::AccumCreate: {
            const int bins = st.ctx.argInt(stmt.args, "bins", 20);
            st.accums[stmt.obj] = std::make_unique<CopulaAccumulator>(bins);
            if (st.rank == 0)
                std::cout << "accumulator " << stmt.obj << " bins=" << bins << "\n";
            break;
        }

        case S2Kind::MethodCall:
            if (!execMethod(stmt, st)) return false;
            break;

        case S2Kind::RepeatBegin: {
            if (st.rank == 0)
                std::cout << "\n--- repeat " << stmt.loop_var
                          << " = " << stmt.loop_start
                          << ":" << stmt.loop_end << " ---\n";

            for (int i = stmt.loop_start; i <= stmt.loop_end; ++i) {
                st.ctx.set(stmt.loop_var, std::to_string(i));

                if (st.rank == 0)
                    std::cout << "\n=== " << stmt.loop_var << " = " << i << " ===\n";

                if (!execStmts(stmt.body, st)) return false;
            }
            break;
        }

        case S2Kind::ForeachBegin: {
            // Expand $vars in the raw "lo:hi:n" string, then parse it now.
            const std::string rng_str = st.ctx.expand(stmt.foreach_rng);

            double lo = 0.0, hi = 0.0;
            int    n  = 1;
            {
                const auto c1 = rng_str.find(':');
                const auto c2 = (c1 == std::string::npos)
                                    ? std::string::npos
                                    : rng_str.find(':', c1 + 1);
                if (c1 == std::string::npos || c2 == std::string::npos)
                    throw std::runtime_error("line " + std::to_string(stmt.lineno)
                                             + ": 'foreach' expects lo:hi:n, got '"
                                             + rng_str + "'");
                try {
                    lo = std::stod(rng_str.substr(0, c1));
                    hi = std::stod(rng_str.substr(c1 + 1, c2 - c1 - 1));
                    n  = std::stoi(rng_str.substr(c2 + 1));
                } catch (...) {
                    throw std::runtime_error("line " + std::to_string(stmt.lineno)
                                             + ": 'foreach' could not parse '"
                                             + rng_str + "'");
                }
                if (n < 1)
                    throw std::runtime_error("line " + std::to_string(stmt.lineno)
                                             + ": 'foreach' n must be >= 1");
                if (lo <= 0.0 || hi <= 0.0)
                    throw std::runtime_error("line " + std::to_string(stmt.lineno)
                                             + ": 'foreach' requires positive lo and hi "
                                               "(geometric spacing)");
            }

            if (st.rank == 0)
                std::cout << "\n--- foreach " << stmt.loop_var
                          << " = " << lo << ":" << hi << ":" << n
                          << "  (geometric) ---\n";

            for (int k = 0; k < n; ++k) {
                const double exponent = (n == 1) ? 0.0
                                                 : static_cast<double>(k) / (n - 1);
                const double d = lo * std::pow(hi / lo, exponent);

                std::ostringstream oss;
                oss << std::setprecision(10) << d;
                st.ctx.set(stmt.loop_var, oss.str());
                st.ctx.set(stmt.loop_var + "_idx", std::to_string(k + 1));

                if (st.rank == 0)
                    std::cout << "\n=== " << stmt.loop_var << " = " << d
                              << "  (idx " << (k + 1) << "/" << n << ") ===\n";

                if (!execStmts(stmt.body, st)) return false;
            }
            break;
        }

        case S2Kind::BlockOpen:
        case S2Kind::BlockClose:
            throw std::runtime_error("execStmts: unexpected block delimiter at line "
                                     + std::to_string(stmt.lineno));
        }
    }
    return true;
}
} // anonymous namespace

// -----------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------
bool executeScript2(const S2Program& program, int rank)
{
    RuntimeState st;
    st.rank = rank;

    try {
        return execStmts(program, st);
    } catch (const std::exception& e) {
        if (rank == 0)
            std::cerr << "[fatal] Script2 execution error: " << e.what() << "\n";
        return false;
    }
}

bool runScript2File(const std::string& path, int rank)
{
    try {
        if (rank == 0)
            std::cout << "=== Script2: " << path << " ===\n\n";

        S2Program prog = parseScript2File(path);

        if (rank == 0)
            std::cout << "Parsed " << prog.size()
                      << " top-level statements\n\n";

        return executeScript2(prog, rank);

    } catch (const std::exception& e) {
        if (rank == 0)
            std::cerr << "[fatal] " << e.what() << "\n";
        return false;
    }
}

// -----------------------------------------------------------------------
// CopulaAccumulator::saveMean  (implementation here; declared in header)
// -----------------------------------------------------------------------
void CopulaAccumulator::saveMean(const std::string& path) const
{
    const auto m = mean();
    writeCopulaCSV(m, bins_, path);
    std::cout << "  Saved mean copula (" << bins_ << "x" << bins_
              << ", " << count_ << " realizations) -> " << path << "\n";
}

void BTCAccumulator::saveMean(const std::string& path,
                              double time_eps,
                              int    start_item) const
{
    if (count_ == 0 || sets_.empty())
        throw std::runtime_error("BTCAccumulator::saveMean: no data");

    TimeSeriesSet<double> mean =
        mean_ts_longest_cols<double>(sets_, start_item, time_eps);

    mean.write(path);   // or whatever your TimeSeriesSet write method is named
}

void CopulaAccumulator::saveMeanVTI(const std::string& path) const
{
    if (count_ == 0 || sum_.empty())
        throw std::runtime_error("CopulaAccumulator::saveMeanVTI: no data");

    const std::vector<double> m = mean();   // flat bins_ x bins_, row-major

    // Convert flat row-major vector to CMatrix.
    // Our convention: m[r * bins_ + c] where r = u1-bin, c = u2-bin.
    // write_matrix_as_vti_2d treats M[i][j] as the value at (i, j) in
    // the output grid (i along x, j along y), so we map r -> i, c -> j.
    CMatrix M(bins_, bins_);
    for (int i = 0; i < bins_; ++i) {
        for (int j = 0; j < bins_; ++j) {
            M[i][j] = m[i * bins_ + j];
        }
    }

    if (!write_matrix_as_vti_2d(M, path, /*array_name=*/"copula",
                                /*point_data=*/false))
        throw std::runtime_error("CopulaAccumulator::saveMeanVTI: "
                                 "write_matrix_as_vti_2d failed for '" + path + "'");
}

