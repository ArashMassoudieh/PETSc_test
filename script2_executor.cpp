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

// -----------------------------------------------------------------------
// Runtime state: named grids and accumulators
// -----------------------------------------------------------------------
namespace {

struct RuntimeState
{
    S2Context                                    ctx;
    std::map<std::string, std::unique_ptr<Grid2D>>          grids;
    std::map<std::string, std::unique_ptr<CopulaAccumulator>> accums;
    // Mean CDF accumulator: sum of per-realization CDFs (vector of (u,v) pairs)
    // keyed by a label the user assigns with g.extract_cdf output=NAME
    std::map<std::string, std::vector<TimeSeries<double>>>  cdf_sets;

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
// CopulaAccumulator::saveMean  (implementation here; declared in header)
// -----------------------------------------------------------------------
void CopulaAccumulator::saveMean(const std::string& path) const
{
    const auto m = mean();
    writeCopulaCSV(m, bins_, path);
    std::cout << "  Saved mean copula (" << bins_ << "x" << bins_
              << ", " << count_ << " realizations) -> " << path << "\n";
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
    const std::string& obj    = stmt.obj;
    const std::string& method = stmt.method;
    const auto&        args   = stmt.args;
    S2Context&         ctx    = st.ctx;
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

        throw std::runtime_error("line " + std::to_string(lineno)
            + ": unknown accumulator method '" + method + "'");
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

        // Use the representative delta = geometric mean of range
        const double delta_rep = std::sqrt(rng.lo * rng.hi);

        // Sample correlation vs delta (log-spaced) and build copula at delta_rep
        TimeSeries<double> corr_ts;   // correlation vs delta
        std::vector<double> copula_mat;

        for (int k = 0; k < rng.n; ++k) {
            const double exponent = static_cast<double>(k) / (rng.n - 1);
            const double delta    = rng.lo * std::pow(rng.hi / rng.lo, exponent);

            try {
                TimeSeries<double> pairs = g.sampleGaussianPerturbation(
                    field, kind, n_samples, delta, 0, dir);

                // Correlation of ranks at this delta
                corr_ts.append(delta, pairs.correlation_tc());

                // Bin into copula at representative delta
                if (std::abs(delta - delta_rep) <=
                    std::abs(delta - delta_rep) + 1e-12) {
                    // Pick the delta closest to delta_rep
                    // (simpler: just use the last computed for now,
                    //  and overwrite — caller can use save_corr to see all)
                    copula_mat = binnedCopula(pairs, bins);
                }
            } catch (...) {
                if (rank == 0)
                    std::cerr << "  [warn] extract_copula: failed at delta=" << delta << "\n";
            }
        }

        // Find closest delta index to geometric mean
        {
            double best = 1e30;
            for (int k = 0; k < rng.n; ++k) {
                const double exponent = static_cast<double>(k) / (rng.n - 1);
                const double delta    = rng.lo * std::pow(rng.hi / rng.lo, exponent);
                if (std::abs(std::log(delta) - std::log(delta_rep)) < best) {
                    best = std::abs(std::log(delta) - std::log(delta_rep));
                    // Re-sample at this best delta for the copula
                    try {
                        TimeSeries<double> pairs = g.sampleGaussianPerturbation(
                            field, kind, n_samples, delta, 0, dir);
                        copula_mat = binnedCopula(pairs, bins);
                    } catch (...) {}
                }
            }
        }

        if (rank == 0) {
            std::cout << "  extract_copula " << field << " dir=" << dirstr
                      << " range=[" << rng.lo << ":" << rng.hi << ":" << rng.n << "]"
                      << " bins=" << bins << "\n";

            // Save per-realization copula
            if (!savef.empty()) {
                ensureDir(savef.substr(0, savef.find_last_of("/\\")), rank);
                if (!copula_mat.empty())
                    writeCopulaCSV(copula_mat, bins, savef);
            }

            // Save correlation-vs-delta
            if (!savecorr.empty()) {
                ensureDir(savecorr.substr(0, savecorr.find_last_of("/\\")), rank);
                corr_ts.writefile(savecorr);
            }
        }

        // Accumulate into mean (all ranks do the accumulation so it's consistent,
        // but only rank 0 writes)
        if (!accname.empty() && !copula_mat.empty()) {
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
            // Build mean CDF on a uniform u grid
            const int nu = (int)cdfs[0].size();
            TimeSeries<double> mean_cdf;
            for (int k = 0; k < nu; ++k) {
                double u_sum = 0.0, v_sum = 0.0;
                for (const auto& ts : cdfs) {
                    u_sum += ts.getTime(k);
                    v_sum += ts.getValue(k);
                }
                mean_cdf.append(u_sum / cdfs.size(), v_sum / cdfs.size());
            }
            ensureDir(savef.substr(0, savef.find_last_of("/\\")), rank);
            mean_cdf.writefile(savef);
            std::cout << "  save_mean_cdf " << cdfname << " ("
                      << cdfs.size() << " realizations) -> " << savef << "\n";
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
                // Make loop variable available as $VAR
                st.ctx.set(stmt.loop_var, std::to_string(i));

                if (st.rank == 0)
                    std::cout << "\n=== " << stmt.loop_var << " = " << i << " ===\n";

                if (!execStmts(stmt.body, st)) return false;
            }
            break;
        }

        case S2Kind::BlockOpen:
        case S2Kind::BlockClose:
            // Should never appear in a flat list after parsing
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
