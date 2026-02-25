// sim_runner.cpp
#include "sim_runner.h"
#include "sim_helpers.h"

#include <mpi.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <sstream>
#include <cctype>
#include <cmath>
#include <limits>
#include <iomanip>
#include <filesystem>

#include "grid.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline double du_from_nu(int nu)
{
    return (nu > 0) ? (1.0 / double(nu)) : 1.0;
}

// ============================================================================
// Forward declarations for file-local helpers (static)
// ============================================================================

static bool build_mean_for_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_dir,
    const std::string& resume_run_dir,
    FineScaleOutputs& fine_outputs,
    TimeSeries<double>& invcdf_mean,
    int rank);

static bool run_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    double lc_mean,
    double lx_mean,
    double ly_mean,
    double nu_x_mean,
    double nu_y_mean,
    double dt_mean,
    const TimeSeries<double>& invcdf_mean,
    std::string& up_dir_out,
    std::string& up_btc_path_out,
    std::string& up_btc_deriv_path_out,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs_transport,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_pdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_cdf,
    int rank);

// ============================================================================
// prepare_run_dir_mpi
// ============================================================================

std::string prepare_run_dir_mpi(
    const std::string& output_dir,
    const std::string& resume_run_dir,
    const RunOptions& opts,
    int rank,
    const std::string& run_tag)
{
    std::string run_dir;

    if (rank == 0) {
        createDirectory(output_dir);

        // IMPORTANT:
        // In upscale-only mode, ALWAYS create a NEW output run_dir.
        // resume_run_dir is treated as INPUT SOURCE only.
        if (opts.upscale_only) {

            if (!opts.hardcoded_mean) {
                if (resume_run_dir.empty()) {
                    std::cerr << "ERROR: --upscale-only (without --hardcoded-mean) requires --run-dir=/path/to/input_run\n";
                    MPI_Abort(PETSC_COMM_WORLD, 1);
                }
                if (!dirExists(resume_run_dir)) {
                    std::cerr << "ERROR: input resume_run_dir does not exist: " << resume_run_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 2);
                }
            }

            std::string name = "run_" + makeTimestamp();
            if (!run_tag.empty()) name += "_" + run_tag;
            run_dir = joinPath(output_dir, name);
            createDirectory(run_dir);

            if (opts.hardcoded_mean) {
                std::cout << "Hardcoded-mean upscaled-only: created NEW run_dir: "
                          << run_dir << "\n";
            } else {
                std::cout << "Upscale-only: INPUT source = " << resume_run_dir << "\n";
                std::cout << "Upscale-only: OUTPUT run_dir = " << run_dir << "\n";
            }

        } else {
            std::string name = "run_" + makeTimestamp();
            if (!run_tag.empty()) name += "_" + run_tag;
            run_dir = joinPath(output_dir, name);
            createDirectory(run_dir);
            std::cout << "Run directory: " << run_dir << "\n";
        }
    }

    // broadcast run_dir
    int len = 0;
    if (rank == 0) len = (int)run_dir.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    run_dir.resize(len);
    MPI_Bcast(run_dir.data(), len, MPI_CHAR, 0, PETSC_COMM_WORLD);

    if (!run_dir.empty() && run_dir.back() != '/' && run_dir.back() != '\\') run_dir += "/";
    MPI_Barrier(PETSC_COMM_WORLD);

    return run_dir;
}

// ============================================================================
// small helpers
// ============================================================================

static inline void build_flat_invcdf(TimeSeries<double>& invcdf_mean, double qx_const)
{
    invcdf_mean.clear();
    invcdf_mean.append(0.0, qx_const);
    invcdf_mean.append(1.0, qx_const);
}

static inline std::string trim_copy(std::string s)
{
    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

static inline std::string to_lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (unsigned char)std::tolower(c); });
    return s;
}

static bool try_load_hardcoded_invcdf(TimeSeries<double>& invcdf_mean, const std::string& path)
{
    if (path.empty()) return false;
    if (!fileExists(path)) return false;

    TimeSeries<double> tmp;
    if (!read_inverse_cdf_any_format(path, tmp)) return false;
    if (tmp.size() < 2) return false;

    invcdf_mean = tmp;
    return true;
}

static inline std::string dirname_of(std::string path)
{
    const auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return std::string();
    return path.substr(0, pos);
}

static bool build_mean_qx_inverse_cdf_from_multi(
    const std::string& in_multi_path,
    const std::string& out_mean_path
){
    std::ifstream f(in_multi_path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    const char delim = detect_delim(header);

    std::vector<double> u_out;
    std::vector<double> v_out;
    u_out.reserve(4096);
    v_out.reserve(4096);

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        const auto tok = split_line_delim(line, delim);
        if (tok.size() < 2) continue;

        const size_t npairs = tok.size() / 2;
        if (npairs < 1) continue;

        double t_ref = 0.0;
        double sum_q = 0.0;
        int    cnt_q = 0;

        for (size_t k = 0; k < npairs; ++k) {
            double t = 0.0, q = 0.0;
            if (!try_parse_double(tok[2*k + 0], t)) continue;
            if (!try_parse_double(tok[2*k + 1], q)) continue;

            if (k == 0) t_ref = t;
            if (std::abs(t - t_ref) > 1e-10) continue;

            if (is_finite_number(q)) {
                sum_q += q;
                cnt_q++;
            }
        }

        if (cnt_q <= 0) continue;

        u_out.push_back(t_ref);
        v_out.push_back(sum_q / (double)cnt_q);
    }

    if (u_out.size() < 2) return false;

    std::ofstream o(out_mean_path);
    if (!o) return false;

    o << "u,v\n";
    o << std::setprecision(15);
    for (size_t i = 0; i < u_out.size(); ++i) {
        o << u_out[i] << "," << v_out[i] << "\n";
    }
    return true;
}

// ============================================================================
// Robust realization reject helpers (delete folder + log)
// ============================================================================

static bool remove_dir_recursive(const std::string& p)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    if (p.empty()) return true;
    if (!fs::exists(fs::path(p), ec)) return true;
    fs::remove_all(fs::path(p), ec);
    return !ec;
}

static void append_reject_log(const std::string& run_dir,
                              int realization,
                              int attempt,
                              const std::string& stage,
                              const std::string& what)
{
    const std::string path = joinPath(run_dir, "rejected_realizations.csv");
    const bool exists = fileExists(path);

    std::ofstream f(path, std::ios::app);
    if (!f) return;

    if (!exists) f << "realization,attempt,stage,message\n";

    auto esc = [](std::string s){
        std::replace(s.begin(), s.end(), '\n', ' ');
        return s;
    };

    f << realization << "," << attempt << "," << esc(stage) << "," << esc(what) << "\n";
}

// ============================================================================
// PT mean-arrival helpers
// ============================================================================

static inline double trapz_integral(const TimeSeries<double>& ts)
{
    const int n = (int)ts.size();
    if (n < 2) return 0.0;

    double area = 0.0;
    for (int i = 1; i < n; ++i) {
        const double t0 = ts.getTime(i - 1);
        const double t1 = ts.getTime(i);
        const double y0 = ts.getValue(i - 1);
        const double y1 = ts.getValue(i);
        const double dt = (t1 - t0);
        if (!std::isfinite(dt) || dt <= 0.0) continue;
        if (!std::isfinite(y0) || !std::isfinite(y1)) continue;
        area += 0.5 * (y0 + y1) * dt;
    }
    return area;
}

static inline double trapz_integral_t_times_y(const TimeSeries<double>& ts)
{
    const int n = (int)ts.size();
    if (n < 2) return 0.0;

    double area = 0.0;
    for (int i = 1; i < n; ++i) {
        const double t0 = ts.getTime(i - 1);
        const double t1 = ts.getTime(i);
        const double y0 = ts.getValue(i - 1);
        const double y1 = ts.getValue(i);
        const double dt = (t1 - t0);
        if (!std::isfinite(dt) || dt <= 0.0) continue;
        if (!std::isfinite(y0) || !std::isfinite(y1)) continue;
        area += 0.5 * (t0 * y0 + t1 * y1) * dt;
    }
    return area;
}

static inline double mean_time_from_pdf(const TimeSeries<double>& pdf_ts)
{
    const double den = trapz_integral(pdf_ts);
    if (!(den > 0.0) || !std::isfinite(den)) return std::numeric_limits<double>::quiet_NaN();
    const double num = trapz_integral_t_times_y(pdf_ts);
    return num / den;
}

static bool read_pt_mean_csv(const std::string& path, std::vector<double>& x, std::vector<double>& mu)
{
    x.clear(); mu.clear();

    std::ifstream f(path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    const char delim = detect_delim(header);

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        const auto tok = split_line_delim(line, delim);
        if (tok.size() < 2) continue;

        double xv = 0.0, muv = 0.0;
        if (!try_parse_double(tok[0], xv)) continue;
        if (!try_parse_double(tok[1], muv)) continue;

        x.push_back(xv);
        mu.push_back(muv);
    }

    return !x.empty();
}

// ============================================================================
// Fine-scale helpers
// ============================================================================

static void generateKField(Grid2D& g, const SimParams& P, unsigned long seed)
{
    if (P.CorrelationModel == SimParams::correlationmode::gaussian) {
        g.generateGaussianRandomField("K_normal_score",
                                      P.correlation_ls_x,
                                      P.correlation_ls_y,
                                      0.0, 1.0, seed);
    }
    else {
        g.makeGaussianFieldSGS("K_normal_score", P.correlation_ls_x, P.correlation_ls_y, 10, seed);
    }
}

// ============================================================================
// Wiener diffusion runner (kept)
// ============================================================================

int runDiffusionSimulation(const RunOptions &opts, int realization, const std::string &output_dir)
{
    std::string m = to_lower_copy(trim_copy(opts.wiener_mode));
    WienerParams wp;
    wp.dt   = opts.wiener_dt;
    wp.D    = opts.wiener_Dx;
    wp.rx   = opts.wiener_rx;
    wp.ry   = opts.wiener_ry;
    wp.seed = opts.wiener_seed + (unsigned long)realization;
    if (m == "1dx")      wp.mode = WienerMode::W1D_X;
    else if (m == "1dy") wp.mode = WienerMode::W1D_Y;
    else                 wp.mode = WienerMode::W2D;

    std::cout << "Running diffusion simulation\n"
              << "  Mode: " << opts.wiener_mode << "\n"
              << "  D = " << wp.D << ", dt = " << wp.dt << "\n"
              << "  Ellipse: rx = " << wp.rx << ", ry = " << wp.ry << "\n"
              << "  Particles: " << realization << "\n";

    PathwaySet pathwaySet;
    pathwaySet.InitializeAtOrigin(realization);

    TimeSeries<double> times = pathwaySet.trackDiffusion(wp.dt, wp.rx, wp.ry, wp.D);
    times.writefile(joinPath(output_dir, "diffusiontimes.csv"));

    double mean     = times.mean();
    double stdev    = times.stddev();
    double meanlog  = times.mean_log();
    double stdevlog = times.log().stddev();
    auto [tmin, tmax] = pathwaySet.travelTimeRange();

    TimeSeries<double> fptDist = times.distribution(100, 0);
    fptDist.writefile(joinPath(output_dir, "FirstPassageTimeDistribution.csv"));

    {
        std::ofstream f(joinPath(output_dir, "diffusion_stats.txt"));
        if (f.is_open()) {
            f << "# First Passage Time Statistics\n"
              << "mode=" << opts.wiener_mode << "\n"
              << "D=" << wp.D << "\n"
              << "dt=" << wp.dt << "\n"
              << "rx=" << wp.rx << "\n"
              << "ry=" << wp.ry << "\n"
              << "seed=" << wp.seed << "\n"
              << "num_particles=" << realization << "\n"
              << "mean=" << mean << "\n"
              << "stdev=" << stdev << "\n"
              << "mean_log=" << meanlog << "\n"
              << "stdev_log=" << stdevlog << "\n"
              << "min=" << tmin << "\n"
              << "max=" << tmax << "\n";
        }
    }

    return 0;
}

// ============================================================================
// Mean correlations + stats helpers
// ============================================================================

static void writeMeanCorrelations(
    const std::string& run_dir,
    const SimParams& P,
    FineScaleOutputs& outputs)
{
    outputs.advective_correlations.write(joinPath(run_dir, "advective_correlations.txt"));
    outputs.velocity_x_correlations.write(joinPath(run_dir, "diffusion_x_correlations.txt"));
    outputs.velocity_y_correlations.write(joinPath(run_dir, "diffusion_y_correlations.txt"));
    outputs.K_x_correlations.write(joinPath(run_dir, "K_x_correlations.txt"));
    outputs.K_y_correlations.write(joinPath(run_dir, "K_y_correlations.txt"));

    // IMPORTANT: union-t mean (NEW)
    auto mean_corr_a    = mean_ts_union_t(outputs.advective_correlations);
    auto mean_corr_x    = mean_ts_union_t(outputs.velocity_x_correlations);
    auto mean_corr_y    = mean_ts_union_t(outputs.velocity_y_correlations);
    auto mean_K_corr_x  = mean_ts_union_t(outputs.K_x_correlations);
    auto mean_K_corr_y  = mean_ts_union_t(outputs.K_y_correlations);

    mean_corr_a.writefile(joinPath(run_dir, "advective_correlations_mean.txt"));
    mean_corr_x.writefile(joinPath(run_dir, "diffusion_x_correlations_mean.txt"));
    mean_corr_y.writefile(joinPath(run_dir, "diffusion_y_correlations_mean.txt"));
    mean_K_corr_x.writefile(joinPath(run_dir, "K_x_correlations_mean.txt"));
    mean_K_corr_y.writefile(joinPath(run_dir, "K_y_correlations_mean.txt"));

    double l_a = mean_corr_a.fitExponentialDecay();
    TimeSeries<double> fitted_a;
    for (size_t i = 0; i < mean_corr_a.size(); ++i) {
        double r = mean_corr_a.getTime(i);
        fitted_a.append(r, std::exp(-r / l_a));
    }
    fitted_a.writefile(joinPath(run_dir, "advective_correlations_mean_fitted.txt"));

    auto fitAndWrite = [&](const TimeSeries<double>& mean_corr,
                           const std::string& filename,
                           double& nu_out) {
        double ell = 0, nu = 0.0;
        if (P.CorrelationModel == SimParams::correlationmode::exponentialfit)
            ell = mean_corr.fitExponentialDecay();
        else if (P.CorrelationModel == SimParams::correlationmode::derivative)
            ell = mean_corr.getTime(0) / (1.0 - mean_corr.getValue(0));
        else if (P.CorrelationModel == SimParams::correlationmode::gaussian)
            ell = mean_corr.fitGaussianDecay();
        else {
            auto pr = mean_corr.fitMaternDecay();
            nu = pr.first; ell = pr.second;
        }

        TimeSeries<double> fitted;
        for (size_t i = 0; i < mean_corr.size(); ++i) {
            double r = mean_corr.getTime(i);
            double rho =
                (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                    std::exp(-(r / ell) * (r / ell)) :
                (P.CorrelationModel == SimParams::correlationmode::matern) ?
                    matern(r, nu, ell) :
                    std::exp(-r / ell);
            fitted.append(r, rho);
        }
        fitted.writefile(joinPath(run_dir, filename));
        nu_out = nu;
        return ell;
    };

    fitAndWrite(mean_corr_x, "diffusion_x_correlations_mean_fitted.txt", outputs.nu_x_mean);
    fitAndWrite(mean_corr_y, "diffusion_y_correlations_mean_fitted.txt", outputs.nu_y_mean);

    auto fitKAndWrite = [&](const TimeSeries<double>& mean_K_corr, const std::string& filename) {
        double ell_K =
            (P.CorrelationModel == SimParams::correlationmode::gaussian)       ? mean_K_corr.fitGaussianDecay() :
            (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? mean_K_corr.fitExponentialDecay() :
            (P.CorrelationModel == SimParams::correlationmode::derivative)      ? mean_K_corr.getTime(0) / (1.0 - mean_K_corr.getValue(0)) :
                                                                                mean_K_corr.fitMaternDecay().second;

        TimeSeries<double> fitted_K;
        for (size_t i = 0; i < mean_K_corr.size(); ++i) {
            double r = mean_K_corr.getTime(i);
            double rho =
                (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                    std::exp(-(r / ell_K) * (r / ell_K)) :
                    std::exp(-r / ell_K);
            fitted_K.append(r, rho);
        }
        fitted_K.writefile(joinPath(run_dir, filename));
        return ell_K;
    };

    outputs.lambda_K_x_mean = fitKAndWrite(mean_K_corr_x, "K_x_correlations_mean_fitted.txt");
    outputs.lambda_K_y_mean = fitKAndWrite(mean_K_corr_y, "K_y_correlations_mean_fitted.txt");
}

static void computeFinalMeans(FineScaleOutputs& outputs)
{
    outputs.lc_mean          = mean_of(outputs.lc_all);
    outputs.velocity_lx_mean = mean_of(outputs.velocity_lx_all);
    outputs.velocity_ly_mean = mean_of(outputs.velocity_ly_all);
    outputs.dt_mean          = mean_of(outputs.dt_all);
}

// ============================================================================
// Fine loop (RETRY/REJECT) — collects transport derivative BTCs + PT pdf/cdf stacks
// ============================================================================

static bool run_fine_loop_collect(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    FineScaleOutputs& outputs,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;
    const double du = du_from_nu(nu);

    const std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,lambda_K_x,lambda_K_y,dt_opt\n";
    }

    // Ensure stacks exist
    outputs.BTCs.resize(P.xLocations.size());
    outputs.PT_pdfs.resize(P.xLocations.size());
    outputs.PT_cdfs.resize(P.xLocations.size());

    const int nReal = P.nReal_default;

    // cap to avoid infinite loop if PETSc is fundamentally failing
    const int MAX_RETRY_PER_REAL = 50;

    int success = 0;
    while (success < nReal) {

        const int r = success + 1; // realization number to produce
        int attempt = 0;

        bool done_this_real = false;

        while (!done_this_real) {
            attempt++;

            const std::string rlab = makeRealLabel(r);
            std::string fine_dir = joinPath(run_dir, makeFineFolder(r));

            if (rank == 0) {
                if (dirExists(fine_dir)) remove_dir_recursive(fine_dir);
                createDirectory(fine_dir);
            }
            MPI_Barrier(PETSC_COMM_WORLD);

            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

            const std::string pfx = rlab + "_";

            // ---- locals (do NOT pollute outputs unless we succeed) ----
            TimeSeries<double> K_corr_x, K_corr_y;
            TimeSeries<double> vel_corr_x, vel_corr_y;
            TimeSeries<double> qx_corr_adv;

            double lambda_K_x_emp = std::numeric_limits<double>::quiet_NaN();
            double lambda_K_y_emp = std::numeric_limits<double>::quiet_NaN();
            double lambda_x_emp   = std::numeric_limits<double>::quiet_NaN();
            double lambda_y_emp   = std::numeric_limits<double>::quiet_NaN();
            double lc_emp         = std::numeric_limits<double>::quiet_NaN();
            double dt_optimal     = std::numeric_limits<double>::quiet_NaN();

            TimeSeries<double> qx_inverse_cdf;
            TimeSeries<double> qx_pdf;

            TimeSeriesSet<double> BTCs_FineScaled;
            TimeSeriesSet<double> PT_pdf_local;
            TimeSeriesSet<double> PT_cdf_local;

            try {
                Grid2D g(nx, ny, Lx, Ly);

                // IMPORTANT: change seed on retries so you don't repeat the same failing field
                unsigned long seed =
                    P.run_seed
                    + 1000UL * (unsigned long)r
                    + 1000000UL * (unsigned long)attempt
                    + (unsigned long)rank;

                generateKField(g, P, seed);

                // ---- K correlation sampling (local) ----
                {
                    const int num_deltas = 30;
                    const int num_samples_per_delta = 10000;

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_x_range.first *
                                       std::pow(P.correlation_x_range.second / P.correlation_x_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "K_normal_score", Grid2D::ArrayKind::Cell,
                                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                            K_corr_x.append(delta, samples.correlation_tc());
                        } catch (...) {}
                    }
                    K_corr_x.writefile(joinPath(fine_dir, pfx + "K_correlation_x.txt"));

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_y_range.first *
                                       std::pow(P.correlation_y_range.second / P.correlation_y_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "K_normal_score", Grid2D::ArrayKind::Cell,
                                num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                            K_corr_y.append(delta, samples.correlation_tc());
                        } catch (...) {}
                    }
                    K_corr_y.writefile(joinPath(fine_dir, pfx + "K_correlation_y.txt"));

                    lambda_K_x_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? K_corr_x.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? K_corr_x.getTime(0)/(1.0-K_corr_x.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? K_corr_x.fitGaussianDecay() :
                                                                                             K_corr_x.fitMaternDecay().second;

                    lambda_K_y_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? K_corr_y.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? K_corr_y.getTime(0)/(1.0-K_corr_y.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? K_corr_y.fitGaussianDecay() :
                                                                                             K_corr_y.fitMaternDecay().second;
                }

                g.normalizeField("K_normal_score", Grid2D::ArrayKind::Cell);
                g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
                g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.vti"));

                g.createExponentialField("K_normal_score", P.stdev, P.g_mean, "K");
                g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));

                // ---- PETSc solve: can throw (KSP did not converge) ----
                g.DarcySolve(Lx, 0, "K", "K");

                g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.vti"));
                g.writeNamedVTI("qx",   Grid2D::ArrayKind::Fx,   joinPath(fine_dir, pfx + "qx.vti"));
                g.writeNamedVTI("qy",   Grid2D::ArrayKind::Fy,   joinPath(fine_dir, pfx + "qy.vti"));

                // ---- velocity normal score + correlations (local) ----
                {
                    TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
                    TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
                    g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);

                    const int num_deltas = 30;
                    const int num_samples_per_delta = 10000;

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_x_range.first *
                                       std::pow(P.correlation_x_range.second / P.correlation_x_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "qx_normal_score", Grid2D::ArrayKind::Fx,
                                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                            vel_corr_x.append(delta, samples.correlation_tc());
                        } catch (...) {}
                    }
                    vel_corr_x.writefile(joinPath(fine_dir, pfx + "velocity_correlation_x.txt"));

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_y_range.first *
                                       std::pow(P.correlation_y_range.second / P.correlation_y_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "qx_normal_score", Grid2D::ArrayKind::Fx,
                                num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                            vel_corr_y.append(delta, samples.correlation_tc());
                        } catch (...) {}
                    }
                    vel_corr_y.writefile(joinPath(fine_dir, pfx + "velocity_correlation_y.txt"));

                    lambda_x_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? vel_corr_x.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? vel_corr_x.getTime(0)/(1.0-vel_corr_x.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? vel_corr_x.fitGaussianDecay() :
                                                                                             vel_corr_x.fitMaternDecay().second;

                    lambda_y_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? vel_corr_y.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? vel_corr_y.getTime(0)/(1.0-vel_corr_y.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? vel_corr_y.fitGaussianDecay() :
                                                                                             vel_corr_y.fitMaternDecay().second;
                }

                qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100, 1e-6);
                qx_pdf         = g.extractFieldPDF("qx", Grid2D::ArrayKind::Fx, 50,  1e-6, true);

                qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
                qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));
                qx_pdf.writefile(joinPath(fine_dir, pfx + "qx_pdf.txt"));

                dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;

                // ---- Transport (can also throw; reject if it does) ----
                g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
                g.SetVal("diffusion", P.Diffusion_coefficient);
                g.SetVal("porosity", 1.0);
                g.SetVal("c_left", 1.0);

                g.setBTCLocations(P.xLocations);

                if (opts.solve_fine_scale_transport) {
                    g.SolveTransport(P.t_end_pdf, std::min(dt_optimal, 0.5 / 10.0),
                                     "transport_", 500, fine_dir, "C", &BTCs_FineScaled, r);

                    BTCs_FineScaled.write(joinPath(fine_dir, pfx + "BTC_FineScaled.csv"));
                    BTCs_FineScaled.derivative().write(joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv"));
                }

                // ---- Advective correlation + PT locals ----
                {
                    PathwaySet pathways;
                    pathways.Initialize(10000, PathwaySet::Weighting::FluxWeighted, &g);
                    pathways.trackAllPathwaysWithDiffusion(&g, 0.01, g.getDiffusion());

                    const std::vector<double>& btc_locations = g.getBTCLocations();
                    if (!btc_locations.empty()) {
                        TimeSeriesSet<double> btc_pdf = pathways.getBreakthroughCurve(
                            btc_locations, true, PathwaySet::BTCType::PDF, 200);
                        btc_pdf.write(joinPath(fine_dir, pfx + "btc_pdf.csv"));

                        TimeSeriesSet<double> btc_cdf = pathways.getBreakthroughCurve(
                            btc_locations, true, PathwaySet::BTCType::CDF);
                        btc_cdf.write(joinPath(fine_dir, pfx + "btc_cdf.csv"));

                        if (opts.perform_particle_tracking) {
                            PT_pdf_local = btc_pdf;
                            PT_cdf_local = btc_cdf;

                            const std::string mean_csv = joinPath(fine_dir, pfx + "PT_mean_arrival_time.csv");
                            std::ofstream mf(mean_csv);
                            mf << "x,mean_arrival_time\n";
                            mf << std::setprecision(15);
                            for (int i = 0; i < (int)btc_locations.size() && i < (int)btc_pdf.size(); ++i) {
                                const double mu = mean_time_from_pdf(btc_pdf[i]);
                                mf << btc_locations[i] << "," << mu << "\n";
                            }
                        }
                    }

                    double Delta_x_min = 0.001, Delta_x_max = 0.5;
                    int num_Delta_x = 30;
                    int num_samples_per_Delta_x = 10000;

                    for (int i = 0; i < num_Delta_x; ++i) {
                        double exponent = static_cast<double>(i) / (num_Delta_x - 1);
                        double Delta_x = Delta_x_min * std::pow(Delta_x_max / Delta_x_min, exponent);
                        try {
                            PathwaySet particle_pairs = pathways.sampleParticlePairs(Delta_x, num_samples_per_Delta_x);
                            double correlation = particle_pairs.calculateCorrelation(0, 1, "qx");
                            qx_corr_adv.append(Delta_x, correlation);
                        } catch (...) {}
                    }

                    qx_corr_adv.writefile(joinPath(fine_dir, pfx + "qx_correlation_vs_distance.txt"));
                    lc_emp = qx_corr_adv.fitExponentialDecay();

                    pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
                    pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"), 1000);
                }

                // =========================
                // SUCCESS -> COMMIT to outputs
                // =========================
                outputs.K_x_correlations.append(K_corr_x);
                outputs.K_y_correlations.append(K_corr_y);

                outputs.velocity_x_correlations.append(vel_corr_x, "Realization" + aquiutils::numbertostring(r));
                outputs.velocity_y_correlations.append(vel_corr_y, "Realization" + aquiutils::numbertostring(r));
                outputs.advective_correlations.append(qx_corr_adv, "Realization" + aquiutils::numbertostring(r));

                outputs.velocity_lx_all.push_back(lambda_x_emp);
                outputs.velocity_ly_all.push_back(lambda_y_emp);
                outputs.K_lx_all.push_back(lambda_K_x_emp);
                outputs.K_ly_all.push_back(lambda_K_y_emp);
                outputs.dt_all.push_back(dt_optimal);
                outputs.lc_all.push_back(lc_emp);

                if (opts.solve_fine_scale_transport) {
                    for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_FineScaled.size(); ++i) {
                        BTCs_FineScaled.setSeriesName(i, fmt_x(P.xLocations[i]));
                        outputs.BTCs[i].append(BTCs_FineScaled[i].derivative(), pfx + fmt_x(P.xLocations[i]));
                    }
                }

                if (opts.perform_particle_tracking) {
                    const int nmin = std::min((int)P.xLocations.size(), (int)PT_pdf_local.size());
                    for (int i = 0; i < nmin; ++i)
                        outputs.PT_pdfs[i].append(PT_pdf_local[i], pfx + fmt_x(P.xLocations[i]));

                    const int nmin2 = std::min((int)P.xLocations.size(), (int)PT_cdf_local.size());
                    for (int i = 0; i < nmin2; ++i)
                        outputs.PT_cdfs[i].append(PT_cdf_local[i], pfx + fmt_x(P.xLocations[i]));
                }

                if (rank == 0) {
                    outputs.inverse_qx_cdfs.append(qx_inverse_cdf, "qx_inverse_cdf" + aquiutils::numbertostring(r));
                    outputs.qx_pdfs.append(qx_pdf, "qx_pdf" + aquiutils::numbertostring(r));

                    std::ofstream f(stats_csv, std::ios::app);
                    f << r << "," << lc_emp << "," << lambda_x_emp << "," << lambda_y_emp << ","
                      << lambda_K_x_emp << "," << lambda_K_y_emp << "," << dt_optimal << "\n";
                }

                done_this_real = true;

            } catch (const std::exception& e) {

                if (rank == 0) {
                    std::cerr << "Rejected realization " << r
                              << " (attempt " << attempt << "): " << e.what() << "\n";
                    append_reject_log(run_dir, r, attempt, "exception", e.what());
                    remove_dir_recursive(fine_dir);
                }

                MPI_Barrier(PETSC_COMM_WORLD);

                if (attempt >= MAX_RETRY_PER_REAL) {
                    if (rank == 0) {
                        std::cerr << "ERROR: realization " << r << " failed "
                                  << attempt << " times. Aborting.\n";
                    }
                    MPI_Abort(PETSC_COMM_WORLD, 123);
                    return false;
                }

                continue;

            } catch (...) {

                if (rank == 0) {
                    std::cerr << "Rejected realization " << r
                              << " (attempt " << attempt << "): unknown exception\n";
                    append_reject_log(run_dir, r, attempt, "unknown", "unknown exception");
                    remove_dir_recursive(fine_dir);
                }

                MPI_Barrier(PETSC_COMM_WORLD);

                if (attempt >= MAX_RETRY_PER_REAL) {
                    MPI_Abort(PETSC_COMM_WORLD, 124);
                    return false;
                }

                continue;
            }

            MPI_Barrier(PETSC_COMM_WORLD);
        }

        // only count after success
        success++;
    }

    // Post-processing
    if (rank == 0) {
        outputs.inverse_qx_cdfs.write(joinPath(run_dir, "qx_inverse_cdfs.txt"));
        outputs.qx_pdfs.write(joinPath(run_dir, "qx_pdfs.txt"));

        // IMPORTANT: union-t mean (NEW)
        mean_ts_union_t(outputs.qx_pdfs).writefile(joinPath(run_dir, "qx_mean_pdf.txt"));

        writeMeanCorrelations(run_dir, P, outputs);
        computeFinalMeans(outputs);

        // mean PT arrival time over realizations
        if (opts.perform_particle_tracking) {
            std::vector<double> sum_mu(P.xLocations.size(), 0.0);
            std::vector<int>    cnt_mu(P.xLocations.size(), 0);

            for (int r = 1; r <= nReal; ++r) {
                const std::string rlab = makeRealLabel(r);
                const std::string pfx  = rlab + "_";
                const std::string fine_dir = joinPath(run_dir, makeFineFolder(r));

                const std::string mean_csv = joinPath(fine_dir, pfx + "PT_mean_arrival_time.csv");
                if (!fileExists(mean_csv)) continue;

                std::vector<double> x, mu;
                if (!read_pt_mean_csv(mean_csv, x, mu)) continue;

                const int nmin = std::min((int)P.xLocations.size(), (int)mu.size());
                for (int i = 0; i < nmin; ++i) {
                    if (std::isfinite(mu[i])) {
                        sum_mu[i] += mu[i];
                        cnt_mu[i] += 1;
                    }
                }
            }

            const std::string out_mean_csv = joinPath(run_dir, "PT_mean_arrival_time_mean.csv");
            std::ofstream outm(out_mean_csv);
            outm << "x,mean_arrival_time\n";
            outm << std::setprecision(15);

            for (int i = 0; i < (int)P.xLocations.size(); ++i) {
                const double m = (cnt_mu[i] > 0) ? (sum_mu[i] / (double)cnt_mu[i])
                                                 : std::numeric_limits<double>::quiet_NaN();
                outm << P.xLocations[i] << "," << m << "\n";
            }
        }
    }

    return true;
}

// ============================================================================
// Mean-building helpers
// ============================================================================

static bool loadHardcodedInverseCDF(
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_cdf_copy_path,
    TimeSeries<double>& invcdf_mean)
{
    bool loaded = false;

    loaded = try_load_hardcoded_invcdf(invcdf_mean, opts.hardcoded_qx_cdf_path);

    if (!loaded && !opts.hardcoded_qx_cdf_path.empty()) {
        const std::string src_dir = dirname_of(opts.hardcoded_qx_cdf_path);
        if (!src_dir.empty()) {
            const std::string meanf = joinPath(src_dir, "mean_qx_inverse_cdf.txt");
            const std::string multi = joinPath(src_dir, "qx_inverse_cdfs.txt");

            if (!fileExists(meanf) && fileExists(multi)) {
                if (build_mean_qx_inverse_cdf_from_multi(multi, meanf)) {
                    loaded = try_load_hardcoded_invcdf(invcdf_mean, meanf);
                }
            }
            if (!loaded && fileExists(meanf)) {
                loaded = try_load_hardcoded_invcdf(invcdf_mean, meanf);
            }
        }
    }

    if (!loaded) {
        build_flat_invcdf(invcdf_mean, H.qx_const);
    }

    invcdf_mean.writefile(run_cdf_copy_path);
    return loaded;
}

static void writeHardcodedMeanParams(
    const std::string& run_dir,
    const FineScaleOutputs& outputs,
    const RunOptions& opts,
    const std::string& run_cdf_copy_path,
    bool cdf_loaded,
    double qx_const)
{
    std::ofstream f(joinPath(run_dir, "mean_params_used.txt"));
    f << "lc_mean=" << outputs.lc_mean << "\n";
    f << "lambda_x_mean=" << outputs.velocity_lx_mean << "\n";
    f << "lambda_y_mean=" << outputs.velocity_ly_mean << "\n";
    f << "lambda_K_x_mean=" << outputs.lambda_K_x_mean << "\n";
    f << "lambda_K_y_mean=" << outputs.lambda_K_y_mean << "\n";
    f << "dt_mean=" << outputs.dt_mean << "\n";
    f << "nu_x_mean=" << outputs.nu_x_mean << "\n";
    f << "nu_y_mean=" << outputs.nu_y_mean << "\n";
    f << "mean_qx_inverse_cdf_written_to=" << run_cdf_copy_path << "\n";
    f << "qx_cdf_preferred_path=" << opts.hardcoded_qx_cdf_path << "\n";
    if (!cdf_loaded) f << "qx_const=" << qx_const << "\n";
}

static void writeComputedMeanParams(
    const std::string& run_dir,
    const SimParams& P,
    const FineScaleOutputs& outputs)
{
    std::ofstream f(joinPath(run_dir, "mean_params.txt"));
    f << "nReal=" << (int)outputs.lc_all.size() << "\n";
    f << "lc_mean=" << outputs.lc_mean << "\n";
    f << "lambda_x_mean=" << outputs.velocity_lx_mean << "\n";
    f << "lambda_y_mean=" << outputs.velocity_ly_mean << "\n";
    f << "lambda_K_x_mean=" << outputs.lambda_K_x_mean << "\n";
    f << "lambda_K_y_mean=" << outputs.lambda_K_y_mean << "\n";
    f << "dt_mean=" << outputs.dt_mean << "\n";
    f << "matern_nu_x_mean=" << outputs.nu_x_mean << "\n";
    f << "matern_nu_y_mean=" << outputs.nu_y_mean << "\n";
    f << "du=" << du_from_nu(P.nu) << "\n";
}

static bool loadMeanParamsForUpscaleOnly_FromInputSource(
    const std::string& input_dir,
    FineScaleOutputs& outputs,
    TimeSeries<double>& invcdf_mean)
{
    const std::string mean_params_path = joinPath(input_dir, "mean_params.txt");
    const std::string mean_cdf_path    = joinPath(input_dir, "mean_qx_inverse_cdf.txt");

    bool ok_params = fileExists(mean_params_path) &&
                     read_mean_params_txt(mean_params_path,
                                          outputs.lc_mean,
                                          outputs.velocity_lx_mean,
                                          outputs.velocity_ly_mean,
                                          outputs.dt_mean,
                                          outputs.nu_x_mean,
                                          outputs.nu_y_mean);

    TimeSeries<double> mean_qx_cdf;
    bool ok_cdf = fileExists(mean_cdf_path) &&
                  read_inverse_cdf_any_format(mean_cdf_path, mean_qx_cdf);

    if (!ok_params || !ok_cdf) {
        std::cerr << "ERROR: --upscale-only requires existing mean_params.txt and mean_qx_inverse_cdf.txt in INPUT folder:\n"
                  << "  " << input_dir << "\n"
                  << "If you don't have them, run with --hardcoded-mean.\n";
        return false;
    }

    invcdf_mean = mean_qx_cdf;
    return true;
}

static bool build_mean_for_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_dir,
    const std::string& resume_run_dir,
    FineScaleOutputs& fine_outputs,
    TimeSeries<double>& invcdf_mean,
    int rank)
{
    if (rank != 0) return true; // rank0 builds; broadcast later

    const std::string run_cdf_copy_path = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

    // Case 1: hardcoded mean
    if (opts.hardcoded_mean) {
        fine_outputs.lc_mean          = H.lc_mean;
        fine_outputs.velocity_lx_mean = H.lx_mean;
        fine_outputs.velocity_ly_mean = H.ly_mean;
        fine_outputs.dt_mean          = H.dt_mean;
        fine_outputs.nu_x_mean        = H.nu_x;
        fine_outputs.nu_y_mean        = H.nu_y;

        bool loaded = loadHardcodedInverseCDF(opts, H, run_cdf_copy_path, invcdf_mean);
        writeHardcodedMeanParams(run_dir, fine_outputs, opts, run_cdf_copy_path, loaded, H.qx_const);
        return true;
    }

    // Case 2: computed from fine loop
    if (!opts.upscale_only) {
        // IMPORTANT: union-t mean (NEW)
        invcdf_mean = mean_ts_union_t(fine_outputs.inverse_qx_cdfs);
        invcdf_mean.writefile(run_cdf_copy_path);
        writeComputedMeanParams(run_dir, P, fine_outputs);
        return true;
    }

    // Case 3: strict upscale-only load from INPUT SOURCE (resume_run_dir)
    if (resume_run_dir.empty()) {
        std::cerr << "ERROR: strict --upscale-only requires resume_run_dir as INPUT source.\n";
        return false;
    }

    if (!loadMeanParamsForUpscaleOnly_FromInputSource(resume_run_dir, fine_outputs, invcdf_mean)) {
        MPI_Abort(PETSC_COMM_WORLD, 99);
        return false;
    }

    // copy mean CDF into OUTPUT run folder
    invcdf_mean.writefile(run_cdf_copy_path);

    // provenance
    {
        std::ofstream f(joinPath(run_dir, "input_source.txt"));
        f << "input_source_run_dir=" << resume_run_dir << "\n";
    }

    return true;
}

static void computeMeanBTCs(
    const std::vector<double>& xLocations,
    const std::vector<TimeSeriesSet<double>>& stacks,
    TimeSeriesSet<double>& mean_out)
{
    mean_out.clear();
    for (int i = 0; i < (int)xLocations.size(); ++i) {
        // IMPORTANT: union-t mean (NEW)
        mean_out.append(mean_ts_union_t(stacks[i]), fmt_x(xLocations[i]));
    }
}

static void broadcastMeanParameters(
    FineScaleOutputs& outputs,
    TimeSeries<double>& invcdf_mean,
    int rank)
{
    MPI_Bcast(&outputs.lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.velocity_lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.velocity_ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.nu_x_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.nu_y_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.lambda_K_x_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.lambda_K_y_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
}

// ============================================================================
// Upscaled run (transport + PT separate)
// ============================================================================

static bool run_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    double lc_mean,
    double lx_mean,
    double ly_mean,
    double nu_x_mean,
    double nu_y_mean,
    double dt_mean,
    const TimeSeries<double>& invcdf_mean,
    std::string& up_dir_out,
    std::string& up_btc_path_out,
    std::string& up_btc_deriv_path_out,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs_transport,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_pdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_cdf,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;

    std::string up_dir = joinPath(run_dir, "upscaled_mean");
    if (rank == 0) createDirectory(up_dir);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (!up_dir.empty() && up_dir.back() != '/' && up_dir.back() != '\\') up_dir += "/";

    const std::string up_pfx = "upscaled_";

    Grid2D g_u(nx, ny, Lx, Ly);

    if (P.CorrelationModel == SimParams::correlationmode::exponentialfit)
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::exponential_vdep;
    else if (P.CorrelationModel == SimParams::correlationmode::gaussian)
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::gaussian;
    else if (P.CorrelationModel == SimParams::correlationmode::matern)
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::matern;
    else
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::oneoversum;

    g_u.setDiffusionFactor(P.diffusion_factor);

    int qx_size = nu * (nx + 1);
    int qy_size = (nu + 1) * nx;

    auto& qx_u = g_u.flux("qx");
    auto& qy_u = g_u.flux("qy");

    qx_u.resize(qx_size, 0.0);
    qy_u.resize(qy_size, 0.0);

    for (int j = 0; j < nu; ++j) {
        double u = (nu <= 1) ? 0.0 : static_cast<double>(j) / (nu - 1);
        double v_at_u = invcdf_mean.interpol(u);

        for (int i = 0; i < nx + 1; ++i) {
            int id = j * (nx + 1) + i;
            if (id >= qx_size) return false;
            qx_u[id] = v_at_u;
        }
    }

    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));

    g_u.setMixingParams(lc_mean, lx_mean, ly_mean, nu_x_mean, nu_y_mean);

    g_u.SetVal("diffusion", P.Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
    g_u.setBTCLocations(P.xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    if (opts.solve_upscale_transport) {
        g_u.SolveTransport(P.t_end_pdf, dt_mean, "transport_", 500, up_dir, "Cu", &BTCs_Upscaled);

        for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_Upscaled.size(); ++i) {
            BTCs_Upscaled.setSeriesName(i, fmt_x(P.xLocations[i]));
            Fine_Scale_BTCs_transport[i].append(BTCs_Upscaled[i].derivative(), "Upscaled" + fmt_x(P.xLocations[i]));
        }

        BTCs_Upscaled.write(up_btc_path);
        BTCs_Upscaled.derivative().write(up_btc_deriv_path);
    }

    // PT (SEPARATE)
    if (opts.perform_upscaled_PT) {
        const int n_particles = 100000;
        const double dx_step = 0.01;

        PathwaySet pathways;
        pathways.Initialize(n_particles, PathwaySet::Weighting::FluxWeighted, &g_u);
        pathways.trackAllPathwaysWithDiffusion(&g_u, dx_step, g_u.getDiffusion());

        const std::vector<double>& btc_locations = g_u.getBTCLocations();
        if (!btc_locations.empty()) {
            TimeSeriesSet<double> btc_pdf = pathways.getBreakthroughCurve(
                btc_locations, true, PathwaySet::BTCType::PDF, 200);
            btc_pdf.write(joinPath(up_dir, up_pfx + "PT_btc_pdf.csv"));

            TimeSeriesSet<double> btc_cdf = pathways.getBreakthroughCurve(
                btc_locations, true, PathwaySet::BTCType::CDF);
            btc_cdf.write(joinPath(up_dir, up_pfx + "PT_btc_cdf.csv"));

            const int nmin = std::min((int)P.xLocations.size(), (int)btc_pdf.size());
            for (int i = 0; i < nmin; ++i)
                Fine_Scale_PT_pdf[i].append(btc_pdf[i], "Upscaled_PT" + fmt_x(P.xLocations[i]));

            const int nmin2 = std::min((int)P.xLocations.size(), (int)btc_cdf.size());
            for (int i = 0; i < nmin2; ++i)
                Fine_Scale_PT_cdf[i].append(btc_cdf[i], "Upscaled_PT" + fmt_x(P.xLocations[i]));
        }
    }

    up_dir_out = up_dir;
    up_btc_path_out = up_btc_path;
    up_btc_deriv_path_out = up_btc_deriv_path;
    return true;
}

// ============================================================================
// Main orchestrator
// ============================================================================

bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& resume_run_dir,
    RunOutputs& out,
    int rank)
{
    // Transport compare
    out.Fine_Scale_BTCs.clear();
    out.Fine_Scale_BTCs.resize(P.xLocations.size());

    // PT compare (separate)
    out.Fine_Scale_PT_pdf.clear();
    out.Fine_Scale_PT_cdf.clear();
    out.Fine_Scale_PT_pdf.resize(P.xLocations.size());
    out.Fine_Scale_PT_cdf.resize(P.xLocations.size());

    FineScaleOutputs fine_outputs;
    fine_outputs.BTCs.resize(P.xLocations.size());
    fine_outputs.PT_pdfs.resize(P.xLocations.size());
    fine_outputs.PT_cdfs.resize(P.xLocations.size());

    if (!opts.upscale_only) {
        const bool ok = run_fine_loop_collect(P, opts, out.run_dir, fine_outputs, rank);
        if (!ok) return false;
    }

    out.Fine_Scale_BTCs   = fine_outputs.BTCs;
    out.Fine_Scale_PT_pdf = fine_outputs.PT_pdfs;
    out.Fine_Scale_PT_cdf = fine_outputs.PT_cdfs;

    TimeSeries<double> invcdf_mean;
    const bool ok_mean = build_mean_for_upscaled(P, opts, H, out.run_dir, resume_run_dir, fine_outputs, invcdf_mean, rank);
    if (!ok_mean) return false;

    broadcastMeanParameters(fine_outputs, invcdf_mean, rank);

    std::string up_dir, up_btc_path, up_btc_deriv_path;
    const bool ok_up = run_upscaled(
        P, opts, out.run_dir,
        fine_outputs.lc_mean,
        fine_outputs.velocity_lx_mean,
        fine_outputs.velocity_ly_mean,
        fine_outputs.nu_x_mean,
        fine_outputs.nu_y_mean,
        fine_outputs.dt_mean,
        invcdf_mean,
        up_dir, up_btc_path, up_btc_deriv_path,
        out.Fine_Scale_BTCs,
        out.Fine_Scale_PT_pdf,
        out.Fine_Scale_PT_cdf,
        rank
    );
    if (!ok_up) return false;

    out.up_dir = up_dir;
    out.up_btc_path = up_btc_path;
    out.up_btc_deriv_path = up_btc_deriv_path;

    // Means (SEPARATED) — union-t mean used inside computeMeanBTCs now
    computeMeanBTCs(P.xLocations, out.Fine_Scale_BTCs,   out.mean_transport_full);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_PT_pdf, out.mean_pt_pdf);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_PT_cdf, out.mean_pt_cdf);

    out.mean_BTCs = out.mean_transport_full;

    // Write per-x PT compare files
    if (rank == 0) {
        for (int i = 0; i < (int)P.xLocations.size(); ++i) {
            out.Fine_Scale_PT_pdf[i].write(joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "PT_Compare_pdf.csv"));
            out.Fine_Scale_PT_cdf[i].write(joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "PT_Compare_cdf.csv"));
        }
    }

    return true;
}
