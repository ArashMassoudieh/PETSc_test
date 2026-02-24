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

#include "grid.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline double du_from_nu(int nu)
{
    return (nu > 0) ? (1.0 / double(nu)) : 1.0;
}

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

        // -------------------------------------------------------------
        // IMPORTANT CHANGE:
        //   In upscale-only mode, ALWAYS create a NEW output run_dir.
        //   resume_run_dir is treated as INPUT SOURCE only.
        // -------------------------------------------------------------
        if (opts.upscale_only) {

            // In strict upscale-only (no hardcoded mean), REQUIRE existing input folder.
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

// Load an inverse CDF from a file, tolerant to header/no-header.
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

// Build mean inverse-CDF from the multi-series output (qx_inverse_cdfs.txt)
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

        // Expect pairs: (t,q), (t,q), ...
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
// NEW: Particle-tracking mean arrival time helpers (per-realization + mean)
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

        // trapezoid on (t*y)
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

// Reads: x,mean_arrival_time (delimiter-robust)
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
// Helper functions for fine-scale loop
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

static std::pair<double, double> computeAndFitKCorrelations(
    Grid2D& g,
    const SimParams& P,
    const std::string& fine_dir,
    const std::string& pfx,
    FineScaleOutputs& outputs)
{
    int num_deltas = 30;
    int num_samples_per_delta = 10000;

    TimeSeries<double> K_corr_x, K_corr_y;

    // X-direction
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

    double lambda_K_x = (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ?
                            K_corr_x.fitExponentialDecay() :
                            (P.CorrelationModel == SimParams::correlationmode::derivative) ?
                                K_corr_x.getTime(0)/(1.0-K_corr_x.getValue(0)) :
                                (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                                    K_corr_x.fitGaussianDecay() :
                                    K_corr_x.fitMaternDecay().second;

    outputs.K_x_correlations.append(K_corr_x);

    // Y-direction
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

    double lambda_K_y = (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ?
                            K_corr_y.fitExponentialDecay() :
                            (P.CorrelationModel == SimParams::correlationmode::derivative) ?
                                K_corr_y.getTime(0)/(1.0-K_corr_y.getValue(0)) :
                                (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                                    K_corr_y.fitGaussianDecay() :
                                    K_corr_y.fitMaternDecay().second;

    outputs.K_y_correlations.append(K_corr_y);

    return {lambda_K_x, lambda_K_y};
}

static std::pair<double, double> computeAndFitVelocityCorrelations(
    Grid2D& g,
    const SimParams& P,
    const std::string& fine_dir,
    const std::string& pfx,
    FineScaleOutputs& outputs,
    int r)
{
    int num_deltas = 30;
    int num_samples_per_delta = 10000;

    TimeSeries<double> corr_x, corr_y;

    // X-direction
    for (int i = 0; i < num_deltas; ++i) {
        double exponent = static_cast<double>(i) / (num_deltas - 1);
        double delta = P.correlation_x_range.first *
                       std::pow(P.correlation_x_range.second / P.correlation_x_range.first, exponent);
        try {
            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                "qx_normal_score", Grid2D::ArrayKind::Fx,
                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
            corr_x.append(delta, samples.correlation_tc());
        } catch (...) {}
    }
    corr_x.writefile(joinPath(fine_dir, pfx + "velocity_correlation_x.txt"));
    outputs.velocity_x_correlations.append(corr_x, "Realization" + aquiutils::numbertostring(r + 1));

    double lambda_x = (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ?
                          corr_x.fitExponentialDecay() :
                          (P.CorrelationModel == SimParams::correlationmode::derivative) ?
                              corr_x.getTime(0)/(1.0-corr_x.getValue(0)) :
                              (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                                  corr_x.fitGaussianDecay() :
                                  corr_x.fitMaternDecay().second;

    // Y-direction
    for (int i = 0; i < num_deltas; ++i) {
        double exponent = static_cast<double>(i) / (num_deltas - 1);
        double delta = P.correlation_y_range.first *
                       std::pow(P.correlation_y_range.second / P.correlation_y_range.first, exponent);
        try {
            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                "qx_normal_score", Grid2D::ArrayKind::Fx,
                num_samples_per_delta, delta, 0, PerturbDir::YOnly);
            corr_y.append(delta, samples.correlation_tc());
        } catch (...) {}
    }
    corr_y.writefile(joinPath(fine_dir, pfx + "velocity_correlation_y.txt"));
    outputs.velocity_y_correlations.append(corr_y, "Realization" + aquiutils::numbertostring(r + 1));

    double lambda_y = (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ?
                          corr_y.fitExponentialDecay() :
                          (P.CorrelationModel == SimParams::correlationmode::derivative) ?
                              corr_y.getTime(0)/(1.0-corr_y.getValue(0)) :
                              (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                                  corr_y.fitGaussianDecay() :
                                  corr_y.fitMaternDecay().second;

    return {lambda_x, lambda_y};
}

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

    //pathwaySet.writeCombinedVTK(joinPath(output_dir, "diffusionpaths.vtk"));
    times.writefile(joinPath(output_dir, "diffusiontimes.csv"));

    double mean     = times.mean();
    double stdev    = times.stddev();
    double meanlog  = times.mean_log();
    double stdevlog = times.log().stddev();
    auto [tmin, tmax] = pathwaySet.travelTimeRange();

    TimeSeries<double> fptDist = times.distribution(100, 0);
    fptDist.writefile(joinPath(output_dir, "FirstPassageTimeDistribution.csv"));

    // Write statistics file
    {
        std::ofstream f(joinPath(output_dir, "diffusion_stats.txt"));
        if (!f.is_open()) {
            std::cerr << "WARNING: could not write diffusion_stats.txt\n";
        } else {
            f << "# First Passage Time Statistics\n"
              << "# Brownian diffusion escape from ellipse\n"
              << "#\n"
              << "# Parameters\n"
              << "mode            = " << opts.wiener_mode << "\n"
              << "D               = " << wp.D << "\n"
              << "dt              = " << wp.dt << "\n"
              << "rx              = " << wp.rx << "\n"
              << "ry              = " << wp.ry << "\n"
              << "seed            = " << wp.seed << "\n"
              << "num_particles   = " << realization << "\n"
              << "#\n"
              << "# Results\n"
              << "mean            = " << mean << "\n"
              << "stdev           = " << stdev << "\n"
              << "mean_log        = " << meanlog << "\n"
              << "stdev_log       = " << stdevlog << "\n"
              << "min             = " << tmin << "\n"
              << "max             = " << tmax << "\n";
        }
    }

    std::cout << "  Mean FPT: " << mean << "  StdDev: " << stdev << "\n"
              << "  Mean(log): " << meanlog << "  StdDev(log): " << stdevlog << "\n"
              << "  Range: [" << tmin << ", " << tmax << "]\n"
              << "  Stats written to: " << joinPath(output_dir, "diffusion_stats.txt") << "\n";

    return 0;
}

static double computeAndFitAdvectiveCorrelation(
    Grid2D& g,
    const SimParams& P,
    const RunOptions& opts,
    const std::string& fine_dir,
    const std::string& pfx,
    FineScaleOutputs& outputs,
    int r)
{
    PathwaySet pathways;
    pathways.Initialize(10000, PathwaySet::Weighting::FluxWeighted, &g);
    pathways.trackAllPathwaysWithDiffusion(&g, 0.01, g.getDiffusion());

    // --- Breakthrough curves ---
    const std::vector<double>& btc_locations = g.getBTCLocations();
    if (!btc_locations.empty()) {
        TimeSeriesSet<double> btc_pdf = pathways.getBreakthroughCurve(
            btc_locations, true, PathwaySet::BTCType::PDF, 200);
        btc_pdf.write(joinPath(fine_dir, pfx + "btc_pdf.csv"));

        // --- NEW: per-realization PT mean arrival time (from PDF) ---
        {
            const std::string mean_csv = joinPath(fine_dir, pfx + "PT_mean_arrival_time.csv");
            std::ofstream mf(mean_csv);
            mf << "x,mean_arrival_time\n";
            mf << std::setprecision(15);

            for (int i = 0; i < (int)btc_locations.size() && i < (int)btc_pdf.size(); ++i) {
                const double mu = mean_time_from_pdf(btc_pdf[i]);
                mf << btc_locations[i] << "," << mu << "\n";
            }
        }

        TimeSeriesSet<double> btc_cdf = pathways.getBreakthroughCurve(
            btc_locations, true, PathwaySet::BTCType::CDF);
        btc_cdf.write(joinPath(fine_dir, pfx + "btc_cdf.csv"));
    }

    // --- Correlation analysis ---
    double Delta_x_min = 0.001, Delta_x_max = 0.5;
    int num_Delta_x = 30;
    int num_samples_per_Delta_x = 10000;
    TimeSeries<double> qx_correlation;
    for (int i = 0; i < num_Delta_x; ++i) {
        double exponent = static_cast<double>(i) / (num_Delta_x - 1);
        double Delta_x = Delta_x_min * std::pow(Delta_x_max / Delta_x_min, exponent);
        try {
            PathwaySet particle_pairs = pathways.sampleParticlePairs(Delta_x, num_samples_per_Delta_x);
            double correlation = particle_pairs.calculateCorrelation(0, 1, "qx");
            qx_correlation.append(Delta_x, correlation);
        } catch (...) {}
    }
    qx_correlation.writefile(joinPath(fine_dir, pfx + "qx_correlation_vs_distance.txt"));
    outputs.advective_correlations.append(qx_correlation, "Realization" + aquiutils::numbertostring(r + 1));
    pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
    pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"),1000);
    return qx_correlation.fitExponentialDecay();
}

static void writeRealizationMeta(
    const std::string& fine_dir,
    const std::string& pfx,
    const SimParams& P,
    int r,
    int nx, int ny, int nu,
    double Lx, double Ly,
    double lc_emp,
    double lambda_x_emp,
    double lambda_y_emp,
    double lambda_K_x_emp,
    double lambda_K_y_emp,
    double dt_optimal,
    unsigned long seed)
{
    std::ofstream m(joinPath(fine_dir, pfx + "meta.txt"));
    m << "realization=" << r << "\n";
    m << "nx=" << nx << "\nny=" << ny << "\nnu=" << nu << "\n";
    m << "Lx=" << Lx << "\nLy=" << Ly << "\n";
    m << "D=" << P.Diffusion_coefficient << "\n";
    m << "correlation_ls_x=" << P.correlation_ls_x << "\n";
    m << "correlation_ls_y=" << P.correlation_ls_y << "\n";
    m << "correlation_model=" << (int)P.CorrelationModel << "\n";
    m << "stdev=" << P.stdev << "\n";
    m << "g_mean=" << P.g_mean << "\n";
    m << "lc=" << lc_emp << "\n";
    m << "dt_opt=" << dt_optimal << "\n";
    m << "seed=" << seed << "\n";

    m << "\n# Fitted K correlation lengths:\n";
    m << "lambda_K_x=" << lambda_K_x_emp << "\n";
    m << "lambda_K_y=" << lambda_K_y_emp << "\n";

    m << "\n# Fitted velocity correlation lengths:\n";
    m << "fitted_velocity_lambda_x=" << lambda_x_emp << "\n";
    m << "fitted_velocity_lambda_y=" << lambda_y_emp << "\n";

    m << "\n# Input correlation parameters:\n";
    m << "input_lambda_x=" << P.correlation_ls_x << "\n";
    m << "input_lambda_y=" << P.correlation_ls_y << "\n";
}

static void writeMeanCorrelations(
    const std::string& run_dir,
    const SimParams& P,
    FineScaleOutputs& outputs)
{
    // Write all correlation time series
    outputs.advective_correlations.write(joinPath(run_dir, "advective_correlations.txt"));
    outputs.velocity_x_correlations.write(joinPath(run_dir, "diffusion_x_correlations.txt"));
    outputs.velocity_y_correlations.write(joinPath(run_dir, "diffusion_y_correlations.txt"));
    outputs.K_x_correlations.write(joinPath(run_dir, "K_x_correlations.txt"));
    outputs.K_y_correlations.write(joinPath(run_dir, "K_y_correlations.txt"));

    // Mean correlation curves (raw)
    auto mean_corr_a = outputs.advective_correlations.mean_ts();
    auto mean_corr_x = outputs.velocity_x_correlations.mean_ts();
    auto mean_corr_y = outputs.velocity_y_correlations.mean_ts();
    auto mean_K_corr_x = outputs.K_x_correlations.mean_ts();
    auto mean_K_corr_y = outputs.K_y_correlations.mean_ts();

    mean_corr_a.writefile(joinPath(run_dir, "advective_correlations_mean.txt"));
    mean_corr_x.writefile(joinPath(run_dir, "diffusion_x_correlations_mean.txt"));
    mean_corr_y.writefile(joinPath(run_dir, "diffusion_y_correlations_mean.txt"));
    mean_K_corr_x.writefile(joinPath(run_dir, "K_x_correlations_mean.txt"));
    mean_K_corr_y.writefile(joinPath(run_dir, "K_y_correlations_mean.txt"));

    // Fit and write advective (always exponential)
    double l_a = mean_corr_a.fitExponentialDecay();
    TimeSeries<double> fitted_a;
    for (size_t i = 0; i < mean_corr_a.size(); ++i) {
        double r = mean_corr_a.getTime(i);
        fitted_a.append(r, std::exp(-r / l_a));
    }
    fitted_a.writefile(joinPath(run_dir, "advective_correlations_mean_fitted.txt"));

    // Fit and write velocity correlations
    auto fitAndWrite = [&](const TimeSeries<double>& mean_corr, const std::string& filename, double& nu_out) {
        double ell = 0, nu = 0.0;
        if (P.CorrelationModel == SimParams::correlationmode::exponentialfit)
            ell = mean_corr.fitExponentialDecay();
        else if (P.CorrelationModel == SimParams::correlationmode::derivative)
            ell = mean_corr.getTime(0) / (1.0 - mean_corr.getValue(0));
        else if (P.CorrelationModel == SimParams::correlationmode::gaussian)
            ell = mean_corr.fitGaussianDecay();
        else {
            auto [n, e] = mean_corr.fitMaternDecay();
            nu = n; ell = e;
        }

        TimeSeries<double> fitted;
        for (size_t i = 0; i < mean_corr.size(); ++i) {
            double r = mean_corr.getTime(i);
            double rho = (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
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

    // Fit and write K correlations
    auto fitKAndWrite = [&](const TimeSeries<double>& mean_K_corr, const std::string& filename) {
        double ell_K = (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                           mean_K_corr.fitGaussianDecay() :
                           (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ?
                               mean_K_corr.fitExponentialDecay() :
                               (P.CorrelationModel == SimParams::correlationmode::derivative) ?
                                   mean_K_corr.getTime(0) / (1.0 - mean_K_corr.getValue(0)) :
                                   mean_K_corr.fitMaternDecay().second;

        TimeSeries<double> fitted_K;
        for (size_t i = 0; i < mean_K_corr.size(); ++i) {
            double r = mean_K_corr.getTime(i);
            double rho = (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
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
    outputs.lc_mean = mean_of(outputs.lc_all);
    outputs.velocity_lx_mean = mean_of(outputs.velocity_lx_all);
    outputs.velocity_ly_mean = mean_of(outputs.velocity_ly_all);
    outputs.dt_mean = mean_of(outputs.dt_all);
}

static void writeCorrelationSummary(
    const std::string& run_dir,
    const SimParams& P,
    const FineScaleOutputs& outputs)
{
    std::ofstream summary(joinPath(run_dir, "correlation_summary.txt"));
    summary << "=== Correlation Length Scale Summary ===\n\n";

    summary << "Input (K field generation):\n";
    summary << "  lambda_x_input = " << P.correlation_ls_x << "\n";
    summary << "  lambda_y_input = " << P.correlation_ls_y << "\n";
    summary << "  correlation_model = " << (int)P.CorrelationModel << "\n";

    summary << "\nFitted from K field:\n";
    summary << "  lambda_K_x_mean = " << outputs.lambda_K_x_mean << "\n";
    summary << "  lambda_K_y_mean = " << outputs.lambda_K_y_mean << "\n";
    summary << "  K_x ratio (fitted/input) = " << outputs.lambda_K_x_mean / P.correlation_ls_x << "\n";
    summary << "  K_y ratio (fitted/input) = " << outputs.lambda_K_y_mean / P.correlation_ls_y << "\n";

    summary << "\nFitted from velocity field:\n";
    summary << "  lambda_velocity_x_mean = " << outputs.velocity_lx_mean << "\n";
    summary << "  lambda_velocity_y_mean = " << outputs.velocity_ly_mean << "\n";

    summary << "\nAdvective:\n";
    summary << "  lc_mean = " << outputs.lc_mean << "\n";

    summary << "\nGrid resolution:\n";
    summary << "  nx=" << P.nx << ", ny=" << P.ny << "\n";
    summary << "  dx=" << P.Lx/P.nx << ", dy=" << P.Ly/P.ny << "\n";
    summary << "  Points per correlation length (input):\n";
    summary << "    x: " << P.correlation_ls_x / (P.Lx/P.nx) << "\n";
    summary << "    y: " << P.correlation_ls_y / (P.Ly/P.ny) << "\n";
}

// ============================================================================
// Main loop (refactored)
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

    const int nReal = P.nReal_default;

    for (int r = 1; r <= nReal; ++r) {
        const std::string rlab = makeRealLabel(r);
        std::string fine_dir = joinPath(run_dir, makeFineFolder(r));
        if (rank == 0) createDirectory(fine_dir);
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        const std::string pfx = rlab + "_";
        Grid2D g(nx, ny, Lx, Ly);

        PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;
        unsigned long seed = P.run_seed + 1000UL*(unsigned long)r + (unsigned long)rank;

        // Generate K field
        generateKField(g, P, seed);

        // Compute K correlations
        auto [lambda_K_x_emp, lambda_K_y_emp] = computeAndFitKCorrelations(g, P, fine_dir, pfx, outputs);

        // Normalize and create exponential field
        g.normalizeField("K_normal_score", Grid2D::ArrayKind::Cell);
        PetscTime(&t_asm0);
        PetscTime(&t_total0);

        g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
        g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.vti"));

        g.createExponentialField("K_normal_score", P.stdev, P.g_mean, "K");
        g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));

        // Solve Darcy flow
        PetscTime(&t_solve0);
        g.DarcySolve(Lx, 0, "K", "K");
        std::cout << "Darcy solved ... " << std::endl;
        PetscTime(&t_solve1);

        g.computeMassBalanceError("MassBalanceError");
        g.writeNamedMatrix("MassBalanceError", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "error.txt"));

        g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.vti"));
        g.writeNamedMatrix("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.txt"));
        g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.vti"));
        g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.vti"));

        // Compute velocity normal scores
        TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
        TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
        g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);
        g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx_normal_score.vti"));

        // Compute velocity correlations
        auto [lambda_x_emp, lambda_y_emp] = computeAndFitVelocityCorrelations(g, P, fine_dir, pfx, outputs, r);

        // Store velocity and K fitted values
        outputs.velocity_lx_all.push_back(lambda_x_emp);
        outputs.velocity_ly_all.push_back(lambda_y_emp);
        outputs.K_lx_all.push_back(lambda_K_x_emp);
        outputs.K_ly_all.push_back(lambda_K_y_emp);

        // Inverse CDF + PDF
        TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100, 1e-6);
        TimeSeries<double> qx_pdf = g.extractFieldPDF("qx", Grid2D::ArrayKind::Fx, 50, 1e-6, true);
        qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
        qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));
        qx_pdf.writefile(joinPath(fine_dir, pfx + "qx_pdf.txt"));

        // Optimal dt
        double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;
        outputs.dt_all.push_back(dt_optimal);

        // Transport
        g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
        g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.txt"));
        g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.txt"));
        g.writeNamedMatrix("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.txt"));

        g.SetVal("diffusion", P.Diffusion_coefficient);
        g.assignConstant("D_y", Grid2D::ArrayKind::Fy, P.Diffusion_coefficient);
        g.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "D_y.vti"));
        g.SetVal("porosity", 1.0);
        g.SetVal("c_left", 1.0);

        TimeSeriesSet<double> BTCs_FineScaled;
        g.setBTCLocations(P.xLocations);

        if (opts.solve_fine_scale_transport) {
            g.SolveTransport(P.t_end_pdf, std::min(dt_optimal, 0.5 / 10.0),
                             "transport_", 500, fine_dir, "C", &BTCs_FineScaled, r);

            for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_FineScaled.size(); ++i) {
                BTCs_FineScaled.setSeriesName(i, fmt_x(P.xLocations[i]));
                outputs.BTCs[i].append(BTCs_FineScaled[i].derivative(), pfx + fmt_x(P.xLocations[i]));
            }

            BTCs_FineScaled.write(joinPath(fine_dir, pfx + "BTC_FineScaled.csv"));
            BTCs_FineScaled.derivative().write(joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv"));
        }

        PetscTime(&t_total1);

        // Compute advective correlation (also writes PT mean per realization)
        double lc_emp = computeAndFitAdvectiveCorrelation(g, P, opts, fine_dir, pfx, outputs, r);
        outputs.lc_all.push_back(lc_emp);

        // Write meta file
        if (rank == 0) {
            writeRealizationMeta(fine_dir, pfx, P, r, nx, ny, nu, Lx, Ly,
                                 lc_emp, lambda_x_emp, lambda_y_emp,
                                 lambda_K_x_emp, lambda_K_y_emp, dt_optimal, seed);
        }

        // Accumulate for CSV
        if (rank == 0) {
            outputs.inverse_qx_cdfs.append(qx_inverse_cdf, "qx_inverse_cdf" + aquiutils::numbertostring(r + 1));
            outputs.qx_pdfs.append(qx_pdf, "qx_pdf" + aquiutils::numbertostring(r + 1));

            std::ofstream f(stats_csv, std::ios::app);
            f << r << "," << lc_emp << "," << lambda_x_emp << "," << lambda_y_emp << ","
              << lambda_K_x_emp << "," << lambda_K_y_emp << "," << dt_optimal << "\n";
        }

        if (rank == 0) {
            std::cout << "[Fine " << rlab << "] Assembly time: " << (t_asm1 - t_asm0) << " s\n";
            std::cout << "[Fine " << rlab << "] Solve time:    " << (t_solve1 - t_solve0) << " s\n";
            std::cout << "[Fine " << rlab << "] Total time:   " << (t_total1 - t_total0) << " s\n";
            std::cout << "[Fine " << rlab << "] Outputs saved to: " << fine_dir << "\n";
        }

        MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Post-processing: write all mean correlations
    if (rank == 0) {
        outputs.inverse_qx_cdfs.write(joinPath(run_dir, "qx_inverse_cdfs.txt"));
        outputs.qx_pdfs.write(joinPath(run_dir, "qx_pdfs.txt"));
        outputs.qx_pdfs.mean_ts().writefile(joinPath(run_dir, "qx_mean_pdf.txt"));

        writeMeanCorrelations(run_dir, P, outputs);
        computeFinalMeans(outputs);
        writeCorrelationSummary(run_dir, P, outputs);

        // --- NEW: mean PT arrival time over realizations (fine loop) ---
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

            std::cout << "Wrote mean PT arrival times over realizations: " << out_mean_csv << "\n";
        }
    }

    return true;
}

// ============================================================================
// Helper functions for mean parameter building
// ============================================================================

static bool loadHardcodedInverseCDF(
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_cdf_copy_path,
    TimeSeries<double>& invcdf_mean)
{
    bool loaded = false;

    // (1) Try load preferred mean file
    loaded = try_load_hardcoded_invcdf(invcdf_mean, opts.hardcoded_qx_cdf_path);

    // (2) If mean missing, compute from qx_inverse_cdfs.txt (if exists)
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
            // If mean file exists now (maybe created earlier), try load it anyway
            if (!loaded && fileExists(meanf)) {
                loaded = try_load_hardcoded_invcdf(invcdf_mean, meanf);
            }
        }
    }

    // (3) Flat fallback
    if (!loaded) {
        build_flat_invcdf(invcdf_mean, H.qx_const);
    }

    // Always write a copy into this NEW run_dir
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

    f << "\n# Input parameters:\n";
    f << "input_correlation_ls_x=" << P.correlation_ls_x << "\n";
    f << "input_correlation_ls_y=" << P.correlation_ls_y << "\n";
    f << "correlation_model=" << (int)P.CorrelationModel << "\n";
}

// -----------------------------------------------------------------------------
// STRICT upscale-only mean loading must use INPUT SOURCE folder, not run_dir.
// -----------------------------------------------------------------------------
static bool loadMeanParamsForUpscaleOnly_FromInputSource(
    const std::string& input_dir,
    FineScaleOutputs& outputs,
    TimeSeries<double>& invcdf_mean)
{
    const std::string mean_params_path = joinPath(input_dir, "mean_params.txt");
    const std::string mean_cdf_path = joinPath(input_dir, "mean_qx_inverse_cdf.txt");

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
    std::cout << "Loaded mean_params.txt and mean_qx_inverse_cdf.txt from INPUT folder\n";
    return true;
}

// ============================================================================
// Main mean-building function (refactored)
// ============================================================================

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
    if (rank != 0) return true; // only rank0 builds; broadcast later

    const std::string run_cdf_copy_path = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

    // ========================================================================
    // Case 1: Using hardcoded mean values
    // ========================================================================
    if (opts.hardcoded_mean) {
        // Set mean parameters from hardcoded values
        fine_outputs.lc_mean = H.lc_mean;
        fine_outputs.velocity_lx_mean = H.lx_mean;
        fine_outputs.velocity_ly_mean = H.ly_mean;
        fine_outputs.dt_mean = H.dt_mean;
        fine_outputs.nu_x_mean = H.nu_x;
        fine_outputs.nu_y_mean = H.nu_y;
        // Note: lambda_K_x_mean and lambda_K_y_mean not in HardcodedMean struct

        // Load or build inverse CDF
        bool loaded = loadHardcodedInverseCDF(opts, H, run_cdf_copy_path, invcdf_mean);

        // Write parameters to file
        writeHardcodedMeanParams(run_dir, fine_outputs, opts, run_cdf_copy_path, loaded, H.qx_const);

        // Log to console
        if (loaded) {
            std::cout << "Using HARD-CODED mean params + qx inverse CDF (loaded/derived)\n";
            std::cout << "Copied to: " << run_cdf_copy_path << "\n";
        } else {
            std::cout << "Using HARD-CODED mean params + CONSTANT qx (no mean/multi files).\n";
            std::cout << "Wrote flat CDF to: " << run_cdf_copy_path << "\n";
        }

        return true;
    }

    // ========================================================================
    // Case 2: Compute mean from current fine run
    // ========================================================================
    if (!opts.upscale_only) {
        // Means already computed in run_fine_loop_collect
        // (fine_outputs.lc_mean, velocity_lx_mean, etc. are already set)

        // Compute mean inverse CDF
        invcdf_mean = fine_outputs.inverse_qx_cdfs.mean_ts();
        invcdf_mean.writefile(run_cdf_copy_path);

        // Write mean parameters to file
        writeComputedMeanParams(run_dir, P, fine_outputs);

        return true;
    }

    // ========================================================================
    // Case 3: Upscale-only mode - load existing mean files from INPUT SOURCE
    // ========================================================================
    if (resume_run_dir.empty()) {
        std::cerr << "ERROR: strict --upscale-only requires resume_run_dir as INPUT source.\n";
        return false;
    }

    if (!loadMeanParamsForUpscaleOnly_FromInputSource(resume_run_dir, fine_outputs, invcdf_mean)) {
        MPI_Abort(PETSC_COMM_WORLD, 99);
        return false;
    }

    // Copy mean CDF into OUTPUT run folder (self-contained)
    invcdf_mean.writefile(run_cdf_copy_path);

    // Record provenance (optional)
    {
        std::ofstream f(joinPath(run_dir, "input_source.txt"));
        f << "input_source_run_dir=" << resume_run_dir << "\n";
    }

    return true;
}

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
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;

    std::string up_dir = joinPath(run_dir, "upscaled_mean");
    if (rank == 0) createDirectory(up_dir);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (!up_dir.empty() && up_dir.back() != '/' && up_dir.back() != '\\') up_dir += "/";

    const std::string up_pfx = std::string("upscaled_");

    if (rank == 0) {
        std::cout << "\n=== UPSCALED RUN USING MEAN PARAMETERS ===\n";
        std::cout << "lc_mean=" << lc_mean << "\n";
        std::cout << "lambda_x_mean=" << lx_mean << "\n";
        std::cout << "lambda_y_mean=" << ly_mean << "\n";
        std::cout << "dt_mean=" << dt_mean << "\n";
        std::cout << "upscaled output: " << up_dir << "\n";
    }

    Grid2D g_u(nx, ny, Lx, Ly);

    if (P.CorrelationModel == SimParams::correlationmode::exponentialfit) {
        g_u.VelocityCorrelationModel =
            Grid2D::velocity_correlation_model::exponential;
    }
    else if (P.CorrelationModel == SimParams::correlationmode::gaussian) {
        g_u.VelocityCorrelationModel =
            Grid2D::velocity_correlation_model::gaussian;
    }
    else if (P.CorrelationModel == SimParams::correlationmode::matern) {
        g_u.VelocityCorrelationModel =
            Grid2D::velocity_correlation_model::matern;
    }
    else {
        g_u.VelocityCorrelationModel =
            Grid2D::velocity_correlation_model::oneoversum;
    }

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
            if (id >= qx_size) {
                std::cerr << "ERROR: qx index out of bounds: " << id << " >= " << qx_size << "\n";
                return false;
            }
            qx_u[id] = v_at_u;
        }
    }

    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.txt"));

    g_u.setMixingParams(lc_mean, lx_mean, ly_mean, nu_x_mean, nu_y_mean);

    const double dt_pdf = dt_mean;
    const int output_interval_pdf = 500;

    g_u.SetVal("diffusion", P.Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));
    g_u.writeNamedMatrix("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.txt"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
    g_u.setBTCLocations(P.xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    if (opts.solve_upscale_transport) {
        g_u.SolveTransport(P.t_end_pdf, dt_pdf, "transport_", output_interval_pdf, up_dir, "Cu", &BTCs_Upscaled);

        for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_Upscaled.size(); ++i) {
            BTCs_Upscaled.setSeriesName(i, fmt_x(P.xLocations[i]));
            Fine_Scale_BTCs[i].append(BTCs_Upscaled[i].derivative(), "Upscaled" + fmt_x(P.xLocations[i]));
        }

        BTCs_Upscaled.write(up_btc_path);
        BTCs_Upscaled.derivative().write(up_btc_deriv_path);

        g_u.writeNamedVTI_Auto("C", joinPath(up_dir, "Cu.vti"));
        g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell, joinPath(up_dir, up_pfx + "Cu.txt"));
    }

    if (opts.perform_upscaled_PT) {
        const int n_particles = 100000;
        const double dx_step = 0.01;

        PathwaySet pathways;
        pathways.Initialize(n_particles, PathwaySet::Weighting::FluxWeighted, &g_u);
        pathways.trackAllPathwaysWithDiffusion(&g_u, dx_step, g_u.getDiffusion());

        // Breakthrough curves
        const std::vector<double>& btc_locations = g_u.getBTCLocations();
        if (!btc_locations.empty()) {
            TimeSeriesSet<double> btc_pdf = pathways.getBreakthroughCurve(
                btc_locations, true, PathwaySet::BTCType::PDF, 200);
            btc_pdf.write(joinPath(up_dir, up_pfx + "PT_btc_pdf.csv"));

            TimeSeriesSet<double> btc_cdf = pathways.getBreakthroughCurve(
                btc_locations, true, PathwaySet::BTCType::CDF);
            btc_cdf.write(joinPath(up_dir, up_pfx + "PT_btc_cdf.csv"));

            // Append PT derivative BTCs to fine-scale comparison set
            for (int i = 0; i < (int)btc_locations.size() && i < (int)btc_pdf.size(); ++i) {
                Fine_Scale_BTCs[i].append(btc_pdf[i], "Upscaled_PT" + fmt_x(btc_locations[i]));
            }
        }

        // Save 1000 pathways to VTK
        pathways.writeCombinedVTK(joinPath(up_dir, up_pfx + "pathways.vtk"), 1000);

        if (rank == 0)
            std::cout << "Upscaled particle tracking: " << n_particles
                      << " particles, " << pathways.size() << " completed\n";
    }

    up_dir_out = up_dir;
    up_btc_path_out = up_btc_path;
    up_btc_deriv_path_out = up_btc_deriv_path;
    return true;
}

// ============================================================================
// Helper functions for simulation blocks
// ============================================================================

static void computeMeanBTCs(
    const std::vector<double>& xLocations,
    const std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs,
    TimeSeriesSet<double>& mean_BTCs)
{
    mean_BTCs.clear();
    for (int i = 0; i < (int)xLocations.size(); ++i) {
        mean_BTCs.append(Fine_Scale_BTCs[i].mean_ts(), fmt_x(xLocations[i]));
    }
}

static void broadcastMeanParameters(
    FineScaleOutputs& outputs,
    TimeSeries<double>& invcdf_mean,
    int rank)
{
    // Broadcast scalar mean parameters
    MPI_Bcast(&outputs.lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.velocity_lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.velocity_ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.nu_x_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.nu_y_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.lambda_K_x_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.lambda_K_y_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // Broadcast inverse CDF
    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
}

// ============================================================================
// Main simulation orchestrator (refactored)
// ============================================================================

bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& resume_run_dir,
    RunOutputs& out,
    int rank)
{
    // Initialize outputs
    out.Fine_Scale_BTCs.clear();
    out.Fine_Scale_BTCs.resize(P.xLocations.size());

    FineScaleOutputs fine_outputs;
    fine_outputs.BTCs.resize(P.xLocations.size());

    // Run fine-scale simulations
    if (!opts.upscale_only) {
        const bool ok = run_fine_loop_collect(P, opts, out.run_dir, fine_outputs, rank);
        if (!ok) return false;
    }

    // Copy BTCs to output
    out.Fine_Scale_BTCs = fine_outputs.BTCs;

    // Compute mean BTCs
    computeMeanBTCs(P.xLocations, out.Fine_Scale_BTCs, out.mean_BTCs);

    // Build mean parameters for upscaled model
    TimeSeries<double> invcdf_mean;
    const bool ok_mean = build_mean_for_upscaled(P, opts, H, out.run_dir, resume_run_dir, fine_outputs, invcdf_mean, rank);
    if (!ok_mean) return false;

    // Broadcast parameters to all MPI ranks
    broadcastMeanParameters(fine_outputs, invcdf_mean, rank);

    // Run upscaled simulation
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
        rank
        );

    if (!ok_up) return false;

    // Store upscaled outputs
    out.up_dir = up_dir;
    out.up_btc_path = up_btc_path;
    out.up_btc_deriv_path = up_btc_deriv_path;

    return true;
}
