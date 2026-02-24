// sim_runner.h
#pragma once

#include <string>
#include <vector>
#include <utility>   // std::pair

// Keep sim_runner.h LIGHT.
// Do NOT include grid.h here (it causes heavy include chains / incomplete type issues).
#include "TimeSeries.h"
#include "TimeSeriesSet.h"

// -----------------------------
// Simulation switches / options
// -----------------------------
struct RunOptions
{
    // run modes
    bool upscale_only   = false;   // skip fine loop; do only upscaled (requires input dir unless hardcoded_mean)
    bool hardcoded_mean = false;   // use hardcoded lc/lx/ly/dt; and optionally load qx inverse-CDF from file

    // transport toggles
    bool solve_fine_scale_transport = true;
    bool solve_upscale_transport    = true;

    // particle tracking toggles
    bool perform_particle_tracking  = true; // fine-scale PT (per realization)
    bool perform_upscaled_PT        = true; // upscaled PT

    // When hardcoded_mean=true, you can supply a qx inverse-CDF file (u,v).
    // Priority used by sim_runner.cpp (hardcoded_mean case):
    //   1) load mean_qx_inverse_cdf.txt from this path (if exists / valid)
    //   2) else, if "qx_inverse_cdfs.txt" exists next to it, compute mean_qx_inverse_cdf.txt and load it
    //   3) else fallback to HardcodedMean::qx_const
    //
    // NOTE: sim_runner.cpp will ALSO copy/write mean_qx_inverse_cdf.txt into the new run_dir.
    std::string hardcoded_qx_cdf_path;

    // -----------------------------
    // Wiener particle tracker (1D/2D)
    // -----------------------------
    bool wiener_enable = false;          // default OFF
    std::string wiener_mode = "2d";      // "1dx" | "1dy" | "2d"
    double wiener_Dx = 0.0;
    double wiener_rx = 0.0;
    double wiener_ry = 0.0;
    double wiener_dt = 1e-3;
    unsigned long wiener_seed = 12345UL;
    std::string wiener_release = "left-flux"; // "left-uniform" | "left-flux" | "center"
};

struct FineScaleOutputs
{
    // Correlation time series sets (all realizations)
    TimeSeriesSet<double> velocity_x_correlations;
    TimeSeriesSet<double> velocity_y_correlations;
    TimeSeriesSet<double> advective_correlations;
    TimeSeriesSet<double> K_x_correlations;
    TimeSeriesSet<double> K_y_correlations;

    // Inverse CDFs and PDFs
    TimeSeriesSet<double> inverse_qx_cdfs;
    TimeSeriesSet<double> qx_pdfs;

    // Fitted parameters per realization
    std::vector<double> lc_all;
    std::vector<double> velocity_lx_all;
    std::vector<double> velocity_ly_all;
    std::vector<double> K_lx_all;
    std::vector<double> K_ly_all;
    std::vector<double> dt_all;

    // Mean parameters
    double nu_x_mean = 1.5;
    double nu_y_mean = 1.5;
    double lambda_K_x_mean = 0.0;
    double lambda_K_y_mean = 0.0;
    double lc_mean = 0.0;
    double velocity_lx_mean = 0.0;
    double velocity_ly_mean = 0.0;
    double dt_mean = 0.0;

    // TRANSPORT compare stacks (per-x), stores derivative curves
    std::vector<TimeSeriesSet<double>> BTCs;

    // PT compare stacks (SEPARATE)
    std::vector<TimeSeriesSet<double>> PT_pdfs;
    std::vector<TimeSeriesSet<double>> PT_cdfs;
};

struct HardcodedMean
{
    double lc_mean = 0.0;
    double lx_mean = 0.0;
    double ly_mean = 0.0;
    double dt_mean = 0.0;
    double nu_x = 0.0;
    double nu_y = 0.0;

    // fallback only when RunOptions::hardcoded_qx_cdf_path is not provided or load fails
    double qx_const = 0.0;
};

struct SimParams
{
    int nx = 0;
    int ny = 0;
    int nu = 0;

    double Lx = 0.0;
    double Ly = 0.0;

    int nReal_default = 20;

    // random field params
    double correlation_ls_x = 0.0;
    double correlation_ls_y = 0.0;
    double stdev = 0.0;
    double g_mean = 0.0;
    double diffusion_factor = 1.0;

    unsigned long run_seed = 0UL;

    // transport / BTC params
    double Diffusion_coefficient = 0.0;
    double t_end_pdf = 0.0;

    std::vector<double> xLocations;

    enum class correlationmode { exponentialfit, derivative, gaussian, matern, oneoversums }
        CorrelationModel = correlationmode::derivative;

    std::pair<double,double> correlation_x_range = {0.001, 0.02};
    std::pair<double,double> correlation_y_range = {0.001, 0.02};
};

struct RunOutputs
{
    std::string run_dir;

    // TRANSPORT compare stacks (per-x)
    std::vector<TimeSeriesSet<double>> Fine_Scale_BTCs;

    // Backward-compatible mean container (TRANSPORT mean)
    TimeSeriesSet<double> mean_BTCs;

    // Preferred explicit transport mean
    TimeSeriesSet<double> mean_transport_full;

    // PT compare stacks (SEPARATE)
    std::vector<TimeSeriesSet<double>> Fine_Scale_PT_pdf;
    std::vector<TimeSeriesSet<double>> Fine_Scale_PT_cdf;

    // PT means (per-x)
    TimeSeriesSet<double> mean_pt_pdf;
    TimeSeriesSet<double> mean_pt_cdf;

    // upscaled outputs (paths)
    std::string up_dir;
    std::string up_btc_path;
    std::string up_btc_deriv_path;
};

// Creates/chooses a run directory (rank0 decides; broadcasts to all ranks).
// NOTE: In upscale-only mode, sim_runner.cpp ALWAYS creates a NEW output run_dir.
//       resume_run_dir is treated as INPUT source only.
std::string prepare_run_dir_mpi(
    const std::string& output_dir,
    const std::string& resume_run_dir,
    const RunOptions& opts,
    int rank,
    const std::string& run_tag = "");

// Runs fine loop (unless upscale_only) + mean building + upscaled run.
// IMPORTANT: For strict upscale-only (hardcoded_mean==false), mean files are loaded
//            from resume_run_dir (INPUT SOURCE), outputs are written to out.run_dir (NEW folder).
bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& resume_run_dir,
    RunOutputs& out,
    int rank);

// Optional Wiener diffusion runner
int runDiffusionSimulation(const RunOptions &opts, int realization, const std::string &output_dir);
