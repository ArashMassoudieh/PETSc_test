// sim_runner.h

#pragma once

#include <string>
#include <vector>

// Keep sim_runner.h LIGHT.
// Do NOT include grid.h here (it causes heavy include chains / incomplete type issues).
// Only include what we must for types used in signatures.
#include "TimeSeries.h"
#include "TimeSeriesSet.h"

// NOTE:
// Do NOT include sim_helpers.h here.
// sim_runner.cpp can include sim_helpers.h as needed.
// This avoids circular / heavy include chains (plotter.cpp includes sim_helpers.h).
// sim_runner.h is meant to stay lightweight.

// -----------------------------
// Simulation switches / options
// -----------------------------
struct RunOptions
{
    // run modes
    bool upscale_only   = false;   // skip fine loop; do only upscaled (requires run_dir unless hardcoded_mean)
    bool hardcoded_mean = false;   // use hardcoded lc/lx/ly/dt; and optionally load qx inverse-CDF from file

    // transport toggles
    bool solve_fine_scale_transport = true;
    bool solve_upscale_transport    = true;

    // When hardcoded_mean=true, you can supply a qx inverse-CDF file (u,v).
    // Priority used by sim_runner.cpp (hardcoded_mean case):
    //   1) load mean_qx_inverse_cdf.txt from this path (if exists / valid)
    //   2) else, if "qx_inverse_cdfs.txt" exists next to it, compute mean_qx_inverse_cdf.txt and load it
    //   3) else fallback to HardcodedMean::qx_const
    //
    // NOTE: sim_runner.cpp will ALSO copy/write mean_qx_inverse_cdf.txt into the new run_dir.
    std::string hardcoded_qx_cdf_path;
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
    std::vector<double> lc_all;           // advective correlation length
    std::vector<double> velocity_lx_all;  // velocity lambda_x
    std::vector<double> velocity_ly_all;  // velocity lambda_y
    std::vector<double> K_lx_all;         // K lambda_x
    std::vector<double> K_ly_all;         // K lambda_y
    std::vector<double> dt_all;           // optimal timesteps


    // Mean parameters
    double nu_x_mean = 1.5;
    double nu_y_mean = 1.5;
    double lambda_K_x_mean = 0.0;
    double lambda_K_y_mean = 0.0;
    double lc_mean = 0.0;              // advective correlation mean
    double velocity_lx_mean = 0.0;     // velocity lambda_x mean
    double velocity_ly_mean = 0.0;     // velocity lambda_y mean
    double dt_mean;

    // BTCs per location
    std::vector<TimeSeriesSet<double>> BTCs;
};

// -----------------------------
// Hardcoded mean parameters
// -----------------------------
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

// -----------------------------
// Simulation parameters (input)
// -----------------------------
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
    double diffusion_factor = 1;

    unsigned long run_seed = 0UL;

    // transport / BTC params
    double Diffusion_coefficient = 0.0;
    double t_end_pdf = 0.0;

    std::vector<double> xLocations;
    enum class correlationmode {exponentialfit, derivative, gaussian, matern} CorrelationModel = correlationmode::derivative;
    pair<double,double> correlation_x_range = {0.001,0.02};
    pair<double,double> correlation_y_range = {0.001,0.02};
};

// -----------------------------
// Outputs (written by runner)
// -----------------------------
struct RunOutputs
{
    std::string run_dir;

    // per-location stacks (each element is a TimeSeriesSet holding all realizations + upscaled)
    std::vector<TimeSeriesSet<double>> Fine_Scale_BTCs;

    // mean BTC per location (one series per x)
    TimeSeriesSet<double> mean_BTCs;

    // upscaled outputs (paths)
    std::string up_dir;
    std::string up_btc_path;
    std::string up_btc_deriv_path;
};

// -----------------------------
// API
// -----------------------------

// Creates/chooses a run directory (rank0 decides; broadcasts to all ranks).
std::string prepare_run_dir_mpi(
    const std::string& output_dir,
    const std::string& resume_run_dir,
    const RunOptions& opts,
    int rank,
    const std::string& run_tag = "");

// Runs fine loop (unless upscale_only) + mean building + upscaled run.
// Requires out.run_dir already set (usually from prepare_run_dir_mpi).
bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    RunOutputs& out,
    int rank);

