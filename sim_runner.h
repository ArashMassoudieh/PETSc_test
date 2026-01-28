// sim_runner.h

#pragma once

#include <string>
#include <vector>

// Keep sim_runner.h LIGHT.
// Do NOT include grid.h here (it causes heavy include chains / incomplete type issues).
// Only include what we must for types used in signatures.
#include "TimeSeries.h"     // for TimeSeries<T> and TimeSeriesSet<T>
#include "sim_helpers.h"    // fileExists, joinPath, createDirectory, read_mean_* helpers, mean_of, fmt_x, etc.

// -----------------------------
// Simulation switches / options
// -----------------------------
struct RunOptions
{
    // run modes
    bool upscale_only = false;              // skip fine loop; do only upscaled (requires run_dir unless hardcoded_mean)
    bool hardcoded_mean = true;            // use hardcoded lc/lx/ly/dt; and optionally load qx CDF from file

    // transport toggles
    bool solve_fine_scale_transport = true;
    bool solve_upscale_transport    = true;

    // NEW: when hardcoded_mean=true, you can supply a qx inverse-CDF file (u,v).
    // If empty or load fails -> fallback to constant qx_const.
    std::string hardcoded_qx_cdf_path;
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

    // fallback only when hardcoded_qx_cdf_path is not provided or load fails
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

    int nReal_default = 1;

    // random field params
    double correlation_ls_x = 0.0;
    double correlation_ls_y = 0.0;
    double stdev = 0.0;
    double g_mean = 0.0;

    unsigned long run_seed = 0UL;

    // transport / BTC params
    double Diffusion_coefficient = 0.0;
    double t_end_pdf = 0.0;

    std::vector<double> xLocations;
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
    int rank);

// Runs fine loop (unless upscale_only) + mean building + upscaled run.
// Requires out.run_dir already set (usually from prepare_run_dir_mpi).
bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    RunOutputs& out,
    int rank);
