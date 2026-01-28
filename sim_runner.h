#pragma once

#include <string>
#include <vector>

#include "TimeSeries.h"
#include <TimeSeriesSet.h>
#include "sim_helpers.h" // joinPath/createDirectory/makeTimestamp/... + TimeSeriesSet
#include "plotter.h"

// Forward declarations to avoid heavy includes here
class Grid2D;

// All run-time options (flags) go here
struct RunOptions
{
    bool upscale_only = false; // loading from fine_scale folders
    bool solve_fine_scale_transport = true;
    bool solve_upscale_transport = true;

    bool use_timeseriesset_mean = true;

    bool plotter = false;          // handled in main (calls plotter module)
    bool hardcoded_mean = true;   // if true: skip fine and skip loading/reconstructing

    // compare time switches (from plotter.h)
    TBaseMode tbase_mode = TBaseMode::Fixed;
    AlignMode align_mode = AlignMode::Resample;
};

// All numerical parameters you want to keep in main
struct SimParams
{
    // domain / grid
    int nx = 300;
    int ny = 100;
    int nu = 100;
    double Lx = 3.0;
    double Ly = 1.0;

    // random field
    double correlation_ls_x = 1.0;
    double correlation_ls_y = 0.1;
    double stdev = 2.0;
    double g_mean = 0.0;

    // transport
    double Diffusion_coefficient = 0.01;
    double t_end_pdf = 20.0;

    // realizations
    int nReal_default = 20;
    unsigned long run_seed = 20260115UL;

    // BTC locations
    std::vector<double> xLocations{0.5, 1.5, 2.5};
};

// Hardcoded mean params (only used when opts.hardcoded_mean = true)
struct HardcodedMean
{
    double lc_mean   = 0.39753;
    double lx_mean   = 1.25394;
    double ly_mean   = 0.125017;
    double dt_mean   = 7.31122e-05;
    double qx_const  = 1.0; // constant qx -> makes a flat inverse CDF
};

// Outputs you might want back in main
struct RunOutputs
{
    std::string run_dir;
    std::string up_dir;

    std::string up_btc_path;
    std::string up_btc_deriv_path;

    // Optional: “wide” per-location stacks for writing CSVs in main (same behavior as your code)
    std::vector<TimeSeriesSet<double>> Fine_Scale_BTCs;
    TimeSeriesSet<double> mean_BTCs; // from Fine_Scale_BTCs mean_ts()
};

// Build or pick run directory and broadcast it
std::string prepare_run_dir_mpi(
    const std::string& output_dir,
    const std::string& resume_run_dir,
    const RunOptions& opts,
    int rank);

// Main entry for running simulation parts.
// - If opts.hardcoded_mean: skip fine AND skip all loading/reconstruction.
// - Else if opts.upscale_only: expects mean files exist (no folder reconstruction inside runner).
// - Else: runs fine loop and then upscaled using computed means.
// Returns true on success.
bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    RunOutputs& out,
    int rank);
