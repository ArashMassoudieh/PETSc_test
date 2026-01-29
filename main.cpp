// main.cpp

#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"

#include <iostream>
#include <string>
#include <mpi.h>

#include "sim_helpers.h"
#include "plotter.h"     // TBaseMode, AlignMode, run_final_aggregation_and_plots
#include "sim_runner.h"  // SimParams, RunOptions, HardcodedMean, RunOutputs

int main(int argc, char** argv)
{
    PETScInit petsc(argc, argv);

#ifdef Arash
    std::string output_dir = "/home/arash/Projects/UpscalingResults";
#elif PowerEdge
    std::string output_dir = "/mnt/3rd900/Projects/PETSc_test/Results";
#elif Behzad
    std::string output_dir = "/home/behzad/Projects/PETSc_test/Results";
#elif SligoCreek
    std::string output_dir = "/media/arash/E/Projects/PETSc_test/Results";
#elif WSL
    std::string output_dir = "/home/behzad/Projects/PETSc_test/Results";
#else
    std::string output_dir = "./Results";
#endif

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // -----------------------------
    // Simulation options (passed to sim_runner)
    // -----------------------------
    RunOptions opts;
    opts.upscale_only   = false;
    opts.hardcoded_mean = false; // IMPORTANT: do NOT default true (otherwise you always skip fine loop)
    opts.solve_fine_scale_transport = true;
    opts.solve_upscale_transport    = true;

    // If you want a default resume folder for upscale-only:
    std::string resume_run_dir = joinPath(output_dir, "run_20260129_093954_params");

    // If you want to hardcode the qx inverse-CDF path from main (OPTIONAL):
    // (This file exists in the run folder you showed.)
     opts.hardcoded_qx_cdf_path = joinPath(resume_run_dir, "mean_qx_inverse_cdf.txt");

    // -----------------------------
    // Plot options (kept in main only)
    // -----------------------------
    bool plotter = false;
    TBaseMode tbase_mode = TBaseMode::Fixed;
    AlignMode align_mode = AlignMode::MakeUniform;
    bool use_timeseriesset_mean = true;

    // -----------------------------
    // CLI
    // -----------------------------
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        // --- run selection ---
        if (a == "--upscale-only") opts.upscale_only = true;
        else if (a == "--hardcoded-mean") opts.hardcoded_mean = true;
        else if (a.rfind("--run-dir=", 0) == 0) resume_run_dir = a.substr(std::string("--run-dir=").size());

        // --- hardcoded qx inverse-CDF path (u,v csv; header optional) ---
        else if (a.rfind("--qx-cdf=", 0) == 0) opts.hardcoded_qx_cdf_path = a.substr(std::string("--qx-cdf=").size());

        // --- transport toggles ---
        else if (a == "--no-fine-transport") opts.solve_fine_scale_transport = false;
        else if (a == "--fine-transport")    opts.solve_fine_scale_transport = true;
        else if (a == "--no-up-transport")   opts.solve_upscale_transport = false;
        else if (a == "--up-transport")      opts.solve_upscale_transport = true;

        // --- plotter switches ---
        else if (a == "--plotter") plotter = true;

        else if (a == "--tbase")      tbase_mode = TBaseMode::Fixed;
        else if (a == "--t-upscaled") tbase_mode = TBaseMode::FromUpscaled;
        else if (a == "--t-fine")     tbase_mode = TBaseMode::FromFirstFine;

        else if (a == "--resample")     align_mode = AlignMode::Resample;
        else if (a == "--make-uniform") align_mode = AlignMode::MakeUniform;

        else if (a == "--mean-ts")    use_timeseriesset_mean = true;
        else if (a == "--no-mean-ts") use_timeseriesset_mean = false;
    }

    // Your rule: hardcoded mean implies upscale-only
    if (opts.hardcoded_mean) opts.upscale_only = true;

    // -----------------------------
    // Params (kept here)
    // -----------------------------
    SimParams P;
    P.nx = 300;
    P.ny = 100;
    P.nu = 100;
    P.Lx = 3.0;
    P.Ly = 1.0;

    P.correlation_ls_x = 1.0;
    P.correlation_ls_y = 0.1;
    P.stdev = 2.0;
    P.g_mean = 0.0;

    P.Diffusion_coefficient = 0.01;
    P.t_end_pdf = 20.0;

    P.nReal_default = 1;
    P.run_seed = 20260115UL;

    P.xLocations = {0.5, 1.5, 2.5};

    // Hardcoded means (used only when --hardcoded-mean)
    HardcodedMean H;
    H.lc_mean  = 0.39753;
    H.lx_mean  = 1.25394;
    H.ly_mean  = 0.125017;
    H.dt_mean  = 7.31122e-05;

    // Optional overwrite from file (if present)
    {
        const std::string mean_path = joinPath(resume_run_dir, "mean_params.txt");
        if (fileExists(mean_path)) {
            if (load_hardcoded_mean_from_file(mean_path, H) && rank == 0) {
                std::cout << "Loaded mean params: " << mean_path << "\n";
            }
        }
    }

    // fallback only if --qx-cdf not provided or load fails
    H.qx_const = 1.0;

    // -----------------------------
    // Prepare run_dir (MPI-safe)
    // -----------------------------
    RunOutputs out;
    out.run_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank);

    // -----------------------------
    // Run simulation blocks (fine loop + mean + upscaled)
    // -----------------------------
    const bool ok = run_simulation_blocks(P, opts, H, out, rank);
    if (!ok) {
        if (rank == 0) std::cerr << "ERROR: simulation runner failed.\n";
        MPI_Abort(PETSC_COMM_WORLD, 123);
    }

    // -----------------------------
    // Write compare + mean CSVs
    // -----------------------------
    for (int i = 0; i < (int)P.xLocations.size(); ++i) {
        const std::string out_cmp = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare.csv");
        out.Fine_Scale_BTCs[i].write(out_cmp);
        if (rank == 0) std::cout << "Wrote: " << out_cmp << "\n";
    }
    {
        const std::string out_mean = joinPath(out.run_dir, "BTC_mean.csv");
        out.mean_BTCs.write(out_mean);
        if (rank == 0) std::cout << "Wrote: " << out_mean << "\n";
    }

    // -----------------------------
    // Plotter (rank 0 only)
    // -----------------------------
    if (plotter && rank == 0) {
        const bool okp = run_final_aggregation_and_plots(
            out.run_dir,
            out.up_btc_path,
            out.up_btc_deriv_path,
            tbase_mode,
            align_mode,
            use_timeseriesset_mean,
            /*t_end_cmp=*/10.0,
            /*dt_cmp=*/0.001
        );

        if (!okp) std::cerr << "WARNING: final aggregation/plotting reported failure.\n";
    }

    if (rank == 0) {
        std::cout << "\nMixing PDF simulation complete!\n";
        std::cout << "All outputs saved to: " << out.run_dir << "\n";

        if (opts.hardcoded_mean) {
            if (!opts.hardcoded_qx_cdf_path.empty()) {
                std::cout << "Hardcoded mean mode: qx inverse-CDF path = "
                          << opts.hardcoded_qx_cdf_path << "\n";
            } else {
                std::cout << "Hardcoded mean mode: no --qx-cdf provided; using constant qx fallback.\n";
            }
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
