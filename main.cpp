// main.cpp

#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"

#include <iostream>
#include <string>
#include <mpi.h>

#include "sim_helpers.h"
#include "plotter.h"
#include "sim_runner.h"

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
    // Options (CLI)
    // -----------------------------
    RunOptions opts;
    std::string resume_run_dir = output_dir + "/run_20260115_132010";

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        if (a == "--upscale-only") opts.upscale_only = true;
        else if (a == "--mean-ts") opts.use_timeseriesset_mean = true;
        else if (a == "--no-mean-ts") opts.use_timeseriesset_mean = false;
        else if (a.rfind("--run-dir=", 0) == 0) resume_run_dir = a.substr(std::string("--run-dir=").size());

        else if (a == "--tbase")      opts.tbase_mode = TBaseMode::Fixed;
        else if (a == "--t-upscaled") opts.tbase_mode = TBaseMode::FromUpscaled;
        else if (a == "--t-fine")     opts.tbase_mode = TBaseMode::FromFirstFine;

        else if (a == "--resample")     opts.align_mode = AlignMode::Resample;
        else if (a == "--make-uniform") opts.align_mode = AlignMode::MakeUniform;

        else if (a == "--plotter") opts.plotter = true;
        else if (a == "--hardcoded-mean") opts.hardcoded_mean = true;
    }

    // hardcoded mean implies upscaled-only
    if (opts.hardcoded_mean) opts.upscale_only = true;

    // -----------------------------
    // Params (keep here)
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

    P.nReal_default = 20;
    P.run_seed = 20260115UL;

    P.xLocations = {0.5, 1.5, 2.5};

    // Hardcoded means (edit here; only used when --hardcoded-mean)
    HardcodedMean H;
    H.lc_mean  = 0.39753;
    H.lx_mean  = 1.25394;
    H.ly_mean  = 0.125017;
    H.dt_mean  = 7.31122e-05;
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
    // Write your “simple” CSV outputs (kept behavior)
    // -----------------------------
    for (int i = 0; i < (int)P.xLocations.size(); ++i) {
        std::string out_cmp_d = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare.csv");
        out.Fine_Scale_BTCs[i].write(out_cmp_d);
        if (rank == 0) std::cout << "Wrote: " << out_cmp_d << "\n";
    }
    {
        std::string out_mean = joinPath(out.run_dir, "BTC_mean.csv");
        out.mean_BTCs.write(out_mean);
        if (rank == 0) std::cout << "Wrote: " << out_mean << "\n";
    }

    // -----------------------------
    // Plotter (rank 0 only)
    // -----------------------------
    if (opts.plotter && rank == 0) {
        const bool okp = run_final_aggregation_and_plots(
            out.run_dir,
            out.up_btc_path,
            out.up_btc_deriv_path,
            opts.tbase_mode,
            opts.align_mode,
            opts.use_timeseriesset_mean,
            /*t_end_cmp=*/10.0,
            /*dt_cmp=*/0.001
        );

        if (!okp) std::cerr << "WARNING: final aggregation/plotting reported failure.\n";
    }

    if (rank == 0) {
        std::cout << "\nMixing PDF simulation complete!\n";
        std::cout << "All outputs saved to: " << out.run_dir << "\n";
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
