// main.cpp

#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"

#include <iostream>
#include <string>
#include <mpi.h>
#include <cmath>        // std::round, std::abs

#include "sim_helpers.h"
#include "plotter.h"     // TBaseMode, AlignMode, run_final_aggregation_and_plots
#include "sim_runner.h"  // SimParams, RunOptions, HardcodedMean, RunOutputs

int main(int argc, char** argv)
{
    PETScInit petsc(argc, argv);

#ifdef Arash
    std::string output_dir = "/home/arash/Projects/UpscalingResults";
#elif PowerEdge
    std::string output_dir = "/mnt/3rd900/Projects/PETSc_test/UpscalingResults";
#elif Behzad
    std::string output_dir = "/home/behzad/Projects/PETSc_test/UpscalingResults";
#elif SligoCreek
    std::string output_dir = "/media/arash/E/Projects/PETSc_test/UpscalingResults";
#elif Jason
    std::string output_dir = "/home/arash/Projects/UpscalingResults";
#elif WSL
    std::string output_dir = "/home/behzad/Projects/PETSc_test/UpscalingResults";
#else
    std::string output_dir = "./UpscalingResults";
#endif

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // -----------------------------
    // Simulation options (passed to sim_runner)
    // -----------------------------
    RunOptions opts;

    opts.upscale_only   = false;
    opts.hardcoded_mean = false; // (your current setting)
    opts.solve_fine_scale_transport = false;
    opts.solve_upscale_transport    = false;

    // Resume folder: existing source folder (mean, qx, ...)
    bool user_set_qx_cdf  = false;
    bool user_set_run_dir = false;

    // Your existing resume folder naming
    std::string resume_run_dir = joinPath(output_dir, "run_20260131_112925_std2_D0.1_aniso");

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
        else if (a.rfind("--run-dir=", 0) == 0) {
            resume_run_dir = a.substr(std::string("--run-dir=").size());
            user_set_run_dir = true;
        }

        // --- hardcoded qx inverse-CDF path (u,v csv; header optional) ---
        else if (a.rfind("--qx-cdf=", 0) == 0) {
            opts.hardcoded_qx_cdf_path = a.substr(std::string("--qx-cdf=").size());
            user_set_qx_cdf = true;
        }

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

    P.correlation_ls_x = 1;
    P.correlation_ls_y = 0.1;
    P.diffusion_factor = 0.15; // Calibration coefficient
    P.stdev = 2.0;
    P.g_mean = 0.0;
    P.CorrelationModel = SimParams::correlationmode::exponentialfit;
    P.correlation_x_range = {0.001,0.3};
    P.correlation_y_range = {0.001,0.3};
    // "D" in your naming = diffusion coefficient
    P.Diffusion_coefficient = 0.01;

    P.t_end_pdf = 20.0;

    P.nReal_default = 20;
    P.run_seed = 20260115UL;

    P.xLocations = {0.5, 1.5, 2.5};

    // -----------------------------
    // stdev integer check (as before)
    // -----------------------------
    const int std_int = static_cast<int>(std::round(P.stdev));
    if (std::abs(P.stdev - std_int) > 1e-12) {
        if (rank == 0) {
            std::cerr << "ERROR: P.stdev must be an integer value (e.g., 1 or 2). Got "
                      << P.stdev << "\n";
        }
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Resume folder should not be empty in your workflow
    if (resume_run_dir.empty()) {
        if (rank == 0) {
            std::cerr << "ERROR: resume_run_dir is empty. "
                      << "Set it in code or pass --run-dir=/path/to/resume_folder\n";
        }
        MPI_Abort(PETSC_COMM_WORLD, 66);
    }

    // If user did NOT provide --qx-cdf, point to the expected mean file in resume_run_dir.
    if (!user_set_qx_cdf) {
        opts.hardcoded_qx_cdf_path = joinPath(resume_run_dir, "mean_qx_inverse_cdf.txt");
    }

    // -----------------------------
    // WARNING-ONLY sanity check:
    // resume folder is an *input source* (mean, qx, ...),
    // run folder is new and based on P.
    // So mismatch is not fatal: print warning only.
    // Only do this warning when resume_run_dir was not explicitly overridden.
    // -----------------------------
    if (opts.upscale_only || opts.hardcoded_mean) {
        std::string err;
        if (!validate_resume_run_dir(resume_run_dir, P, err)) {
            if (rank == 0) {
                std::cerr << "\nWARNING: resume_run_dir name does not match current parameters.\n"
                          << err
                          << "Continuing anyway (resume folder used as input source).\n\n";
            }
        }
    }

    // -----------------------------
    // Build run tag for NEW run_dir folder name
    // Example: std2_D1e-3_aniso
    // -----------------------------
    const std::string run_tag = make_run_tag_std_D_aniso(P);

    // -----------------------------
    // Hardcoded means (used only when --hardcoded-mean)
    // -----------------------------
    HardcodedMean H;
    H.lc_mean  = 0.396413;
    H.lx_mean  = 1.24365;
    H.ly_mean  = 0.124846;
    H.dt_mean  = 7.31122e-05;
    H.nu_x = 1.5;
    H.nu_y = 1.5;

    // Optional overwrite from mean_params.txt (if present)
    {
        const std::string mean_path = joinPath(resume_run_dir, "mean_params.txt");
        if (fileExists(mean_path)) {
            double lc = H.lc_mean, lx = H.lx_mean, ly = H.ly_mean, dt = H.dt_mean; double nu_x=H.nu_x; double nu_y=H.nu_y;
            if (read_mean_params_txt(mean_path, lc, lx, ly, dt, nu_x, nu_y)) {
                H.lc_mean = lc;
                H.lx_mean = lx;
                H.ly_mean = ly;
                H.dt_mean = dt;
                H.nu_x = nu_x;
                H.nu_y = nu_y;

                if (rank == 0) std::cout << "Loaded mean params: " << mean_path << "\n";
            }
        }
    }

    // constant fallback if no mean qx could be loaded / computed
    H.qx_const = 1.0;

    // -----------------------------
    // Prepare run_dir (MPI-safe)
    // NOTE: This creates/chooses the OUTPUT run folder.
    // It must NOT overwrite resume_run_dir.
    // It should append run_tag only when creating a new run directory.
    // -----------------------------
    RunOutputs out;
    out.run_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank, run_tag);

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
        std::cout << "Resume folder (input source): " << resume_run_dir << "\n";
        std::cout << "All outputs saved to (new run dir): " << out.run_dir << "\n";

        if (opts.hardcoded_mean) {
            if (!opts.hardcoded_qx_cdf_path.empty()) {
                std::cout << "Hardcoded mean mode: qx inverse-CDF path (preferred mean file) = "
                          << opts.hardcoded_qx_cdf_path << "\n";
                std::cout << "If missing, runner will compute from qx_inverse_cdfs.txt or use constant fallback.\n";
            } else {
                std::cout << "Hardcoded mean mode: no qx path set; using constant qx fallback.\n";
            }
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
