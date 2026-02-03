// main.cpp

#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"

#include <iostream>
#include <string>
#include <mpi.h>
#include <cmath>        // std::round, std::abs
#include <limits>

#include "sim_helpers.h"
#include "plotter.h"     // TBaseMode, AlignMode, run_final_aggregation_and_plots
#include "sim_runner.h"  // SimParams, RunOptions, HardcodedMean, RunOutputs

static bool parse_range3(const std::string& s, double& a, double& b, double& c)
{
    // expects "min:max:step"
    auto p1 = s.find(':');
    if (p1 == std::string::npos) return false;
    auto p2 = s.find(':', p1 + 1);
    if (p2 == std::string::npos) return false;

    std::string s1 = s.substr(0, p1);
    std::string s2 = s.substr(p1 + 1, p2 - (p1 + 1));
    std::string s3 = s.substr(p2 + 1);

    return try_parse_double(s1, a) && try_parse_double(s2, b) && try_parse_double(s3, c);
}

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
    std::string output_dir = "/home/arash/Projects/PETSc_test/UpscalingResults";
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
    opts.hardcoded_mean = true; // your current setting (keeps calibration fast)
    opts.solve_fine_scale_transport = false;
    opts.solve_upscale_transport    = true;

    // Resume folder: existing source folder (mean, qx, ...)
    bool user_set_qx_cdf  = false;
    bool user_set_run_dir = false;

    // Existing resume folder naming (input source)
    std::string resume_run_dir = joinPath(output_dir, "100Realizations_20260202_003241_std2_D0.1_aniso");

    // -----------------------------
    // Plot options (kept in main only)
    // -----------------------------
    bool plotter = false;
    TBaseMode tbase_mode = TBaseMode::Fixed;
    AlignMode align_mode = AlignMode::MakeUniform;
    bool use_timeseriesset_mean = true;

    // -----------------------------
    // Calibration options
    // -----------------------------
    bool do_calib = true;
    double calib_min = 0.1, calib_max = 0.5, calib_step = 0.05;
    std::string black_mean_csv = joinPath(resume_run_dir, "BTC_mean.csv"); // can be overridden

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

        // --- calibration ---
        else if (a == "--calib-df") {
            do_calib = true;
        }
        else if (a.rfind("--calib-df=", 0) == 0) {
            do_calib = true;
            const std::string spec = a.substr(std::string("--calib-df=").size());
            if (!parse_range3(spec, calib_min, calib_max, calib_step)) {
                if (rank == 0) std::cerr << "ERROR: --calib-df=min:max:step expected. Got: " << spec << "\n";
                MPI_Abort(PETSC_COMM_WORLD, 77);
            }
        }
        else if (a.rfind("--black-mean=", 0) == 0) {
            black_mean_csv = a.substr(std::string("--black-mean=").size());
        }
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

    // THIS is the calibration target (we will sweep it if do_calib)
    P.diffusion_factor = 0.15; // Calibration coefficient

    // "D" in naming = diffusion coefficient (physics diffusion)
    P.Diffusion_coefficient = 0.1;

    P.stdev = 2.0;
    P.g_mean = 0.0;
    P.CorrelationModel = SimParams::correlationmode::exponentialfit;
    P.correlation_x_range = {0.001, P.correlation_ls_x};
    P.correlation_y_range = {0.001, P.correlation_ls_y};

    P.t_end_pdf = 20.0;

    P.nReal_default = 100;
    P.run_seed = 20260115UL;

    P.xLocations = {0.5, 1.5, 2.5};

    // -----------------------------
    // stdev integer check
    // -----------------------------
    const int std_int = static_cast<int>(std::round(P.stdev));
    if (std::abs(P.stdev - std_int) > 1e-12) {
        if (rank == 0) {
            std::cerr << "ERROR: P.stdev must be an integer value (e.g., 1 or 2). Got "
                      << P.stdev << "\n";
        }
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Resume folder should not be empty
    if (resume_run_dir.empty()) {
        if (rank == 0) {
            std::cerr << "ERROR: resume_run_dir is empty. "
                      << "Set it in code or pass --run-dir=/path/to/resume_folder\n";
        }
        MPI_Abort(PETSC_COMM_WORLD, 66);
    }

    // If user did NOT provide --qx-cdf, point to expected mean file in resume_run_dir
    if (!user_set_qx_cdf) {
        opts.hardcoded_qx_cdf_path = joinPath(resume_run_dir, "mean_qx_inverse_cdf.txt");
    }

    // warning-only sanity check
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
    // Hardcoded means (used only when --hardcoded-mean)
    // -----------------------------
    HardcodedMean H;
    H.lc_mean  = 0.396413;
    H.lx_mean  = 1.24365;
    H.ly_mean  = 0.124846;
    H.dt_mean  = 7.31122e-05;
    H.nu_x = 1.5;
    H.nu_y = 1.5;
    H.qx_const = 1.0;

    // Optional overwrite from mean_params.txt (if present)
    {
        const std::string mean_path = joinPath(resume_run_dir, "mean_params.txt");
        if (fileExists(mean_path)) {
            double lc = H.lc_mean, lx = H.lx_mean, ly = H.ly_mean, dt = H.dt_mean;
            double nu_x = H.nu_x, nu_y = H.nu_y;
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

    // -----------------------------
    // CALIBRATION MODE
    // -----------------------------
    if (do_calib) {
        if (rank == 0) {
            std::cout << "\n=== CALIBRATING diffusion_factor ===\n";
            std::cout << "Range: [" << calib_min << ", " << calib_max << "] step " << calib_step << "\n";
            std::cout << "Black mean file: " << black_mean_csv << "\n";
            std::cout << "Using red curves from: x=*.BTC_Compare.csv (Upscaled column)\n\n";
        }

        createDirectory(output_dir);

        const std::string calib_csv = joinPath(output_dir, "calibration_summary.csv");
        if (rank == 0) {
            std::ofstream f(calib_csv);
            f << "diffusion_factor,score_rmse_mean,run_dir\n";
        }

        double best_df = std::numeric_limits<double>::quiet_NaN();
        double best_score = std::numeric_limits<double>::infinity();
        std::string best_run_dir;

        // sweep
        for (double df = calib_min; df <= calib_max + 0.5*calib_step; df += calib_step) {
            P.diffusion_factor = df;

            // build run tag including df
            const std::string run_tag = make_run_tag_std_D_aniso_df(P);

            RunOutputs out;
            out.run_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank, run_tag);

            const bool ok = run_simulation_blocks(P, opts, H, out, rank);
            if (!ok) {
                if (rank == 0) std::cerr << "ERROR: simulation failed for df=" << df << "\n";
                MPI_Abort(PETSC_COMM_WORLD, 123);
            }

            // score using RED from BTC_Compare + BLACK from BTC_mean
            double score = std::numeric_limits<double>::infinity();
            if (rank == 0) {
                score = score_upscaled_vs_black_mean_from_compare(black_mean_csv, out.run_dir);

                std::ofstream f(calib_csv, std::ios::app);
                f << df << "," << score << "," << out.run_dir << "\n";

                std::cout << "df=" << df << "  score=" << score << "  run=" << out.run_dir << "\n";

                if (std::isfinite(score) && score < best_score) {
                    best_score = score;
                    best_df = df;
                    best_run_dir = out.run_dir;
                }
            }

            MPI_Barrier(PETSC_COMM_WORLD);
        }

        if (rank == 0) {
            std::cout << "\n=== BEST diffusion_factor ===\n";
            std::cout << "best_df     = " << best_df << "\n";
            std::cout << "best_score  = " << best_score << "\n";
            std::cout << "best_run    = " << best_run_dir << "\n";
            std::cout << "Summary CSV = " << calib_csv << "\n";
        }

        MPI_Barrier(PETSC_COMM_WORLD);
        return 0;
    }

    // -----------------------------
    // NORMAL MODE (single run)
    // -----------------------------
    const std::string run_tag = make_run_tag_std_D_aniso_df(P);

    RunOutputs out;
    out.run_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank, run_tag);

    const bool ok = run_simulation_blocks(P, opts, H, out, rank);
    if (!ok) {
        if (rank == 0) std::cerr << "ERROR: simulation runner failed.\n";
        MPI_Abort(PETSC_COMM_WORLD, 123);
    }

    // Write compare + mean CSVs
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

    // Plotter (rank 0 only)
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
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
