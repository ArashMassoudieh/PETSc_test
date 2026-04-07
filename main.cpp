// File overview: main.cpp is part of the PETSc_test simulation/analysis workflow.
// main.cpp

#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"

#include <iostream>
#include <string>
#include <mpi.h>
#include <cmath>        // std::round, std::abs, pow
#include <limits>
#include <fstream>
#include <algorithm>

#include "sim_helpers.h"
#include "plotter.h"     // TBaseMode, AlignMode, run_final_aggregation_and_plots
#include "sim_runner.h"  // SimParams, RunOptions, HardcodedMean, RunOutputs
#include "Pathway.h"

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

    // Scratch by default (no file loading). Enable via --hardcoded-mean.
    opts.hardcoded_mean = false;

    // NEW: recovery/rebuild mode
    // Read fine/upscaled BTC/PT from files already written in an existing run dir
    // and rebuild BTC_mean... outputs without re-solving.
    opts.read_btc_from_files = true;

    opts.solve_fine_scale_transport = false;
    opts.solve_upscale_transport    = false;

    opts.analyze_qx_ranks = true; // copula ...
    opts.analyze_qx_copula_diagnostics = false; // expensive; opt-in
    opts.qx_copula_bootstrap = 50;
    opts.qx_copula_max_points = 1200;

    opts.perform_particle_tracking = false;
    opts.perform_upscaled_PT = false;

    // -----------------------------
    // NEW: default upscaled mixing model
    // -----------------------------
    opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::ExponentialVDep;

    // -----------------------------
    // NEW: Wiener defaults (OFF unless --wiener)
    // -----------------------------
    opts.wiener_enable  = false;
    opts.wiener_mode    = "2d";         // 1dx | 1dy | 2d
    opts.wiener_Dx      = 0.1;
    opts.wiener_dt      = 1e-3;
    opts.wiener_seed    = 12345UL;
    opts.wiener_release = "left-flux";  // left-uniform | left-flux | center

    // Resume folder: existing source folder (mean, qx, ...)
    bool user_set_qx_cdf        = false;
    bool user_set_run_dir       = false;
    bool user_set_btc_input_dir = false;

    // NEW: BTC input folder for --btc-from-files mode
    std::string btc_input_dir;

    // ---------------------------------------------------------
    // Pick ONE resume folder here (only ONE assignment!)
    //   (You can still override via --run-dir=...)
    // ---------------------------------------------------------
    //std::string resume_run_dir = joinPath(output_dir, "Finished Runs/std=2, D=0, aniso1&0.1");
    //std::string resume_run_dir = joinPath(output_dir, "Finished Runs/100Realizations_std2_D0.01_aniso");
    //std::string resume_run_dir = joinPath(output_dir, "100Realizations_std2_D0.01_aniso");
    //std::string resume_run_dir = joinPath(output_dir, "Finished Runs/100Realizations_20260202_003241_std2_D0.1_aniso1&0.1");
    //std::string resume_run_dir = joinPath(output_dir, "Finished Runs/std=1, D=0, aniso1&0.1");
    //std::string resume_run_dir = joinPath(output_dir, "Finished Runs/100Realizations_20260210_183158_std1_D0.01_aniso1&0.1_df0.15");
    //std::string resume_run_dir = joinPath(output_dir, "Finished Runs/100Realizations_20260211_083055_std1_D0.1_aniso1&0.1_df0.15");
    std::string resume_run_dir = joinPath(output_dir, "Finished Runs/100Realizations_std2_D0.1_aniso_df1_vdep");

    // DEFAULT (your current canonical):
    //resume_run_dir = joinPath(output_dir, "Finished Runs/100Realizations_20260207_111642_std1_D0.1_aniso0.1&1_df0.15");

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
    bool do_calib = false;
    double calib_min = 0.5, calib_max = 3, calib_step = 0.5;

    // NOTE: may be overridden later (in hardcoded_mean mode we default to resume folder)
    std::string black_mean_csv = joinPath(resume_run_dir, "BTC_mean.csv");

    // Score-only mode: DO NOT run simulations; just scan existing run folders and score them.
    bool score_only = false;         // default OFF (enable via --score-only)
    std::string score_root_dir = joinPath(output_dir, "Calibration"); // same as calib_root_dir

    // Calibration root folder (where df-run folders should live)
    std::string calib_root_dir = joinPath(output_dir, "Calibration");

    // -----------------------------
    // CLI
    // -----------------------------
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        // --- run selection ---
        if (a == "--upscale-only") opts.upscale_only = true;
        else if (a == "--hardcoded-mean") opts.hardcoded_mean = true;
        else if (a == "--btc-from-files") {
            opts.read_btc_from_files = true;
        }
        else if (a.rfind("--btc-input-dir=", 0) == 0) {
            btc_input_dir = a.substr(std::string("--btc-input-dir=").size());
            user_set_btc_input_dir = true;
            opts.read_btc_from_files = true;
        }
        else if (a.rfind("--run-dir=", 0) == 0) {
            resume_run_dir = a.substr(std::string("--run-dir=").size());
            user_set_run_dir = true;
        }

        // --- qx inverse-CDF path override ---
        else if (a.rfind("--qx-cdf=", 0) == 0) {
            opts.hardcoded_qx_cdf_path = a.substr(std::string("--qx-cdf=").size());
            user_set_qx_cdf = true;
        }

        // --- transport toggles ---
        else if (a == "--no-fine-transport") opts.solve_fine_scale_transport = false;
        else if (a == "--fine-transport")    opts.solve_fine_scale_transport = true;
        else if (a == "--no-up-transport")   opts.solve_upscale_transport = false;
        else if (a == "--up-transport")      opts.solve_upscale_transport = true;
        else if (a == "--no-qx-ranks")       opts.analyze_qx_ranks = false;
        else if (a == "--qx-ranks")          opts.analyze_qx_ranks = true;
        else if (a == "--copula-analysis")   opts.analyze_qx_copula_diagnostics = true;
        else if (a == "--no-copula-analysis")opts.analyze_qx_copula_diagnostics = false;
        else if (a.rfind("--copula-bootstrap=", 0) == 0) {
            opts.qx_copula_bootstrap = std::max(10, std::stoi(a.substr(std::string("--copula-bootstrap=").size())));
        }
        else if (a.rfind("--copula-max-points=", 0) == 0) {
            opts.qx_copula_max_points = std::max(200, std::stoi(a.substr(std::string("--copula-max-points=").size())));
        }
        else if (a == "--lc-auto")           opts.lc_source = RunOptions::LcSource::AutoPreferRankCopula;
        else if (a == "--lc-raw")            opts.lc_source = RunOptions::LcSource::RawQx;
        else if (a == "--lc-rank")           opts.lc_source = RunOptions::LcSource::RankCopula;

        // --- upscaled mixing model ---
        else if (a == "--mixing-exp") {
            opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::Exponential;
        }
        else if (a == "--mixing-vdep") {
            opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::ExponentialVDep;
        }
        else if (a == "--mixing-gaussian") {
            opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::Gaussian;
        }
        else if (a == "--mixing-matern") {
            opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::Matern;
        }
        else if (a == "--mixing-oneoversum") {
            opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::OneOverSum;
        }

        // --- plotter switches ---
        else if (a == "--plotter") plotter = true;
        else if (a == "--no-plotter") plotter = false;

        else if (a == "--tbase")      tbase_mode = TBaseMode::Fixed;
        else if (a == "--t-upscaled") tbase_mode = TBaseMode::FromUpscaled;
        else if (a == "--t-fine")     tbase_mode = TBaseMode::FromFirstFine;

        else if (a == "--resample")     align_mode = AlignMode::Resample;
        else if (a == "--make-uniform") align_mode = AlignMode::MakeUniform;

        else if (a == "--mean-ts")    use_timeseriesset_mean = true;
        else if (a == "--no-mean-ts") use_timeseriesset_mean = false;

        // --- NEW: Wiener particle diffusion (1D/2D) ---
        else if (a == "--wiener") opts.wiener_enable = true;
        else if (a.rfind("--wmode=", 0) == 0) opts.wiener_mode = a.substr(std::string("--wmode=").size());

        // NOTE: Keeping your exact parsing (even if indices look off) since you said "don't change anything"
        else if (a.rfind("--D=", 0) == 0)     opts.wiener_Dx = std::stod(a.substr(5));
        else if (a.rfind("--rx=", 0) == 0)    opts.wiener_rx = std::stod(a.substr(6));
        else if (a.rfind("--ry=", 0) == 0)    opts.wiener_ry = std::stod(a.substr(7));
        else if (a.rfind("--wdt=", 0) == 0)   opts.wiener_dt = std::stod(a.substr(8));
        else if (a.rfind("--wseed=", 0) == 0) opts.wiener_seed = (unsigned long)std::stoull(a.substr(9));
        else if (a.rfind("--wrelease=", 0) == 0) opts.wiener_release = a.substr(std::string("--wrelease=").size());

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
        else if (a == "--score-only") {
            score_only = true;
            do_calib = true;
        }
        else if (a.rfind("--score-root=", 0) == 0) {
            score_root_dir = a.substr(std::string("--score-root=").size());
            score_only = true;
            do_calib = true;
        }
        else if (a.rfind("--calib-root=", 0) == 0) {
            calib_root_dir = a.substr(std::string("--calib-root=").size());
        }
    }

    // Your rule: hardcoded mean implies upscale-only
    if (opts.hardcoded_mean && !opts.read_btc_from_files) opts.upscale_only = true;

    // NEW: BTC-from-files mode is a rebuild/recovery path, not a solve path
    if (opts.read_btc_from_files) {
        opts.upscale_only = false;
        opts.hardcoded_mean = false;
        opts.solve_fine_scale_transport = false;
        opts.solve_upscale_transport = false;
        opts.perform_particle_tracking = false;
        opts.perform_upscaled_PT = false;
        opts.analyze_qx_ranks = false;

        if (!user_set_btc_input_dir) {
            btc_input_dir = resume_run_dir;
        }

        if (btc_input_dir.empty()) {
            if (rank == 0) {
                std::cerr << "ERROR: --btc-from-files requires an input run folder.\n"
                          << "Use --btc-input-dir=/path/to/existing_run or --run-dir=/path/to/existing_run\n";
            }
            MPI_Abort(PETSC_COMM_WORLD, 301);
        }

        if (!dirExists(btc_input_dir)) {
            if (rank == 0) {
                std::cerr << "ERROR: BTC input directory does not exist:\n"
                          << "  " << btc_input_dir << "\n";
            }
            MPI_Abort(PETSC_COMM_WORLD, 302);
        }

        if (rank == 0) {
            std::cout << "BTC-from-files mode enabled.\n";
            std::cout << "Input/output run folder: " << btc_input_dir << "\n";
        }
    }

    // -----------------------------
    // NEW: mean_ts switch (loaded from main; no extra helpers)
    //   - --mean-ts    => robust mean on TimeSeriesSet (Longest ref / union-safe)
    //   - --no-mean-ts => legacy-style mean (First-series grid)
    // -----------------------------
    opts.mean_ts_mode = use_timeseriesset_mean ? RunOptions::MeanTSMode::Longest
                                              : RunOptions::MeanTSMode::First;

    // -----------------------------
    // Params (kept here)
    // -----------------------------
    SimParams P;
    P.nx = 300; // 300
    P.ny = 100; // 100
    P.nu = 100; // 100
    P.Lx = 3.0;
    P.Ly = 1.0;

    P.correlation_ls_x = 1;
    P.correlation_ls_y = 0.1;

    // calibration target
    P.diffusion_factor = 1;

    // "D" in naming = diffusion coefficient (physics diffusion)
    P.Diffusion_coefficient = 0.1;

    P.stdev = 2.0;
    P.g_mean = 0.0;
    P.CorrelationModel = SimParams::correlationmode::exponentialfit; //or gaussian
    P.correlation_x_range = {0.001, P.correlation_ls_x};
    P.correlation_y_range = {0.001, P.correlation_ls_y};

    P.t_end_pdf = 20.0;

    P.nReal_default = 100; // 100
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

    // -------------------------------------------------
    // RESUME / FILE INPUTS
    //   - Only required when opts.hardcoded_mean == true
    //   - In scratch mode (hardcoded_mean == false): do NOT load anything
    //   - In BTC-from-files mode: also do NOT require qx-cdf / mean_params here
    // -------------------------------------------------
    if (opts.hardcoded_mean && !opts.read_btc_from_files) {

        if (resume_run_dir.empty()) {
            if (rank == 0) {
                std::cerr << "ERROR: resume_run_dir is empty. "
                          << "Set it in code or pass --run-dir=/path/to/resume_folder\n";
            }
            MPI_Abort(PETSC_COMM_WORLD, 66);
        }

        if (!user_set_qx_cdf) {
            opts.hardcoded_qx_cdf_path = joinPath(resume_run_dir, "mean_qx_inverse_cdf.txt");
        }

        if (!fileExists(opts.hardcoded_qx_cdf_path)) {
            if (rank == 0) {
                std::cerr << "ERROR: qx inverse-CDF file not found:\n"
                          << "  " << opts.hardcoded_qx_cdf_path << "\n";
            }
            MPI_Abort(PETSC_COMM_WORLD, 88);
        }

        if (black_mean_csv.empty()) {
            black_mean_csv = joinPath(resume_run_dir, "BTC_mean.csv");
        }

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

    } else {
        // SCRATCH MODE / BTC-FROM-FILES MODE: do not require any qx-cdf file here.
    }

    // -----------------------------
    // Wiener run (early exit)
    // -----------------------------
    if (opts.wiener_enable)
    {
        const std::string run_tag = make_run_tag_std_D_aniso_df(P);
        const std::string wiener_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank, run_tag);

        opts.wiener_rx = P.correlation_ls_x;
        opts.wiener_ry = P.correlation_ls_y;
        opts.wiener_Dx = P.Diffusion_coefficient;
        opts.wiener_dt = 1.0/(opts.wiener_Dx*(1.0/pow(opts.wiener_rx,2) + 1.0/pow(opts.wiener_ry,2)))*0.001;
        return runDiffusionSimulation(opts, 1000000, wiener_dir);
    }

    // -----------------------------
    // HardcodedMean H (only loaded when hardcoded_mean == true)
    // -----------------------------
    HardcodedMean H;

    if (opts.hardcoded_mean && !opts.read_btc_from_files) {

        H.lc_mean  = std::numeric_limits<double>::quiet_NaN();
        H.lx_mean  = std::numeric_limits<double>::quiet_NaN();
        H.ly_mean  = std::numeric_limits<double>::quiet_NaN();
        H.dt_mean  = std::numeric_limits<double>::quiet_NaN();
        H.nu_x     = std::numeric_limits<double>::quiet_NaN();
        H.nu_y     = std::numeric_limits<double>::quiet_NaN();

        H.qx_const = 1.0;

        {
            const std::string mean_path = joinPath(resume_run_dir, "mean_params.txt");
            if (!fileExists(mean_path)) {
                if (rank == 0) {
                    std::cerr << "ERROR: mean_params.txt not found in resume folder:\n"
                              << "  " << mean_path << "\n";
                }
                MPI_Abort(PETSC_COMM_WORLD, 80);
            }

            double lc = 0, lx = 0, ly = 0, dt = 0;
            double nu_x = 0, nu_y = 0;

            if (!read_mean_params_txt(mean_path, lc, lx, ly, dt, nu_x, nu_y)) {
                if (rank == 0) {
                    std::cerr << "ERROR: failed to parse mean_params.txt:\n"
                              << "  " << mean_path << "\n";
                }
                MPI_Abort(PETSC_COMM_WORLD, 81);
            }

            H.lc_mean = lc;
            H.lx_mean = lx;
            H.ly_mean = ly;
            H.dt_mean = dt;
            H.nu_x    = nu_x;
            H.nu_y    = nu_y;

            if (rank == 0) {
                std::cout << "Resume folder (input source): " << resume_run_dir << "\n";
                std::cout << "Loaded mean params: " << mean_path << "\n";
                std::cout << "Using qx inverse-CDF: " << opts.hardcoded_qx_cdf_path << "\n";
            }
        }

    } else {
        // Not needed for scratch or BTC-from-files mode.
        H.lc_mean  = std::numeric_limits<double>::quiet_NaN();
        H.lx_mean  = std::numeric_limits<double>::quiet_NaN();
        H.ly_mean  = std::numeric_limits<double>::quiet_NaN();
        H.dt_mean  = std::numeric_limits<double>::quiet_NaN();
        H.nu_x     = std::numeric_limits<double>::quiet_NaN();
        H.nu_y     = std::numeric_limits<double>::quiet_NaN();
        H.qx_const = 1.0;
    }

    // -----------------------------
    // CALIBRATION / SCORING MODE
    // -----------------------------
    if (opts.read_btc_from_files && do_calib) {
        if (rank == 0) {
            std::cerr << "ERROR: --btc-from-files is currently intended for single-run recovery/rebuild mode, "
                      << "not calibration/scoring mode.\n";
        }
        MPI_Abort(PETSC_COMM_WORLD, 303);
    }

    if (do_calib) {

        if (score_root_dir.empty()) score_root_dir = calib_root_dir;

        // SCORE-ONLY
        if (score_only) {
            if (rank == 0) {
                std::cout << "\n=== SCORE-ONLY calibration ===\n";
                std::cout << "Score root: " << score_root_dir << "\n";
                std::cout << "Black mean: " << black_mean_csv << "\n\n";
            }

            std::vector<double> dfs;
            std::vector<std::string> dirs;

            if (rank == 0) {
                if (!list_calibration_runs_with_df(score_root_dir, dfs, dirs)) {
                    std::cerr << "ERROR: no calibration run folders with df token were found in: "
                              << score_root_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 201);
                }
            }

            if (rank == 0) {
                const std::string calib_csv = joinPath(score_root_dir, "calibration_summary_scored.csv");
                std::ofstream f(calib_csv);
                f << "diffusion_factor,score_logrmse_mean,run_dir\n";

                double best_df = std::numeric_limits<double>::quiet_NaN();
                double best_score = std::numeric_limits<double>::infinity();
                std::string best_run;

                for (size_t i = 0; i < dfs.size(); ++i) {
                    const double dfv = dfs[i];
                    const std::string& rdir = dirs[i];

                    const double score = score_upscaled_vs_black_mean_from_compare(
                        black_mean_csv, rdir, P.xLocations);

                    f << dfv << "," << score << "," << rdir << "\n";
                    std::cout << "df=" << dfv << "  score=" << score << "  run=" << rdir << "\n";

                    if (std::isfinite(score) && score < best_score) {
                        best_score = score;
                        best_df = dfv;
                        best_run = rdir;
                    }
                }

                std::cout << "\n=== BEST diffusion_factor (score-only) ===\n";
                std::cout << "best_df     = " << best_df << "\n";
                std::cout << "best_score  = " << best_score << "\n";
                std::cout << "best_run    = " << best_run << "\n";
                std::cout << "Summary CSV = " << calib_csv << "\n";
            }

            MPI_Barrier(PETSC_COMM_WORLD);
            return 0;
        }

        // RUN-THEN-SCORE
        if (rank == 0) {
            std::cout << "\n=== CALIBRATING diffusion_factor (run-then-score) ===\n";
            std::cout << "Range: [" << calib_min << ", " << calib_max << "] step " << calib_step << "\n";
            std::cout << "Black mean file: " << black_mean_csv << "\n";
            std::cout << "Runs will be saved under: " << calib_root_dir << "\n";
            std::cout << "Scoring uses per-x BTC_Compare: <run_dir>/x=0.50BTC_Compare.csv etc\n\n";
        }

        createDirectory(output_dir);
        createDirectory(calib_root_dir);

        std::vector<double> df_list_rank0;
        std::vector<std::string> run_dirs_rank0;

        for (double df = calib_min; df <= calib_max + 0.5 * calib_step; df += calib_step) {
            P.diffusion_factor = df;

            const std::string run_tag = make_run_tag_std_D_aniso_df(P);

            RunOutputs out;
            out.run_dir = prepare_run_dir_mpi(calib_root_dir, resume_run_dir, opts, rank, run_tag);

            const bool ok = run_simulation_blocks(P, opts, H, resume_run_dir, out, rank);
            if (!ok) {
                if (rank == 0) std::cerr << "ERROR: simulation failed for df=" << df << "\n";
                MPI_Abort(PETSC_COMM_WORLD, 123);
            }

            if (rank == 0) {
                for (int i = 0; i < (int)P.xLocations.size(); ++i) {
                    const std::string out_cmp_legacy = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare.csv");
                    const std::string out_cmp_pdf    = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare_pdf.csv");
                    const std::string out_cmp_cdf    = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare_cdf.csv");

                    // legacy transport compare -> PDF
                    out.Fine_Scale_BTCs[i].write(out_cmp_legacy);

                    // explicit transport compare
                    out.Fine_Scale_BTCs_pdf[i].write(out_cmp_pdf);
                    out.Fine_Scale_BTCs_cdf[i].write(out_cmp_cdf);
                }

                // legacy transport mean -> PDF
                {
                    const std::string out_mean = joinPath(out.run_dir, "BTC_mean.csv");
                    out.mean_BTCs.write(out_mean);
                }

                // explicit transport means
                {
                    const std::string out_mean_pdf = joinPath(out.run_dir, "BTC_mean_pdf.csv");
                    out.mean_transport_pdf.write(out_mean_pdf);
                }
                {
                    const std::string out_mean_cdf = joinPath(out.run_dir, "BTC_mean_cdf.csv");
                    out.mean_transport_cdf.write(out_mean_cdf);
                }

                // keep existing separated transport alias output
                {
                    const std::string out_tr = joinPath(out.run_dir, "BTC_mean_transport_full.csv");
                    out.mean_transport_full.write(out_tr);
                }

                {
                    const std::string out_pt_pdf = joinPath(out.run_dir, "PT_mean_pdf.csv");
                    out.mean_pt_pdf.write(out_pt_pdf);
                }
                {
                    const std::string out_pt_cdf = joinPath(out.run_dir, "PT_mean_cdf.csv");
                    out.mean_pt_cdf.write(out_pt_cdf);
                }

                df_list_rank0.push_back(df);
                run_dirs_rank0.push_back(out.run_dir);

                std::cout << "Finished df=" << df << "  run=" << out.run_dir << "\n";
            }

            MPI_Barrier(PETSC_COMM_WORLD);
        }

        if (rank == 0) {
            const std::string calib_csv = joinPath(calib_root_dir, "calibration_summary.csv");
            std::ofstream f(calib_csv);
            f << "diffusion_factor,score_logrmse_mean,run_dir\n";

            double best_df = std::numeric_limits<double>::quiet_NaN();
            double best_score = std::numeric_limits<double>::infinity();
            std::string best_run_dir;

            for (size_t i = 0; i < df_list_rank0.size(); ++i) {
                const double dfv = df_list_rank0[i];
                const std::string& rdir = run_dirs_rank0[i];

                const double score = score_upscaled_vs_black_mean_from_compare(
                    black_mean_csv, rdir, P.xLocations);

                f << dfv << "," << score << "," << rdir << "\n";
                std::cout << "df=" << dfv << "  score=" << score << "  run=" << rdir << "\n";

                if (std::isfinite(score) && score < best_score) {
                    best_score = score;
                    best_df = dfv;
                    best_run_dir = rdir;
                }
            }

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

    if (opts.read_btc_from_files) {
        // Rebuild means directly inside the existing run folder.
        out.run_dir = btc_input_dir;
        if (!out.run_dir.empty() && out.run_dir.back() != '/' && out.run_dir.back() != '\\')
            out.run_dir += "/";
    } else {
        out.run_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank, run_tag);
    }

    const std::string input_source_dir = opts.read_btc_from_files ? btc_input_dir : resume_run_dir;

    const bool ok = run_simulation_blocks(P, opts, H, input_source_dir, out, rank);
    if (!ok) {
        if (rank == 0) std::cerr << "ERROR: simulation runner failed.\n";
        MPI_Abort(PETSC_COMM_WORLD, 123);
    }

    // Write compare + mean CSVs (RANK 0 ONLY — avoid clobber/races)
    if (rank == 0) {
        for (int i = 0; i < (int)P.xLocations.size(); ++i) {
            const std::string out_cmp_legacy = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare.csv");
            const std::string out_cmp_pdf    = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare_pdf.csv");
            const std::string out_cmp_cdf    = joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "BTC_Compare_cdf.csv");

            // legacy transport compare -> PDF
            out.Fine_Scale_BTCs[i].write(out_cmp_legacy);
            std::cout << "Wrote: " << out_cmp_legacy << "\n";

            // explicit transport compare files
            out.Fine_Scale_BTCs_pdf[i].write(out_cmp_pdf);
            std::cout << "Wrote: " << out_cmp_pdf << "\n";

            out.Fine_Scale_BTCs_cdf[i].write(out_cmp_cdf);
            std::cout << "Wrote: " << out_cmp_cdf << "\n";
        }
        {
            // legacy transport mean -> PDF
            const std::string out_mean = joinPath(out.run_dir, "BTC_mean.csv");
            out.mean_BTCs.write(out_mean);
            std::cout << "Wrote: " << out_mean << "\n";
        }

        // explicit transport means
        {
            const std::string out_mean_pdf = joinPath(out.run_dir, "BTC_mean_pdf.csv");
            out.mean_transport_pdf.write(out_mean_pdf);
            std::cout << "Wrote: " << out_mean_pdf << "\n";
        }
        {
            const std::string out_mean_cdf = joinPath(out.run_dir, "BTC_mean_cdf.csv");
            out.mean_transport_cdf.write(out_mean_cdf);
            std::cout << "Wrote: " << out_mean_cdf << "\n";
        }

        // keep separated transport alias output
        {
            const std::string out_tr = joinPath(out.run_dir, "BTC_mean_transport_full.csv");
            out.mean_transport_full.write(out_tr);
            std::cout << "Wrote: " << out_tr << "\n";
        }
        {
            const std::string out_pt_pdf = joinPath(out.run_dir, "PT_mean_pdf.csv");
            out.mean_pt_pdf.write(out_pt_pdf);
            std::cout << "Wrote: " << out_pt_pdf << "\n";
        }
        {
            const std::string out_pt_cdf = joinPath(out.run_dir, "PT_mean_cdf.csv");
            out.mean_pt_cdf.write(out_pt_cdf);
            std::cout << "Wrote: " << out_pt_cdf << "\n";
        }

        // Plotter (rank 0 only) — AFTER writing outputs
        if (plotter) {
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

        if (opts.read_btc_from_files) {
            std::cout << "\nBTC-from-files rebuild complete!\n";
            std::cout << "Source folder: " << btc_input_dir << "\n";
            std::cout << "Rebuilt outputs written in-place to: " << out.run_dir << "\n";
        } else {
            std::cout << "\nMixing PDF simulation complete!\n";
            std::cout << "Resume folder (input source): " << resume_run_dir << "\n";
            std::cout << "All outputs saved to (new run dir): " << out.run_dir << "\n";
        }

        if (opts.wiener_enable) {
            std::cout << "Wiener enabled: mode=" << opts.wiener_mode
                      << " Dx=" << opts.wiener_Dx
                      << " dt=" << opts.wiener_dt
                      << " release=" << opts.wiener_release
                      << " seed=" << opts.wiener_seed << "\n";
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
