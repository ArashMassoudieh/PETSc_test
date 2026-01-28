// main.cpp
#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <cstdlib>
#include <algorithm>

#include "grid.h"
#include "TimeSeries.h"
#include "Pathway.h"
#include "PathwaySet.h"

#include "sim_helpers.h"
#include "plotter.h"

int main(int argc, char** argv) {
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
    // CLI options
    // -----------------------------
    bool upscale_only = false;
    bool solve_fine_scale_transport = true;
    bool solve_upscale_transport = true;

    std::string resume_run_dir = output_dir + "/run_20260115_132010";

    bool use_timeseriesset_mean = true;

    // compare switches
    TBaseMode tbase_mode = TBaseMode::Fixed;     // --tbase (default), --t-upscaled, --t-fine
    AlignMode align_mode = AlignMode::Resample;  // --resample (default), --make-uniform

    bool plotter = false;

    // NEW: force upscaled-only using hardcoded mean params (no file/folder loading)
    bool hardcoded_mean = true;

    // hardcoded mean params (ONLY used when hardcoded_mean == true)
    // >>> EDIT THESE <<<
    double lc_mean_h   = 0.39753;
    double lx_mean_h   = 1.25394;
    double ly_mean_h   = 0.125017;
    double dt_mean_h   = 7.31122e-05;
    double qx_const_h  = 1.0;   // constant qx used to build invcdf_mean in hardcoded mode

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--upscale-only") upscale_only = true;
        else if (a == "--mean-ts") use_timeseriesset_mean = true;
        else if (a == "--no-mean-ts") use_timeseriesset_mean = false;
        else if (a.rfind("--run-dir=", 0) == 0) resume_run_dir = a.substr(std::string("--run-dir=").size());

        else if (a == "--tbase")      tbase_mode = TBaseMode::Fixed;
        else if (a == "--t-upscaled") tbase_mode = TBaseMode::FromUpscaled;
        else if (a == "--t-fine")     tbase_mode = TBaseMode::FromFirstFine;

        else if (a == "--resample")     align_mode = AlignMode::Resample;
        else if (a == "--make-uniform") align_mode = AlignMode::MakeUniform;

        else if (a == "--plotter") plotter = true;
        else if (a == "--hardcoded-mean") hardcoded_mean = true;
    }

    // If hardcoded_mean is enabled, we must skip fine loop
    if (hardcoded_mean) upscale_only = true;

    std::string run_dir;
    if (rank == 0) {
        createDirectory(output_dir);

        if (upscale_only) {

            if (hardcoded_mean) {
                // hardcoded mode: we DON'T need a previous run folder.
                // If user provided --run-dir and it exists, use it; otherwise create a fresh run folder.
                if (!resume_run_dir.empty() && dirExists(resume_run_dir)) {
                    run_dir = resume_run_dir;
                    std::cout << "Hardcoded-mean upscaled-only: using run_dir: " << run_dir << "\n";
                } else {
                    run_dir = joinPath(output_dir, "run_" + makeTimestamp());
                    createDirectory(run_dir);
                    std::cout << "Hardcoded-mean upscaled-only: created run_dir: " << run_dir << "\n";
                }
            } else {
                // ORIGINAL upscale-only behavior (requires existing run-dir)
                if (resume_run_dir.empty()) {
                    std::cerr << "ERROR: --upscale-only requires --run-dir=/path/to/run_YYYYMMDD_HHMMSS\n";
                    MPI_Abort(PETSC_COMM_WORLD, 1);
                }
                run_dir = resume_run_dir;
                if (!dirExists(run_dir)) {
                    std::cerr << "ERROR: run_dir does not exist: " << run_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 2);
                }
                std::cout << "Upscale-only: using run_dir: " << run_dir << "\n";
            }

        } else {
            run_dir = joinPath(output_dir, "run_" + makeTimestamp());
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

    // -----------------------------
    // Domain / grid resolution
    // -----------------------------
    int nx = 300;
    int nu = 100;
    int ny = 100;
    double Lx = 3.0;
    double Ly = 1.0;
    double correlation_ls_x = 1;
    double correlation_ls_y = 0.1;
    double stdev = 2.0;
    double g_mean = 0;
    double Diffusion_coefficient = 0.01;
    double t_end_pdf = 20;

    // -----------------------------
    // Realizations
    // -----------------------------
    const int nReal_default = 20;
    const double du = 1 / double(nu);

    std::vector<double> lc_all, lx_all, ly_all, dt_all;
    TimeSeriesSet<double> inverse_qx_cdfs;
    TimeSeriesSet<double> qx_pdfs;
    TimeSeriesSet<double> lambda_x_correlations;
    TimeSeriesSet<double> lambda_y_correlations;
    TimeSeriesSet<double> lambda_a_correlations;
    std::vector<TimeSeriesSet<double>> Fine_Scale_BTCs;

    std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (!upscale_only && rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,dt_opt\n";
    }

    // =====================================================================
    // FINE-SCALE LOOP (skipped in --upscale-only or --hardcoded-mean)
    // =====================================================================
    int nReal = nReal_default;

    std::vector<double> xLocations{0.5, 1.5, 2.5};
    Fine_Scale_BTCs.resize(xLocations.size());

    if (!upscale_only) {
        for (int r = 1; r <= nReal; ++r) {
            const std::string rlab = makeRealLabel(r);
            std::string fine_dir = joinPath(run_dir, makeFineFolder(r));
            if (rank == 0) createDirectory(fine_dir);
            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

            const std::string pfx = rlab + "_";

            Grid2D g(nx, ny, Lx, Ly);

            PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

            unsigned long run_seed = 20260115UL;
            unsigned long seed = run_seed + 1000UL * (unsigned long)r + (unsigned long)rank;
            g.makeGaussianFieldSGS("K_normal_score", correlation_ls_x, correlation_ls_y, 10, seed);
            g.normalizeField("K_normal_score", Grid2D::ArrayKind::Cell);

            PetscTime(&t_asm0);
            PetscTime(&t_total0);

            g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
            g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.vti"));

            g.createExponentialField("K_normal_score", stdev, g_mean, "K");
            g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));

            PetscTime(&t_asm1);

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

            TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
            TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
            g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);
            g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx_normal_score.vti"));

            // velocity autocorrelation (X,Y) via perturbation
            double delta_min = 0.001, delta_max = 0.2;
            int num_deltas = 30;
            int num_samples_per_delta = 10000;

            TimeSeries<double> corr_x, corr_y;

            for (int i = 0; i < num_deltas; ++i) {
                double exponent = static_cast<double>(i) / (num_deltas - 1);
                double delta = delta_min * std::pow(delta_max / delta_min, exponent);
                try {
                    TimeSeries<double> samples = g.sampleGaussianPerturbation(
                        "qx_normal_score", Grid2D::ArrayKind::Fx,
                        num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                    corr_x.append(delta, samples.correlation_tc());
                } catch (...) {}
            }
            corr_x.writefile(joinPath(fine_dir, pfx + "velocity_correlation_x.txt"));
            lambda_x_correlations.append(corr_x, "Realization" + aquiutils::numbertostring(r + 1));
            double lambda_x_emp = corr_x.fitExponentialDecay();

            for (int i = 0; i < num_deltas; ++i) {
                double exponent = static_cast<double>(i) / (num_deltas - 1);
                double delta = delta_min * std::pow(delta_max / delta_min, exponent);
                try {
                    TimeSeries<double> samples = g.sampleGaussianPerturbation(
                        "qx_normal_score", Grid2D::ArrayKind::Fx,
                        num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                    corr_y.append(delta, samples.correlation_tc());
                } catch (...) {}
            }
            corr_y.writefile(joinPath(fine_dir, pfx + "velocity_correlation_y.txt"));
            lambda_y_correlations.append(corr_y, "Realization" + aquiutils::numbertostring(r + 1));
            double lambda_y_emp = corr_y.fitExponentialDecay();

            // inverse CDF + pdf
            TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100, 1e-6);
            TimeSeries<double> qx_pdf = g.extractFieldPDF("qx", Grid2D::ArrayKind::Fx, 50, 1e-6, true);
            qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
            qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));
            qx_pdf.writefile(joinPath(fine_dir, pfx + "qx_pdf.txt"));

            // dt
            double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;

            // transport
            g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

            g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.txt"));
            g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.txt"));
            g.writeNamedMatrix("K",  Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.txt"));

            g.SetVal("diffusion", Diffusion_coefficient);
            g.SetVal("porosity", 1);
            g.SetVal("c_left", 1.0);

            TimeSeriesSet<double> BTCs_FineScaled;
            g.setBTCLocations(xLocations);
            if (solve_fine_scale_transport) {
                g.SolveTransport(
                    t_end_pdf,
                    std::min(dt_optimal, 0.5 / 10.0),
                    "transport_", 500,
                    fine_dir,
                    "C",
                    &BTCs_FineScaled,
                    r
                );

                for (int i = 0; i < (int)xLocations.size() && i < (int)BTCs_FineScaled.size(); ++i) {
                    BTCs_FineScaled.setSeriesName(i, fmt_x(xLocations[i]));
                    Fine_Scale_BTCs[i].append(BTCs_FineScaled[i].derivative(), pfx + fmt_x(xLocations[i]));
                }

                const std::string btc_path       = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
                const std::string btc_deriv_path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

                BTCs_FineScaled.write(btc_path);
                BTCs_FineScaled.derivative().write(btc_deriv_path);
            }

            PetscTime(&t_total1);

            // pathway correlation -> lc
            PathwaySet pathways;
            pathways.Initialize(1000, PathwaySet::Weighting::FluxWeighted, &g);
            pathways.trackAllPathways(&g, 0.01);

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
            lambda_a_correlations.append(qx_correlation, "Realization" + aquiutils::numbertostring(r + 1));
            double advection_correlation_length_scale = qx_correlation.fitExponentialDecay();

            pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
            pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"));

            // meta
            if (rank == 0) {
                std::ofstream m(joinPath(fine_dir, pfx + "meta.txt"));
                m << "realization=" << r << "\n";
                m << "nx=" << nx << "\nny=" << ny << "\nnu=" << nu << "\n";
                m << "Lx=" << Lx << "\nLy=" << Ly << "\n";
                m << "D=" << Diffusion_coefficient << "\n";
                m << "correlation_ls_x=" << correlation_ls_x << "\n";
                m << "correlation_ls_y=" << correlation_ls_y << "\n";
                m << "stdev=" << stdev << "\n";
                m << "g_mean=" << g_mean << "\n";
                m << "lc=" << advection_correlation_length_scale << "\n";
                m << "lambda_x=" << lambda_x_emp << "\n";
                m << "lambda_y=" << lambda_y_emp << "\n";
                m << "dt_opt=" << dt_optimal << "\n";
                m << "seed=" << seed << "\n";
            }

            // accumulate means
            if (rank == 0) {
                lc_all.push_back(advection_correlation_length_scale);
                lx_all.push_back(lambda_x_emp);
                ly_all.push_back(lambda_y_emp);
                dt_all.push_back(dt_optimal);

                inverse_qx_cdfs.append(qx_inverse_cdf, "qx_inverse_cdf" + aquiutils::numbertostring(r + 1));
                qx_pdfs.append(qx_pdf, "qx_pdf" + aquiutils::numbertostring(r + 1));
                std::ofstream f(stats_csv, std::ios::app);
                f << r << "," << advection_correlation_length_scale << ","
                  << lambda_x_emp << "," << lambda_y_emp << "," << dt_optimal << "\n";
            }

            if (rank == 0) {
                std::cout << "[Fine " << rlab << "] Assembly time: " << (t_asm1 - t_asm0) << " s\n";
                std::cout << "[Fine " << rlab << "] Solve time:    " << (t_solve1 - t_solve0) << " s\n";
                std::cout << "[Fine " << rlab << "] Total time:   " << (t_total1 - t_total0) << " s\n";
                std::cout << "[Fine " << rlab << "] Outputs saved to: " << fine_dir << "\n";
            }

            MPI_Barrier(PETSC_COMM_WORLD);
        }
    }

    // Means from fine derivatives (kept)
    TimeSeriesSet<double> mean_BTCs;
    for (int i = 0; i < (int)xLocations.size(); i++)
        mean_BTCs.append(Fine_Scale_BTCs[i].mean_ts(), fmt_x(xLocations[i]));

    if (!upscale_only) {
        inverse_qx_cdfs.write(joinPath(run_dir, "qx_inverse_cdfs.txt"));
        qx_pdfs.write(joinPath(run_dir, "qx_pdfs.txt"));
        qx_pdfs.mean_ts().writefile(joinPath(run_dir, "qx_mean_pdf.txt"));

        lambda_a_correlations.write(joinPath(run_dir, "advective_correlations.txt"));
        lambda_x_correlations.write(joinPath(run_dir, "diffusion_x_correlations.txt"));
        lambda_y_correlations.write(joinPath(run_dir, "diffusion_y_correlations.txt"));
    }

    // =====================================================================
    // MEAN PARAMS + UPSCALED RUN
    // =====================================================================
    double lc_mean = 0, lx_mean = 0, ly_mean = 0, dt_mean = 0;
    TimeSeries<double> invcdf_mean;

    if (rank == 0) {

        if (hardcoded_mean) {
            // No file loading, no folder scanning, no reconstruction
            lc_mean = lc_mean_h;
            lx_mean = lx_mean_h;
            ly_mean = ly_mean_h;
            dt_mean = dt_mean_h;

            invcdf_mean.clear();
            invcdf_mean.append(0.0, qx_const_h);
            invcdf_mean.append(1.0, qx_const_h);

            // bookkeeping
            {
                std::ofstream f(joinPath(run_dir, "mean_params_used.txt"));
                f << "lc_mean=" << lc_mean << "\n";
                f << "lambda_x_mean=" << lx_mean << "\n";
                f << "lambda_y_mean=" << ly_mean << "\n";
                f << "dt_mean=" << dt_mean << "\n";
                f << "qx_const=" << qx_const_h << "\n";
            }

            std::cout << "Using HARD-CODED mean params (no loading).\n";
        }
        else if (!upscale_only) {
            nReal = (int)lc_all.size();
            lc_mean = mean_of(lc_all);
            lx_mean = mean_of(lx_all);
            ly_mean = mean_of(ly_all);
            dt_mean = mean_of(dt_all);
            invcdf_mean = inverse_qx_cdfs.mean_ts();
            invcdf_mean.writefile(joinPath(run_dir, "mean_qx_inverse_cdf.txt"));

            std::ofstream f(joinPath(run_dir, "mean_params.txt"));
            f << "nReal=" << nReal << "\n";
            f << "lc_mean=" << lc_mean << "\n";
            f << "lambda_x_mean=" << lx_mean << "\n";
            f << "lambda_y_mean=" << ly_mean << "\n";
            f << "dt_mean=" << dt_mean << "\n";
            f << "du=" << du << "\n";

        } else {
            const std::string mean_params_path = joinPath(run_dir, "mean_params.txt");
            const std::string mean_cdf_path    = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

            bool ok_params = fileExists(mean_params_path) &&
                             read_mean_params_txt(mean_params_path, lc_mean, lx_mean, ly_mean, dt_mean);

            TimeSeries<double> mean_qx_cdf;
            bool ok_cdf = fileExists(mean_cdf_path) &&
                          read_mean_inverse_cdf_csv(mean_cdf_path, mean_qx_cdf);

            if (ok_params && ok_cdf) {
                std::cout << "Loaded mean_params.txt and mean_qx_inverse_cdf.txt\n";
                invcdf_mean = mean_qx_cdf;
            } else {
                std::cout << "Mean files missing/incomplete; reconstructing from fine_* folders...\n";

                auto fine_folders = list_fine_folders(run_dir);
                if (fine_folders.empty()) {
                    std::cerr << "ERROR: No fine_* folders found under: " << run_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 3);
                }
                nReal = (int)fine_folders.size();

                std::vector<double> lc_v, lx_v, ly_v, dt_v;
                bool ok_fine_params = fileExists(stats_csv) &&
                                      read_fine_params_all_csv(stats_csv, lc_v, lx_v, ly_v, dt_v);

                if (ok_fine_params) {
                    lc_mean = mean_of(lc_v);
                    lx_mean = mean_of(lx_v);
                    ly_mean = mean_of(ly_v);
                    dt_mean = mean_of(dt_v);
                    std::cout << "Loaded params from fine_params_all.csv\n";
                } else {
                    lc_v.clear(); lx_v.clear(); ly_v.clear(); dt_v.clear();
                    for (auto& pr : fine_folders) {
                        int rr = pr.first;
                        std::string fine_dir = pr.second;
                        std::string pfx = makeRealLabel(rr) + "_";
                        std::string meta = joinPath(fine_dir, pfx + "meta.txt");

                        std::map<std::string, std::string> kv;
                        if (!parse_keyval_file(meta, kv)) continue;

                        auto getd = [&](const std::string& k, double& out)->bool{
                            auto it = kv.find(k);
                            if (it == kv.end()) return false;
                            out = std::atof(it->second.c_str());
                            return true;
                        };

                        double v1, v2, v3, v4;
                        if (getd("lc", v1) && getd("lambda_x", v2) && getd("lambda_y", v3) && getd("dt_opt", v4)) {
                            lc_v.push_back(v1); lx_v.push_back(v2); ly_v.push_back(v3); dt_v.push_back(v4);
                        }
                    }
                    if (lc_v.empty()) {
                        std::cerr << "ERROR: Could not read params from fine_params_all.csv or meta.txt files.\n";
                        MPI_Abort(PETSC_COMM_WORLD, 4);
                    }
                    lc_mean = mean_of(lc_v);
                    lx_mean = mean_of(lx_v);
                    ly_mean = mean_of(ly_v);
                    dt_mean = mean_of(dt_v);
                    std::cout << "Loaded params from fine meta.txt files\n";
                }

                // Write mean_params.txt
                std::ofstream f(joinPath(run_dir, "mean_params.txt"));
                f << "lc_mean=" << lc_mean << "\n";
                f << "lambda_x_mean=" << lx_mean << "\n";
                f << "lambda_y_mean=" << ly_mean << "\n";
                f << "dt_mean=" << dt_mean << "\n";
                f << "du=" << du << "\n";

                // IMPORTANT: do NOT silently fall back to constant qx in non-hardcoded mode
                std::cerr << "ERROR: mean_qx_inverse_cdf.txt missing and reconstruction of mean CDF is not implemented here.\n"
                          << "Either provide mean_qx_inverse_cdf.txt, or run with --hardcoded-mean.\n";
                MPI_Abort(PETSC_COMM_WORLD, 99);
            }
        }
    }

    // Broadcast means + invcdf_mean to all ranks
    MPI_Bcast(&lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // Upscaled folder + prefix
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
                return 1;
            }
            qx_u[id] = v_at_u;
        }
    }

    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.txt"));

    g_u.setMixingParams(lc_mean, lx_mean, ly_mean);

    double dt_pdf = dt_mean;
    int output_interval_pdf = 500;

    g_u.SetVal("diffusion", Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));
    g_u.writeNamedMatrix("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.txt"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
    g_u.setBTCLocations(xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    if (solve_upscale_transport) {
        g_u.SolveTransport(t_end_pdf, dt_pdf, "transport_", output_interval_pdf, up_dir, "Cu", &BTCs_Upscaled);

        for (int i = 0; i < (int)xLocations.size() && i < (int)BTCs_Upscaled.size(); ++i) {
            BTCs_Upscaled.setSeriesName(i, fmt_x(xLocations[i]));
            Fine_Scale_BTCs[i].append(BTCs_Upscaled[i].derivative(), "Upscaled" + fmt_x(xLocations[i]));
        }

        BTCs_Upscaled.write(up_btc_path);
        BTCs_Upscaled.derivative().write(up_btc_deriv_path);

        g_u.writeNamedVTI_Auto("C", joinPath(up_dir, "Cu.vti"));
        g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell, joinPath(up_dir, up_pfx + "Cu.txt"));
    }

    // Your per-location + mean CSV outputs (kept)
    for (int i = 0; i < (int)xLocations.size(); i++) {
        std::string out_cmp_d = joinPath(run_dir, fmt_x(xLocations[i]) + "BTC_Compare.csv");
        Fine_Scale_BTCs[i].write(out_cmp_d);
        std::cout << "Wrote derivative WIDE comparison CSV (t_base): " << out_cmp_d << "\n";
    }
    {
        std::string out_cmp_d = joinPath(run_dir, "BTC_mean.csv");
        mean_BTCs.write(out_cmp_d);
    }

    // =====================================================================
    // FINAL AGGREGATION CSVs for comparison / plotting
    // =====================================================================
    if (plotter) {
        if (rank == 0) {
            const bool ok = run_final_aggregation_and_plots(
                run_dir,
                up_btc_path,
                up_btc_deriv_path,
                tbase_mode,
                align_mode,
                use_timeseriesset_mean,
                /*t_end_cmp=*/10.0,
                /*dt_cmp=*/0.001
            );

            if (!ok) {
                std::cerr << "WARNING: final aggregation/plotting reported failure.\n";
            }
        }
    }
    // =====================================================================

    std::cout << "\nMixing PDF simulation complete!\n";
    std::cout << "All outputs saved to: " << run_dir << "\n";
    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
