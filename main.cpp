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

// NEW: moved helpers
#include "sim_helpers.h"

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
#else
    std::string output_dir = "./Results";
#endif

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // -----------------------------
    // CLI options
    // -----------------------------
    bool upscale_only = false;          // read fine outputs and run only upscaled
    std::string resume_run_dir = output_dir + "/run_20260115_132010";    // required for --upscale-only

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--upscale-only") upscale_only = true;
        else if (a.rfind("--run-dir=", 0) == 0) resume_run_dir = a.substr(std::string("--run-dir=").size());
    }

    std::string run_dir;
    if (rank == 0) {
        createDirectory(output_dir);

        if (upscale_only) {
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
    double stdev = 1.0;
    double g_mean = 0;
    double Diffusion_coefficient = 0;

    // -----------------------------
    // Realizations
    // -----------------------------
    const int nReal_default = 1; // used only in full mode; in upscale-only we detect from folders
    const double du = 0.01;
    const int nU = (int)std::round(1.0 / du) + 1;

    std::vector<double> lc_all, lx_all, ly_all, dt_all;
    std::vector<double> invcdf_sum(nU, 0.0);

    std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (!upscale_only && rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,dt_opt\n";
    }

    // =====================================================================
    // FINE-SCALE LOOP (skipped in --upscale-only)
    // =====================================================================
    int nReal = nReal_default;

    if (!upscale_only) {
        for (int r = 1; r <= nReal; ++r) {
            const std::string rlab = makeRealLabel(r);

            std::string fine_dir = joinPath(run_dir, makeFineFolder(r));
            if (rank == 0) createDirectory(fine_dir);
            MPI_Barrier(PETSC_COMM_WORLD);
            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

            const std::string pfx = rlab + "_";

            Grid2D g(nx, ny, Lx, Ly);

            PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

            unsigned long run_seed = 20260115UL;
            unsigned long seed = run_seed + 1000UL*(unsigned long)r + (unsigned long)rank;
            g.makeGaussianFieldSGS("K_normal_score", correlation_ls_x, correlation_ls_y, 10, seed);

            PetscTime(&t_asm0);
            PetscTime(&t_total0);

            g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
            g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "NormalScore.vti"));

            g.createExponentialField("K_normal_score", stdev, g_mean, "K");
            PetscTime(&t_asm1);

            PetscTime(&t_solve0);
            g.DarcySolve(Lx, 0, "K", "K");
            std::cout << "Darcy solved ... " << std::endl;
            g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));
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

            std::cout << "Sampling points for derivative ..." << std::endl;
            TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score", Grid2D::ArrayKind::Fx,
                                                                   Grid2D::DerivDir::X, 10000, 0.05);
            curvture.writefile(joinPath(fine_dir, pfx + "2nd_deriv.txt"));

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
            double lambda_y_emp = corr_y.fitExponentialDecay();

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

            std::vector<double> xLocations{0.5, 1.5, 2.5};
            g.setBTCLocations(xLocations);

            TimeSeriesSet<double> BTCs_FineScaled;
            g.SolveTransport(10,
                             std::min(dt_optimal, 0.5 / 10.0),
                             "transport_", 50,
                             fine_dir,
                             "C",
                             &BTCs_FineScaled,
                             r);

            g.writeNamedVTI_Auto("C", joinPath(fine_dir, pfx + "C.vti"));

            const std::string btc_path       = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
            const std::string btc_deriv_path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

            BTCs_FineScaled.write(btc_path);
            BTCs_FineScaled.derivative().write(btc_deriv_path);

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
            double advection_correlation_length_scale = qx_correlation.fitExponentialDecay();

            pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
            pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"));

            // inverse CDF
            TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100);
            qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
            qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));

            // meta
            if (rank == 0) {
                std::ofstream m(joinPath(fine_dir, pfx + "meta.txt"));
                m << "realization=" << r << "\n";
                m << "nx=" << nx << "\nny=" << ny << "\nnu=" << nu << "\n";
                m << "Lx=" << Lx << "\nLy=" << Ly << "\n";
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

                for (int k = 0; k < nU; ++k) {
                    double u = k * du;
                    invcdf_sum[k] += qx_inverse_cdf.interpol(u);
                }

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

    // =====================================================================
    // MEAN PARAMS + UPSCALED RUN
    // =====================================================================
    double lc_mean = 0, lx_mean = 0, ly_mean = 0, dt_mean = 0;
    std::vector<double> invcdf_mean;

    if (rank == 0) {
        if (!upscale_only) {
            nReal = (int)lc_all.size();
            lc_mean = mean_of(lc_all);
            lx_mean = mean_of(lx_all);
            ly_mean = mean_of(ly_all);
            dt_mean = mean_of(dt_all);

            invcdf_mean.resize(nU, 0.0);
            for (int k = 0; k < nU; ++k) invcdf_mean[k] = invcdf_sum[k] / (double)nReal;

            {
                std::ofstream f(joinPath(run_dir, "mean_params.txt"));
                f << "nReal=" << nReal << "\n";
                f << "lc_mean=" << lc_mean << "\n";
                f << "lambda_x_mean=" << lx_mean << "\n";
                f << "lambda_y_mean=" << ly_mean << "\n";
                f << "dt_mean=" << dt_mean << "\n";
                f << "du=" << du << "\n";
            }
            {
                std::ofstream f(joinPath(run_dir, "mean_qx_inverse_cdf.txt"));
                f << "u,v\n";
                for (int k = 0; k < nU; ++k) {
                    double u = k * du;
                    f << u << "," << invcdf_mean[k] << "\n";
                }
            }
        } else {
            // upscale-only: try load mean files, else reconstruct from fine folders
            const std::string mean_params_path = joinPath(run_dir, "mean_params.txt");
            const std::string mean_cdf_path    = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

            bool ok_params = fileExists(mean_params_path) &&
                             read_mean_params_txt(mean_params_path, lc_mean, lx_mean, ly_mean, dt_mean);

            std::vector<double> uvec, vvec;
            bool ok_cdf = fileExists(mean_cdf_path) &&
                          read_mean_inverse_cdf_csv(mean_cdf_path, uvec, vvec);

            if (ok_params && ok_cdf) {
                invcdf_mean.assign(nU, 0.0);
                for (int k = 0; k < nU; ++k) {
                    double u = k * du;
                    invcdf_mean[k] = interp1_linear(uvec, vvec, u);
                }
                std::cout << "Loaded mean_params.txt and mean_qx_inverse_cdf.txt\n";
            } else {
                std::cout << "Mean files missing/incomplete; reconstructing from fine_* folders...\n";

                auto fine_folders = list_fine_folders(run_dir);
                if (fine_folders.empty()) {
                    std::cerr << "ERROR: No fine_* folders found under: " << run_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 3);
                }
                nReal = (int)fine_folders.size();

                // params from fine_params_all.csv if exists
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
                    // fallback: read meta.txt per folder
                    lc_v.clear(); lx_v.clear(); ly_v.clear(); dt_v.clear();
                    for (auto& pr : fine_folders) {
                        int rr = pr.first;
                        std::string fine_dir = pr.second;
                        std::string pfx = makeRealLabel(rr) + "_";
                        std::string meta = joinPath(fine_dir, pfx + "meta.txt");

                        std::map<std::string,std::string> kv;
                        if (!parse_keyval_file(meta, kv)) continue;

                        auto getd = [&](const std::string& k, double& out)->bool{
                            auto it = kv.find(k);
                            if (it == kv.end()) return false;
                            out = std::atof(it->second.c_str());
                            return true;
                        };

                        double v1,v2,v3,v4;
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

                // inverse CDF average from per-realization qx_inverse_cdf.txt
                std::fill(invcdf_sum.begin(), invcdf_sum.end(), 0.0);
                int used = 0;
                for (auto& pr : fine_folders) {
                    int rr = pr.first;
                    std::string fine_dir = pr.second;
                    if (accumulate_inverse_cdf_on_grid(fine_dir, rr, du, nU, invcdf_sum)) used++;
                }
                if (used == 0) {
                    std::cerr << "ERROR: No qx_inverse_cdf.txt files could be read.\n";
                    MPI_Abort(PETSC_COMM_WORLD, 5);
                }

                invcdf_mean.resize(nU, 0.0);
                for (int k = 0; k < nU; ++k) invcdf_mean[k] = invcdf_sum[k] / (double)used;

                // save reconstructed means for next resume
                {
                    std::ofstream f(joinPath(run_dir, "mean_params.txt"));
                    f << "nReal=" << used << "\n";
                    f << "lc_mean=" << lc_mean << "\n";
                    f << "lambda_x_mean=" << lx_mean << "\n";
                    f << "lambda_y_mean=" << ly_mean << "\n";
                    f << "dt_mean=" << dt_mean << "\n";
                    f << "du=" << du << "\n";
                }
                {
                    std::ofstream f(joinPath(run_dir, "mean_qx_inverse_cdf.txt"));
                    f << "u,v\n";
                    for (int k = 0; k < nU; ++k) {
                        double u = k * du;
                        f << u << "," << invcdf_mean[k] << "\n";
                    }
                }

                std::cout << "Reconstruction done from " << used << " fine realizations.\n";
            }
        }
    }

    // broadcast scalar means
    MPI_Bcast(&lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // broadcast mean inverse cdf vector
    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // build mean inverse CDF TimeSeries
    TimeSeries<double> qx_inverse_cdf_mean;
    for (int k = 0; k < nU; ++k) {
        double u = k * du;
        qx_inverse_cdf_mean.append(u, invcdf_mean[k]);
    }

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

    // Mixing PDF grid
    Grid2D g_u(nx, ny, Lx, Ly);

    int qx_size = nu * (nx + 1);
    int qy_size = (nu + 1) * nx;

    auto& qx_u = g_u.flux("qx");
    auto& qy_u = g_u.flux("qy");

    qx_u.resize(qx_size, 0.0);
    qy_u.resize(qy_size, 0.0);

    // Assign v(u) from mean inverse CDF
    for (int j = 0; j < nu; ++j) {
        double u = static_cast<double>(j) / (nu - 1);
        double v_at_u = qx_inverse_cdf_mean.interpol(u);
        for (int i = 0; i < nx + 1; ++i) {
            int id = j * (nx + 1) + i;
            if (id >= qx_size) {
                std::cerr << "ERROR: qx index out of bounds: " << id << " >= " << qx_size << "\n";
                return 1;
            }
            qx_u[id] = v_at_u;
        }
    }

    // Save initial fields with prefix
    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.txt"));

    // Mixing parameters (mean)
    g_u.setMixingParams(lc_mean, lx_mean, ly_mean);

    double t_end_pdf = 10;
    double dt_pdf = dt_mean;
    int output_interval_pdf = 50;

    g_u.SetVal("diffusion", Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

    std::vector<double> xLocations{0.5, 1.5, 2.5};
    g_u.setBTCLocations(xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    g_u.SolveTransport(t_end_pdf, dt_pdf, "transport_", output_interval_pdf, up_dir, "Cu", &BTCs_Upscaled);

    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    BTCs_Upscaled.write(up_btc_path);
    BTCs_Upscaled.derivative().write(up_btc_deriv_path);

    g_u.writeNamedVTI_Auto("C", joinPath(up_dir, "Cu.vti"));
    g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell, joinPath(up_dir, up_pfx + "Cu.txt"));

    // =====================================================================
    // FINAL AGGREGATION CSVs for comparison / plotting
    //   Requested order:
    //     r0001..r0020, FineMean, Upscaled
    // =====================================================================
    if (rank == 0) {

        const double t_end_cmp = 10.0;
        const double dt_cmp    = 0.001;

        std::vector<double> t_base;
        t_base.reserve((size_t)std::ceil(t_end_cmp / dt_cmp) + 1);
        for (double tt = 0.0; tt <= t_end_cmp + 1e-12; tt += dt_cmp) t_base.push_back(tt);

        std::vector<std::string> out_names;
        std::vector<std::vector<double>> out_cols;

        auto ingest_one = [&](const std::string& path, const std::string& series_prefix) -> bool {
            std::vector<double> t;
            std::vector<std::string> names;
            std::vector<std::vector<double>> cols;
            if (!read_time_series_table_csv(path, t, names, cols)) {
                std::cerr << "WARNING: could not read CSV for aggregation: " << path << "\n";
                return false;
            }

            std::vector<std::vector<double>> cols_rs;
            resample_table_linear(t, cols, t_base, cols_rs);

            for (size_t j = 0; j < names.size(); ++j) {
                out_names.push_back(series_prefix + "_" + names[j]);
                out_cols.push_back(std::move(cols_rs[j]));
            }
            return true;
        };

        // scan fine folders
        auto fine_folders = list_fine_folders(run_dir);

        // ------------------------------------------------------------
        // BTC comparison: r0001.., FineMean, Upscaled
        // ------------------------------------------------------------

        // 1) ingest all fine BTCs in order
        for (auto& pr : fine_folders) {
            int rr = pr.first;
            std::string fine_dir = pr.second;
            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

            std::string pfx = makeRealLabel(rr) + "_";
            std::string btc = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
            ingest_one(btc, "Fine_" + makeRealLabel(rr));
        }

        // 2) compute FineMean from disk (after all r####)
        {
            bool initialized = false;
            std::vector<std::string> mean_names;
            std::vector<std::vector<double>> sum_cols;
            std::vector<std::vector<int>>    cnt_cols;
            int used = 0;

            for (auto& pr : fine_folders) {
                int rr = pr.first;
                std::string fine_dir = pr.second;
                if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

                std::string pfx = makeRealLabel(rr) + "_";
                std::string btc = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");

                std::vector<double> t;
                std::vector<std::string> names;
                std::vector<std::vector<double>> cols;
                if (!read_time_series_table_csv(btc, t, names, cols)) continue;

                std::vector<std::vector<double>> cols_rs;
                resample_table_linear(t, cols, t_base, cols_rs);

                if (!initialized) {
                    mean_names = names;
                    sum_cols.assign(mean_names.size(), std::vector<double>{});
                    cnt_cols.assign(mean_names.size(), std::vector<int>{});
                    initialized = true;
                }
                if (names.size() != mean_names.size()) continue;

                for (size_t j = 0; j < mean_names.size(); ++j) {
                    accumulate_sum_count(cols_rs[j], sum_cols[j], cnt_cols[j]);
                }
                used++;
            }

            if (initialized && used > 0) {
                for (size_t j = 0; j < mean_names.size(); ++j) {
                    out_names.push_back(std::string("FineMean_") + mean_names[j]);
                    out_cols.push_back(finalize_mean_vec(sum_cols[j], cnt_cols[j]));
                }
                std::cout << "Added FineMean (BTC) from " << used << " realizations.\n";
            } else {
                std::cerr << "WARNING: FineMean (BTC) not added (no readable fine BTC files).\n";
            }
        }

        // 3) ingest upscaled last
        ingest_one(up_btc_path, "Upscaled_mean");

        const std::string out_cmp = joinPath(run_dir, "BTC_Compare_Fine_vs_Upscaled.csv");
        if (!write_comparison_csv(out_cmp, t_base, out_names, out_cols)) {
            std::cerr << "WARNING: failed to write " << out_cmp << "\n";
        } else {
            std::cout << "Wrote comparison BTC CSV: " << out_cmp << "\n";
        }

        /*
        // python plot
        const std::string py1 = joinPath(run_dir, "plot_BTC_compare.py");
        if (write_btc_compare_plot_py(py1, out_cmp, "BTC_compare", "Concentration")) {
            int rc = run_python_script(py1);
            if (rc != 0) std::cerr << "WARNING: plotting script failed (rc=" << rc << "): " << py1 << "\n";
        } else {
            std::cerr << "WARNING: could not write plot script: " << py1 << "\n";
        }

        // gnuplot
        const std::string gp1 = joinPath(run_dir, "plot_BTC_compare.gp");
        if (write_btc_compare_plot_gnuplot(gp1, out_cmp, "BTC_compare", "Concentration")) {
            int rc = run_gnuplot_script(gp1);
            if (rc != 0)
                std::cerr << "WARNING: gnuplot failed for BTC comparison\n";
        }

        // simple gnuplot
        const std::string gp1 = joinPath(run_dir, "plot_BTC_compare.gp");
        if (write_btc_compare_plot_gnuplot_simple(gp1, out_cmp, "BTC_compare", "Concentration")) {
            int rc = run_gnuplot_script(gp1);
            if (rc != 0) std::cerr << "WARNING: gnuplot failed (rc=" << rc << "): " << gp1 << "\n";
        }

        */
        // gnuplot by name
        const std::string gp1 = joinPath(run_dir, "plot_BTC_compare.gp");
        if (write_btc_compare_plot_gnuplot_by_basename(gp1, out_cmp, "BTC_compare", "Concentration")) {
            int rc = run_gnuplot_script(gp1);
            if (rc != 0) std::cerr << "WARNING: gnuplot failed (rc=" << rc << "): " << gp1 << "\n";
        }

        // ------------------------------------------------------------
        // Derivative comparison: r0001.., FineDerivMean, UpscaledDeriv
        // ------------------------------------------------------------
        out_names.clear();
        out_cols.clear();

        // 1) ingest all fine derivative BTCs
        for (auto& pr : fine_folders) {
            int rr = pr.first;
            std::string fine_dir = pr.second;
            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

            std::string pfx = makeRealLabel(rr) + "_";
            std::string btc = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");
            ingest_one(btc, "FineDeriv_" + makeRealLabel(rr));
        }

        // 2) FineDerivMean from disk
        {
            bool initialized = false;
            std::vector<std::string> mean_names;
            std::vector<std::vector<double>> sum_cols;
            std::vector<std::vector<int>>    cnt_cols;
            int used = 0;

            for (auto& pr : fine_folders) {
                int rr = pr.first;
                std::string fine_dir = pr.second;
                if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

                std::string pfx = makeRealLabel(rr) + "_";
                std::string btc = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

                std::vector<double> t;
                std::vector<std::string> names;
                std::vector<std::vector<double>> cols;
                if (!read_time_series_table_csv(btc, t, names, cols)) continue;

                std::vector<std::vector<double>> cols_rs;
                resample_table_linear(t, cols, t_base, cols_rs);

                if (!initialized) {
                    mean_names = names;
                    sum_cols.assign(mean_names.size(), std::vector<double>{});
                    cnt_cols.assign(mean_names.size(), std::vector<int>{});
                    initialized = true;
                }
                if (names.size() != mean_names.size()) continue;

                for (size_t j = 0; j < mean_names.size(); ++j) {
                    accumulate_sum_count(cols_rs[j], sum_cols[j], cnt_cols[j]);
                }
                used++;
            }

            if (initialized && used > 0) {
                for (size_t j = 0; j < mean_names.size(); ++j) {
                    out_names.push_back(std::string("FineDerivMean_") + mean_names[j]);
                    out_cols.push_back(finalize_mean_vec(sum_cols[j], cnt_cols[j]));
                }
                std::cout << "Added FineDerivMean from " << used << " realizations.\n";
            } else {
                std::cerr << "WARNING: FineDerivMean not added (no readable fine derivative files).\n";
            }
        }

        // 3) upscaled derivative last
        ingest_one(up_btc_deriv_path, "UpscaledDeriv_mean");

        const std::string out_cmp_d = joinPath(run_dir, "BTC_Compare_FineDerivative_vs_UpscaledDerivative.csv");
        if (!write_comparison_csv(out_cmp_d, t_base, out_names, out_cols)) {
            std::cerr << "WARNING: failed to write " << out_cmp_d << "\n";
        } else {
            std::cout << "Wrote comparison BTC-derivative CSV: " << out_cmp_d << "\n";
        }

        std::cout << "\nMixing PDF simulation complete!\n";
        std::cout << "All outputs saved to: " << run_dir << "\n";

        /*
        // python plot
        const std::string py2 = joinPath(run_dir, "plot_BTC_derivative_compare.py");
        if (write_btc_compare_plot_py(py2, out_cmp_d, "BTC_deriv_compare", "dC/dt")) {
            int rc = run_python_script(py2);
            if (rc != 0) std::cerr << "WARNING: plotting script failed (rc=" << rc << "): " << py2 << "\n";
        } else {
            std::cerr << "WARNING: could not write plot script: " << py2 << "\n";
        }

        // gnuplot
        const std::string gp2 = joinPath(run_dir, "plot_BTC_derivative_compare.gp");
        if (write_btc_compare_plot_gnuplot(gp2, out_cmp_d, "BTC_deriv_compare", "dC/dt")) {
            int rc = run_gnuplot_script(gp2);
            if (rc != 0)
                std::cerr << "WARNING: gnuplot failed for BTC derivative comparison\n";
        }

        // simple gnuplot
        const std::string gp2 = joinPath(run_dir, "plot_BTC_derivative_compare.gp");
        if (write_btc_compare_plot_gnuplot_simple(gp2, out_cmp_d, "BTC_deriv_compare", "dC/dt")) {
            int rc = run_gnuplot_script(gp2);
            if (rc != 0) std::cerr << "WARNING: gnuplot failed (rc=" << rc << "): " << gp2 << "\n";
        }

        */
        // gnuplot by name
        const std::string gp2 = joinPath(run_dir, "plot_BTC_derivative_compare.gp");
        if (write_btc_compare_plot_gnuplot_by_basename(gp2, out_cmp_d, "BTC_deriv_compare", "dC/dt")) {
            int rc = run_gnuplot_script(gp2);
            if (rc != 0) std::cerr << "WARNING: gnuplot failed (rc=" << rc << "): " << gp2 << "\n";
        }

    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
