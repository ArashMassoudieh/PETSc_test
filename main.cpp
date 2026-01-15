// main.cpp
#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"
#include <cmath>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <mpi.h>
#include <fstream>
#include "grid.h"
#include "TimeSeries.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx) { return j*nx + i; }

// NOTE: helper to detect boundary cells. Not used in this file, but can be used for BC assignment.
static inline bool on_bc(PetscInt i, PetscInt j, PetscInt nx, PetscInt ny, double x, double y, double Lx, double Ly) {
    const double eps = 1e-14;
    return (i==0) || (j==0) || (std::abs(x - Lx) < eps) || (std::abs(y - Ly) < eps);
}

// Helper function to create directory if it doesn't exist
bool createDirectory(const std::string& path) {
    struct stat info;

    // Check if directory already exists
    if (stat(path.c_str(), &info) == 0) {
        if (info.st_mode & S_IFDIR) {
            return true; // Directory exists
        } else {
            std::cerr << "Error: Path exists but is not a directory: " << path << std::endl;
            return false;
        }
    }

    // Try to create directory
#ifdef _WIN32
    if (mkdir(path.c_str()) != 0) {
#else
    if (mkdir(path.c_str(), 0755) != 0) {
#endif
        if (errno != EEXIST) {
            std::cerr << "Error creating directory " << path << ": " << strerror(errno) << std::endl;
            return false;
        }
    }

    std::cout << "Created output directory: " << path << std::endl;
    return true;
}

// Helper function to join path with filename
std::string joinPath(const std::string& dir, const std::string& filename) {
    if (dir.empty()) return filename;
    if (dir.back() == '/' || dir.back() == '\\') {
        return dir + filename;
    }
    return dir + "/" + filename;
}

// Timestamp: YYYYMMDD_HHMMSS
static std::string makeTimestamp() {
    std::time_t now = std::time(nullptr);
    std::tm tm_now;
    localtime_r(&now, &tm_now);
    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

static std::string makeRealName(int r) {
    std::ostringstream oss;
    oss << "fine_r" << std::setw(3) << std::setfill('0') << r;
    return oss.str();
}

static double mean_of(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / (double)v.size();
}

int main(int argc, char** argv) {
    // PETSc + MPI init RAII wrapper (assumed)
    PETScInit petsc(argc, argv);

    // Set output directory
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

    // --------------------------------------------------------------------
    // Create a run subfolder inside output_dir: output_dir/run_YYYYMMDD_HHMMSS/
    // --------------------------------------------------------------------
    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::string run_dir;
    if (rank == 0) {
        createDirectory(output_dir);
        run_dir = joinPath(output_dir, "run_" + makeTimestamp());
        createDirectory(run_dir);
        std::cout << "Run directory: " << run_dir << std::endl;
    }

    // Broadcast run_dir to all ranks so everyone writes to the same folder
    int len = 0;
    if (rank == 0) len = (int)run_dir.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    run_dir.resize(len);
    MPI_Bcast(run_dir.data(), len, MPI_CHAR, 0, PETSC_COMM_WORLD);

    // Ensure it ends with '/'
    if (!run_dir.empty() && run_dir.back() != '/' && run_dir.back() != '\\') run_dir += "/";

    MPI_Barrier(PETSC_COMM_WORLD);

    // --------------------------------------------------------------------
    // Domain / grid resolution
    // --------------------------------------------------------------------
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
    // Monte Carlo / realizations
    // -----------------------------
    const int nReal = 20;     // change to 10-20
    const double du = 0.01;   // u-grid for averaging inverse CDF
    const int nU = (int)std::round(1.0 / du) + 1;

    std::vector<double> lc_all;
    std::vector<double> lx_all;
    std::vector<double> ly_all;
    std::vector<double> dt_all;

    std::vector<double> invcdf_sum(nU, 0.0);

    std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,dt_opt\n";
    }

    // =====================================================================
    // FINE-SCALE LOOP (repeat nReal times)
    // =====================================================================
    for (int r = 0; r < nReal; ++r) {

        // realization folder
        std::string fine_dir = joinPath(run_dir, makeRealName(r));
        if (rank == 0) createDirectory(fine_dir);
        MPI_Barrier(PETSC_COMM_WORLD);
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        // Primary spatial grid (x,y)
        Grid2D g(nx, ny, Lx, Ly);

        // PETSc timers
        PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

        // Generate heterogeneous field
        g.makeGaussianFieldSGS("K_normal_score", correlation_ls_x, correlation_ls_y, 10);

        PetscTime(&t_asm0);
        PetscTime(&t_total0);

        // Write raw normal-score permeability field
        g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, "K_normal_score.txt"));
        g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, "NormalScore.vti"));

        // Transform normal score to exponential field "K"
        g.createExponentialField("K_normal_score", stdev, g_mean, "K");
        PetscTime(&t_asm1);

        // Darcy solve
        PetscTime(&t_solve0);
        g.DarcySolve(Lx, 0, "K", "K");
        std::cout << "Darcy solved ... " << std::endl;
        g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, "K.vti"));
        PetscTime(&t_solve1);

        // Diagnostics
        g.computeMassBalanceError("MassBalanceError");
        g.writeNamedMatrix("MassBalanceError", Grid2D::ArrayKind::Cell, joinPath(fine_dir, "error.txt"));

        // Save head + flux fields
        g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, "Head.vti"));
        g.writeNamedMatrix("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, "Head.txt"));
        g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, "qx.vti"));
        g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, "qy.vti"));

        // qx normal score
        TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
        TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
        g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);
        g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(fine_dir, "qx_normal_score.vti"));

        // curvature
        std::cout << "Sampling points for derivative ..." << std::endl;
        TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score", Grid2D::ArrayKind::Fx,
                                                               Grid2D::DerivDir::X, 10000, 0.05);
        curvture.writefile(joinPath(fine_dir, "2nd_deriv.txt"));

        // =============================================================================
        // Velocity autocorrelation analysis via Gaussian perturbation sampling
        // =============================================================================
        double delta_min = 0.001;
        double delta_max = 0.2;
        int num_deltas = 30;
        int num_samples_per_delta = 10000;

        TimeSeries<double> corr_x;
        TimeSeries<double> corr_y;

        std::cout << "Computing velocity autocorrelation (x-direction)..." << std::endl;
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
        corr_x.writefile(joinPath(fine_dir, "velocity_correlation_x.txt"));
        double lambda_x_emp = corr_x.fitExponentialDecay();
        std::cout << "X-direction correlation length scale: " << lambda_x_emp << std::endl;

        std::cout << "Computing velocity autocorrelation (y-direction)..." << std::endl;
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
        corr_y.writefile(joinPath(fine_dir, "velocity_correlation_y.txt"));
        double lambda_y_emp = corr_y.fitExponentialDecay();
        std::cout << "Y-direction correlation length scale: " << lambda_y_emp << std::endl;

        // CFL-like timestep estimate based on max qx
        double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;
        std::cout << "Optimal Time-Step: " << dt_optimal << std::endl;

        // Fine-scale transport
        g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

        g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, "qx.txt"));
        g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, "qy.txt"));
        g.writeNamedMatrix("K",  Grid2D::ArrayKind::Cell, joinPath(fine_dir, "K.txt"));

        g.SetVal("diffusion", Diffusion_coefficient);
        g.SetVal("porosity", 1);
        g.SetVal("c_left", 1.0);

        std::vector<double> xLocations{0.5, 1.5, 2.5};
        g.setBTCLocations(xLocations);

        TimeSeriesSet<double> BTCs_FineScaled;
        g.SolveTransport(10, std::min(dt_optimal, 0.5 / 10.0), "transport_", 50, fine_dir, "C", &BTCs_FineScaled);

        g.writeNamedVTI_Auto("C", joinPath(fine_dir, "C.vti"));

        BTCs_FineScaled.write(joinPath(fine_dir, "BTC_FineScaled.csv"));
        BTCs_FineScaled.derivative().write(joinPath(fine_dir, "BTC_FineScaled_derivative.csv"));

        PetscTime(&t_total1);

        // Many particle tracking
        PathwaySet pathways;
        pathways.Initialize(1000, PathwaySet::Weighting::FluxWeighted, &g);
        pathways.trackAllPathways(&g, 0.01);

        // Compute velocity correlation vs distance
        double Delta_x_min = 0.001;
        double Delta_x_max = 0.5;
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

        qx_correlation.writefile(joinPath(fine_dir, "qx_correlation_vs_distance.txt"));
        double advection_correlation_length_scale = qx_correlation.fitExponentialDecay();

        pathways.writeToFile(joinPath(fine_dir, "pathway_summary.txt"));
        pathways.writeCombinedVTK(joinPath(fine_dir, "all_pathways.vtk"));

        // Extract inverse CDF (quantile function)
        TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100);
        qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
        qx_inverse_cdf.writefile(joinPath(fine_dir, "qx_inverse_cdf.txt"));

        // Save quick meta
        if (rank == 0) {
            std::ofstream m(joinPath(fine_dir, "meta.txt"));
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
        }

        // Accumulate
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

        // Report timing per realization (rank 0)
        if (rank == 0) {
            std::cout << "[Fine " << r << "] Assembly time: " << (t_asm1 - t_asm0) << " s\n";
            std::cout << "[Fine " << r << "] Solve time:    " << (t_solve1 - t_solve0) << " s\n";
            std::cout << "[Fine " << r << "] Total time:   " << (t_total1 - t_total0) << " s\n";
            std::cout << "[Fine " << r << "] Outputs saved to: " << fine_dir << std::endl;
        }

        MPI_Barrier(PETSC_COMM_WORLD);
    } // end fine loop

    // =====================================================================
    // MEAN PARAMS + UPSCALED RUN (single)
    // =====================================================================
    double lc_mean = 0, lx_mean = 0, ly_mean = 0, dt_mean = 0;
    std::vector<double> invcdf_mean;

    if (rank == 0) {
        lc_mean = mean_of(lc_all);
        lx_mean = mean_of(lx_all);
        ly_mean = mean_of(ly_all);
        dt_mean = mean_of(dt_all);

        invcdf_mean.resize(nU, 0.0);
        for (int k = 0; k < nU; ++k) invcdf_mean[k] = invcdf_sum[k] / (double)nReal;

        // Write mean params
        {
            std::ofstream f(joinPath(run_dir, "mean_params.txt"));
            f << "nReal=" << nReal << "\n";
            f << "lc_mean=" << lc_mean << "\n";
            f << "lambda_x_mean=" << lx_mean << "\n";
            f << "lambda_y_mean=" << ly_mean << "\n";
            f << "dt_mean=" << dt_mean << "\n";
            f << "du=" << du << "\n";
        }

        // Write mean inverse CDF
        {
            std::ofstream f(joinPath(run_dir, "mean_qx_inverse_cdf.txt"));
            f << "u,v\n";
            for (int k = 0; k < nU; ++k) {
                double u = k * du;
                f << u << "," << invcdf_mean[k] << "\n";
            }
        }
    }

    // Broadcast scalar means
    MPI_Bcast(&lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // Broadcast mean inverse CDF vector
    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // Build TimeSeries mean inverse CDF
    TimeSeries<double> qx_inverse_cdf_mean;
    for (int k = 0; k < nU; ++k) {
        double u = k * du;
        qx_inverse_cdf_mean.append(u, invcdf_mean[k]);
    }

    // Upscaled output folder
    std::string up_dir = joinPath(run_dir, "upscaled_mean");
    if (rank == 0) createDirectory(up_dir);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (!up_dir.empty() && up_dir.back() != '/' && up_dir.back() != '\\') up_dir += "/";

    std::cout << "\n=== UPSCALED RUN USING MEAN PARAMETERS ===\n";
    if (rank == 0) {
        std::cout << "lc_mean=" << lc_mean << "\n";
        std::cout << "lambda_x_mean=" << lx_mean << "\n";
        std::cout << "lambda_y_mean=" << ly_mean << "\n";
        std::cout << "dt_mean=" << dt_mean << "\n";
        std::cout << "upscaled output: " << up_dir << "\n";
    }

    // ====================================================================
    // STEP 2: Create grid for mixing PDF in (x,u) space
    // ====================================================================
    Grid2D g_u(nx, ny, Lx, Ly);

    // ====================================================================
    // STEP 3: Initialize flux arrays
    // ====================================================================
    int qx_size = nu * (nx + 1);
    int qy_size = (nu + 1) * nx;

    auto& qx_u = g_u.flux("qx");
    auto& qy_u = g_u.flux("qy");

    qx_u.resize(qx_size, 0.0);
    qy_u.resize(qy_size, 0.0);

    // ====================================================================
    // STEP 4: Assign velocity field v(u) from mean inverse CDF
    // ====================================================================
    for (int j = 0; j < nu; ++j) {
        double u = static_cast<double>(j) / (nu - 1);
        double v_at_u = qx_inverse_cdf_mean.interpol(u);

        for (int i = 0; i < nx + 1; ++i) {
            int id = j * (nx + 1) + i;
            if (id >= qx_size) {
                std::cerr << "ERROR: qx index out of bounds: " << id << " >= " << qx_size << std::endl;
                return 1;
            }
            qx_u[id] = v_at_u;
        }
    }

    // Save initial fields
    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, "qy_u_initial.txt"));

    // ====================================================================
    // STEP 6: Set mixing parameters (MEANS)
    // ====================================================================
    g_u.setMixingParams(lc_mean, lx_mean, ly_mean);

    // Molecular diffusion + porosity + inlet PDF/BC
    g_u.SetVal("diffusion", Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    // ====================================================================
    // STEP 7: Solve mixing PDF equation
    // ====================================================================
    double t_end_pdf = 10;
    double dt_pdf = dt_mean;      // use mean dt
    int output_interval_pdf = 50;

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, "D_y.vti"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

    std::vector<double> xLocations{0.5, 1.5, 2.5};
    g_u.setBTCLocations(xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    g_u.SolveTransport(t_end_pdf, dt_pdf, "transport_", output_interval_pdf, up_dir, "Cu", &BTCs_Upscaled);

    BTCs_Upscaled.write(joinPath(up_dir, "BTC_Upscaled.csv"));
    BTCs_Upscaled.derivative().write(joinPath(up_dir, "BTC_Upscaled_derivative.csv"));

    g_u.writeNamedVTI_Auto("C", joinPath(up_dir, "Cu.vti"));
    g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell, joinPath(up_dir, "Cu.txt"));

    if (rank == 0) {
        std::cout << "\nMixing PDF simulation complete!\n";
        std::cout << "All outputs saved to: " << run_dir << std::endl;
    }

    return 0;
}
