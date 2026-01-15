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
// IMPORTANT: This is defined but not called in main() below.
// If output_dir doesn't exist, file writes may fail silently depending on your I/O code.
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
// IMPORTANT: Uses "/" when dir does not end with "/" or "\" (works fine on Linux; ok on Windows in many cases too).
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

int main(int argc, char** argv) {
    // PETSc + MPI init RAII wrapper (assumed)
    PETScInit petsc(argc, argv);

    // Set output directory (can be modified via command line argument)
    // IMPORTANT: output_dir is only defined if one of these macros is defined.
    // If none are defined, output_dir will be uninitialized/undefined behavior.
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

    // Domain / grid resolution for primary Darcy/transport solves
    // NOTE: comments say "very coarse for debugging" but nx=300 is not super coarse.
    int nx = 300;   // Very coarse for debugging
    int nu = 100;    // Very coarse for debugging  (used later for u-grid logic)
    int ny = 100;    // Very coarse for debugging
    double Lx = 3.0;
    double Ly = 1.0;
    double correlation_ls_x = 1;
    double correlation_ls_y = 0.1;
    double stdev = 1.0;
    double g_mean = 0;
    double Diffusion_coefficient = 0;

    // Primary spatial grid (x,y)
    Grid2D g(nx, ny, Lx, Ly);

    // PETSc timers (assembly/solve/total segments)
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

    // Generate heterogeneous field: Gaussian random field via SGS in "normal score" space
    g.makeGaussianFieldSGS("K_normal_score", correlation_ls_x, correlation_ls_y, 10);

    // Start timing (NOTE: order here makes assembly time include earlier work depending on what follows)
    PetscTime(&t_asm0);
    PetscTime(&t_total0);

    // Write raw normal-score permeability field
    g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(run_dir, "K_normal_score.txt"));
    g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(run_dir, "NormalScore.vti"));

    // Transform normal score to exponential field "K"
    g.createExponentialField("K_normal_score", stdev, g_mean, "K");

    // End "assembly" timer (name suggests assembly, but it currently covers field generation + writes too)
    PetscTime(&t_asm1);

    // Start solver timer
    PetscTime(&t_solve0);

    // Darcy solve using "K" (inputs/outputs depend on your Grid2D::DarcySolve signature)
    g.DarcySolve(Lx, 0, "K", "K");
    std::cout << "Darcy solved ... " << std::endl;

    // Save K field
    g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(run_dir, "K.vti"));

    // End solver timer
    PetscTime(&t_solve1);

    // Diagnostics: mass balance error field
    g.computeMassBalanceError("MassBalanceError");
    g.writeNamedMatrix("MassBalanceError", Grid2D::ArrayKind::Cell, joinPath(run_dir, "error.txt"));

    // Save head + flux fields
    g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(run_dir, "Head.vti"));
    g.writeNamedMatrix("head", Grid2D::ArrayKind::Cell, joinPath(run_dir, "Head.txt"));
    g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(run_dir, "qx.vti"));
    g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(run_dir, "qy.vti"));

    // Export qx to TimeSeries then convert to normal-score for statistical derivative sampling
    TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
    TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
    g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);

    g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(run_dir, "qx_normal_score.vti"));

    // Sample second derivative of qx_normal_score (curvature proxy) along X direction
    std::cout << "Sampling points for derivative ..." << std::endl;
    TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score", Grid2D::ArrayKind::Fx, Grid2D::DerivDir::X, 10000, 0.05);
    curvture.writefile(joinPath(run_dir, "2nd_deriv.txt"));

    // =============================================================================
    // Velocity autocorrelation analysis via Gaussian perturbation sampling
    // =============================================================================

    // Define range of separation distances (logarithmic spacing)
    double delta_min = 0.001;
    double delta_max = 0.2;
    int num_deltas = 30;
    int num_samples_per_delta = 10000;

    // TimeSeries to store correlation vs distance for each direction
    TimeSeries<double> corr_radial;
    TimeSeries<double> corr_x;
    TimeSeries<double> corr_y;

    std::cout << "Computing velocity autocorrelation (radial)..." << std::endl;
    for (int i = 0; i < num_deltas; ++i) {
        double exponent = static_cast<double>(i) / (num_deltas - 1);
        double delta = delta_min * std::pow(delta_max / delta_min, exponent);
        try {
            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                "qx_normal_score", Grid2D::ArrayKind::Fx,
                num_samples_per_delta, delta, 0, PerturbDir::Radial);
            double correlation = samples.correlation_tc();
            corr_radial.append(delta, correlation);
            std::cout << "  delta = " << delta << ", correlation = " << correlation << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed at delta = " << delta << ": " << e.what() << std::endl;
        }
    }
    corr_radial.writefile(joinPath(run_dir, "velocity_correlation_radial.txt"));
    double lambda_radial = corr_radial.fitExponentialDecay();
    std::cout << "Radial correlation length scale: " << lambda_radial << std::endl;

    std::cout << "Computing velocity autocorrelation (x-direction)..." << std::endl;
    for (int i = 0; i < num_deltas; ++i) {
        double exponent = static_cast<double>(i) / (num_deltas - 1);
        double delta = delta_min * std::pow(delta_max / delta_min, exponent);
        try {
            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                "qx_normal_score", Grid2D::ArrayKind::Fx,
                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
            double correlation = samples.correlation_tc();
            corr_x.append(delta, correlation);
            std::cout << "  delta = " << delta << ", correlation = " << correlation << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed at delta = " << delta << ": " << e.what() << std::endl;
        }
    }
    corr_x.writefile(joinPath(run_dir, "velocity_correlation_x.txt"));
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
            double correlation = samples.correlation_tc();
            corr_y.append(delta, correlation);
            std::cout << "  delta = " << delta << ", correlation = " << correlation << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed at delta = " << delta << ": " << e.what() << std::endl;
        }
    }
    corr_y.writefile(joinPath(run_dir, "velocity_correlation_y.txt"));
    double lambda_y_emp = corr_y.fitExponentialDecay();
    std::cout << "Y-direction correlation length scale: " << lambda_y_emp << std::endl;

    std::cout << "\n=== Velocity Correlation Length Scales ===" << std::endl;
    std::cout << "  Radial: " << lambda_radial << std::endl;
    std::cout << "  X-dir:  " << lambda_x_emp << std::endl;
    std::cout << "  Y-dir:  " << lambda_y_emp << std::endl;

    // CFL-like timestep estimate based on max qx (NOTE: if max qx ~ 0, dt becomes huge)
    double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;
    std::cout << "Optimal Time-Step: " << dt_optimal << std::endl;

    // Initialize transport concentration field
    g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

    // Save raw fields as matrices (for debugging / external scripts)
    g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(run_dir, "qx.txt"));
    g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(run_dir, "qy.txt"));
    g.writeNamedMatrix("K", Grid2D::ArrayKind::Cell, joinPath(run_dir, "K.txt"));

    // Transport parameters
    g.SetVal("diffusion", Diffusion_coefficient);
    g.SetVal("porosity", 1);
    g.SetVal("c_left", 1.0);

    // BTC locations (x positions) where breakthrough curves will be recorded
    std::vector<double> xLocations{0.5, 1.5, 2.5};
    g.setBTCLocations(xLocations);

    // Fine-scale transport solve (uses g’s spatial velocity field)
    TimeSeriesSet<double> BTCs_FineScaled;
    g.SolveTransport(10, std::min(dt_optimal, 0.5/10.0), "transport_", 50, run_dir, "C", &BTCs_FineScaled);

    // Save final concentration
    g.writeNamedVTI_Auto("C", joinPath(run_dir, "C.vti"));

    // Save BTCs
    BTCs_FineScaled.write(joinPath(run_dir, "BTC_FineScaled.csv"));
    BTCs_FineScaled.derivative().write(joinPath(run_dir, "BTC_FineScaled_derivative.csv"));

    // Total timing checkpoint (NOTE: later you overwrite t_total1 again)
    PetscTime(&t_total1);

    // --------------------------------------------------------------------
    // Single particle tracking (diagnostic)
    // --------------------------------------------------------------------
    Pathway path(0);
    path.addParticle(0.0, 0.5, 0.0);  // Start at left boundary (x=0, y=0.5)

    // Track particle with dx_step = 0.01 (step size in x)
    path.trackParticle(&g, 0.01);

    // Check results
    std::cout << "Pathway length: " << path.pathLength() << std::endl;
    std::cout << "Travel time: " << path.duration() << std::endl;
    std::cout << "Exit position: (" << path.last().x() << ", "
              << path.last().y() << ")" << std::endl;

    // Save pathway polyline / points
    path.writeVTK(joinPath(run_dir, "pathway.vtk"));

    // --------------------------------------------------------------------
    // Many particle tracking (PathwaySet) for statistics + correlation
    // --------------------------------------------------------------------
    PathwaySet pathways;

    // Initialize 1000 pathways with flux weighting (inlet sampling weighted by flux)
    pathways.Initialize(1000, PathwaySet::Weighting::FluxWeighted, &g);

    // Track all pathways with dx_step = 0.01
    pathways.trackAllPathways(&g, 0.01);

    // ====================================================================
    // Compute velocity correlation as a function of separation distance
    // ====================================================================
    std::cout << "\n=== Computing velocity correlation function ===" << std::endl;

    // Define range of separation distances (logarithmic spacing)
    double Delta_x_min = 0.001;
    double Delta_x_max = 0.5;
    int num_Delta_x = 30;  // Number of separation distances to sample
    int num_samples_per_Delta_x = 10000;  // Number of particle pairs per distance

    TimeSeries<double> qx_correlation;

    for (int i = 0; i < num_Delta_x; ++i) {
        // Logarithmic spacing: Delta_x = Delta_x_min * (Delta_x_max/Delta_x_min)^(i/(num-1))
        double exponent = static_cast<double>(i) / (num_Delta_x - 1);
        double Delta_x = Delta_x_min * std::pow(Delta_x_max / Delta_x_min, exponent);

        try {
            // Sample particle pairs at this separation distance
            // particle_pairs[0] contains particles at x
            // particle_pairs[1] contains particles at x+Delta_x
            // NOTE: this creates/returns a PathwaySet containing two "groups" (index 0 and 1) by your convention
            PathwaySet particle_pairs = pathways.sampleParticlePairs(Delta_x, num_samples_per_Delta_x);

            // Calculate correlation between qx values at x and x+Delta_x
            // NOTE: calculateCorrelation(0,1,"qx") assumes both groups have matched pairing logic
            double correlation = particle_pairs.calculateCorrelation(0, 1, "qx");

            // Store in TimeSeries: (Delta_x, correlation)
            qx_correlation.append(Delta_x, correlation);

            std::cout << "  Delta_x = " << Delta_x
                      << ", correlation = " << correlation << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to compute correlation at Delta_x = "
                      << Delta_x << ": " << e.what() << std::endl;
        }
    }

    // Save correlation function
    qx_correlation.writefile(joinPath(run_dir, "qx_correlation_vs_distance.txt"));

    // Fit exponential decay to estimate correlation length scale (used later conceptually)
    double advection_correlation_length_scale = qx_correlation.fitExponentialDecay();

    std::cout << "Velocity correlation function saved with "
              << qx_correlation.size() << " points" << std::endl;

    std::cout << "Velocity correlation length scale = " << advection_correlation_length_scale << std::endl;

    // PathwaySet statistics
    std::cout << "Mean path length: " << pathways.meanPathLength() << std::endl;
    std::cout << "Mean travel time: " << pathways.meanTravelTime() << std::endl;

    // NOTE: "Completed pathways" label, but function is countActive() — check semantics in your class.
    std::cout << "Completed pathways: " << pathways.countActive() << std::endl;

    // Write pathway summary and combined VTK
    pathways.writeToFile(joinPath(run_dir,"pathway_summary.txt"));
    pathways.writeCombinedVTK(joinPath(run_dir,"all_pathways.vtk"));

    // ----- Report -----
    // NOTE: timing here uses t_total1 that was set earlier (after fine-scale transport).
    // Later you call PetscTime(&t_total1) again at the end.
    // (rank already defined above; keep code as-is otherwise)
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Assembly time: " << (t_asm1 - t_asm0) << " s\n";
        std::cout << "Solve time:    " << (t_solve1 - t_solve0) << " s\n";
        std::cout << "Total time:    " << (t_total1 - t_total0) << " s\n";
        std::cout << "All outputs saved to: " << run_dir << std::endl;
    }

    // ====================================================================
    // STEP 1: Extract velocity CDF from the flow field
    // ====================================================================
    std::cout << "\n=== Computing velocity CDF from qx field ===" << std::endl;

    // Extract inverse CDF (quantile function) from qx field with 100 bins/points
    TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100);

    // Make uniform in probability axis at dp=0.01 (common for stable interpolation)
    qx_inverse_cdf = qx_inverse_cdf.make_uniform(0.01);

    // Save inverse CDF
    qx_inverse_cdf.writefile(joinPath(run_dir, "qx_inverse_cdf.txt"));

    // Print velocity range using quantiles p=0 and p=1
    std::cout << "qx inverse CDF extracted with " << qx_inverse_cdf.size() << " points" << std::endl;
    std::cout << "Velocity range: [" << qx_inverse_cdf.interpol(0.0)
              << ", " << qx_inverse_cdf.interpol(1.0) << "]" << std::endl;

    // ====================================================================
    // STEP 2: Create grid for mixing PDF in (x,u) space
    // ====================================================================
    std::cout << "\n=== Creating Mixing PDF Grid ===" << std::endl;

    // IMPORTANT: This creates Grid2D(nx, ny, Lx, Ly).
    // But your comment says (x,u) space and prints nx x nu below.
    // That means ny is being used as the second dimension here, not nu.
    Grid2D g_u(nx, ny, Lx, Ly);

    std::cout << "Grid created: " << nx << " x " << nu
              << " (x is spatial, y is u-space)" << std::endl;

    // ====================================================================
    // STEP 3: Initialize flux arrays (IMPORTANT!)
    // ====================================================================
    std::cout << "\n=== Initializing flux arrays ===" << std::endl;

    // Expected sizes for flux arrays for a staggered grid:
    // Fx size ~ (nx+1)*nu ; Fy size ~ nx*(nu+1)
    // NOTE: these sizes are computed with nu, but g_u was created with ny above.
    int qx_size = nu * (nx + 1);
    int qy_size = (nu + 1) * nx;

    // Initialize flux arrays with zeros
    auto& qx_u = g_u.flux("qx");
    auto& qy_u = g_u.flux("qy");

    qx_u.resize(qx_size, 0.0);
    qy_u.resize(qy_size, 0.0);

    std::cout << "Flux arrays initialized: qx(" << qx_size << "), qy(" << qy_size << ")" << std::endl;

    // ====================================================================
    // STEP 4: Assign velocity field v(u) from CDF
    // ====================================================================
    std::cout << "\n=== Setting velocity field v(u) from CDF ===" << std::endl;

    // For each u-level, assign velocity from the inverse CDF
    for (int j = 0; j < nu; ++j) {
        double u = static_cast<double>(j) / (nu - 1);   // u in [0,1]
        double v_at_u = qx_inverse_cdf.interpol(u);      // v(u) = quantile(u)

        // Set qx for all x-positions at this u-level
        // NOTE: uses nx+1, which matches staggered Fx indexing in x direction
        for (int i = 0; i < nx + 1; ++i) {
            int idx = j * (nx + 1) + i; // NOTE: local var "idx" shadows function idx(...) above
            if (idx >= qx_size) {
                std::cerr << "ERROR: qx index out of bounds: " << idx << " >= " << qx_size << std::endl;
                return 1;
            }
            qx_u[idx] = v_at_u;
        }
    }

    // qy is already set to zero from resize
    std::cout << "Velocity field v(u) assigned to grid" << std::endl;

    // ====================================================================
    // STEP 5: Write initial fields for verification
    // ====================================================================
    std::cout << "\n=== Writing initial fields for verification ===" << std::endl;

    // Save velocity fields on (x,u) grid (as VTI + text)
    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(run_dir, "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(run_dir, "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(run_dir, "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(run_dir, "qy_u_initial.txt"));

    std::cout << "Initial velocity fields written" << std::endl;

    // ====================================================================
    // STEP 6: Set mixing parameters
    // ====================================================================
    std::cout << "\n=== Setting mixing parameters ===" << std::endl;

    double lc = advection_correlation_length_scale;        // Correlation length scale (NOTE: you computed one above too) ~ 0.35
    double lambda_x = lambda_x_emp;     // Transverse dispersion length in x direction ~ 1.0
    double lambda_y = lambda_y_emp;   // Transverse dispersion length in y direction ~ 0.10

    // Set mixing model parameters used by your upscaled PDE
    g_u.setMixingParams(lc, lambda_x, lambda_y);

    // Molecular diffusion + porosity + inlet PDF/BC
    g_u.SetVal("diffusion", Diffusion_coefficient);  // Molecular diffusion
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);       // Uniform PDF at inlet

    // ====================================================================
    // STEP 7: Solve mixing PDF equation
    // ====================================================================
    std::cout << "\n=== Solving Mixing PDF Equation ===" << std::endl;

    double t_end_pdf = 10;
    double dt_pdf = dt_optimal;          // NOTE: dt_optimal came from fine-scale grid, reused here
    int output_interval_pdf = 50;

    // Compute effective mixing diffusion coefficient (writes into "D_y" field by your convention)
    std::cout << "\n=== Computing effective D_y ===" << std::endl;
    g_u.computeMixingDiffusionCoefficient();

    // Save D_y for verification
    std::cout << "\n=== Writing D_y into vti ===" << std::endl;
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(run_dir, "D_y.vti"));

    // Solve transport on (x,u) grid
    std::cout << "\n=== Solving ===" << std::endl;

    // Initialize scalar (NOTE: you later call SolveTransport(...,"Cu",...) not "C")
    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

    TimeSeriesSet<double> BTCs_Upscaled;
    g_u.setBTCLocations(xLocations);

    // IMPORTANT: last argument "Cu" is the field name passed to SolveTransport
    // But above you initialized "C", and later you write "C". Ensure your Grid2D expects this.
    g_u.SolveTransport(t_end_pdf, dt_pdf, "transport_", output_interval_pdf, run_dir, "Cu",&BTCs_Upscaled);

    // Save upscaled BTCs
    BTCs_Upscaled.write(joinPath(run_dir, "BTC_Upscaled.csv"));
    BTCs_Upscaled.derivative().write(joinPath(run_dir, "BTC_Upscaled_derivative.csv"));

    // ====================================================================
    // STEP 8: Save final PDF distribution
    // ====================================================================
    std::cout << "\n=== Saving final results ===" << std::endl;

    // NOTE: filenames say Cu.*, but field written is "C"
    g_u.writeNamedVTI_Auto("C", joinPath(run_dir, "Cu.vti"));
    g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell,
                         joinPath(run_dir, "Cu.txt"));

    std::cout << "\nMixing PDF simulation complete!" << std::endl;

    // NOTE: overwrites t_total1 again (timing report earlier won’t reflect this)
    PetscTime(&t_total1);
    return 0;
}
