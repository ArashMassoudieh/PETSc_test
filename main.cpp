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
#include "grid.h"
#include "TimeSeries.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx) { return j*nx + i; }

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

int main(int argc, char** argv) {
    PETScInit petsc(argc, argv);

    // Set output directory (can be modified via command line argument)
    std::string output_dir = "/home/arash/Projects/UpscalingResults";


    int nx = 300;   // Very coarse for debugging
    int nu = 100;    // Very coarse for debugging
    int ny = 100;    // Very coarse for debugging
    double Lx = 3.0;
    double Ly = 1.0;

    Grid2D g(nx, ny, Lx, Ly);
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

    g.makeGaussianFieldSGS("K_normal_score", 1, 0.1, 10);

    PetscTime(&t_asm0);
    PetscTime(&t_total0);

    g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(output_dir, "K_normal_score.txt"));
    g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(output_dir, "NormalScore.vti"));
    g.createExponentialField("K_normal_score", 1, 0, "K");

    PetscTime(&t_asm1);
    PetscTime(&t_solve0);

    g.DarcySolve(1, 0, "K", "K");
    std::cout << "Darcy solved ... " << std::endl;
    g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(output_dir, "K.vti"));

    PetscTime(&t_solve1);

    g.computeMassBalanceError("MassBalanceError");
    g.writeNamedMatrix("MassBalanceError", Grid2D::ArrayKind::Cell, joinPath(output_dir, "error.txt"));
    g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(output_dir, "Head.vti"));
    g.writeNamedMatrix("head", Grid2D::ArrayKind::Cell, joinPath(output_dir, "Head.txt"));
    g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(output_dir, "qx.vti"));
    g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(output_dir, "qy.vti"));

    TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
    TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
    g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);

    std::cout << "Sampling points for derivative ..." << std::endl;
    TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score", Grid2D::ArrayKind::Fx, Grid2D::DerivDir::X, 10000, 0.05);
    curvture.writefile(joinPath(output_dir, "2nd_deriv.txt"));

    TimeSeries<double> diff_purturbed = g.sampleGaussianPerturbation("qx_normal_score", Grid2D::ArrayKind::Fx, 10000, 0.01, 0, PerturbDir::Radial);
    diff_purturbed.writefile(joinPath(output_dir, "Diff_purturbed.txt"));
    std::cout << "Diffusion purturbed correlation dx = 0.01: " << diff_purturbed.correlation_tc() << std::endl;

    TimeSeries<double> diff_purturbed2 = g.sampleGaussianPerturbation("qx_normal_score", Grid2D::ArrayKind::Fx, 10000, 0.005, 0, PerturbDir::Radial);
    diff_purturbed2.writefile(joinPath(output_dir, "Diff_purturbed2.txt"));
    std::cout << "Diffusion purturbed correlation dx = 0.005: " << diff_purturbed2.correlation_tc() << std::endl;

    double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;
    std::cout << "Optimal Time-Step: " << dt_optimal << std::endl;

    g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
    g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(output_dir, "qx.txt"));
    g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(output_dir, "qy.txt"));
    g.writeNamedMatrix("K", Grid2D::ArrayKind::Cell, joinPath(output_dir, "K.txt"));

    g.SetVal("diffusion", 0.01);
    g.SetVal("porosity", 1);
    g.SetVal("c_left", 1.0);
    std::vector<double> xLocations{0.5, 1.5, 2.5};
    g.setBTCLocations(xLocations);
    TimeSeriesSet<double> BTCs_FineScaled;
    g.SolveTransport(10, std::min(dt_optimal, 0.5/10.0), "transport_", 50, output_dir, "C", &BTCs_FineScaled);
    g.writeNamedVTI_Auto("C", joinPath(output_dir, "C.vti"));
    BTCs_FineScaled.write(joinPath(output_dir, "BTC_FileScaled.csv"));
    PetscTime(&t_total1);

    Pathway path(0);
    path.addParticle(0.0, 0.5, 0.0);  // Start at left boundary

    // Track particle with dx_step = 0.01
    path.trackParticle(&g, 0.01);

    // Check results
    std::cout << "Pathway length: " << path.pathLength() << std::endl;
    std::cout << "Travel time: " << path.duration() << std::endl;
    std::cout << "Exit position: (" << path.last().x() << ", "
              << path.last().y() << ")" << std::endl;

    // Save pathway
    path.writeVTK(joinPath(output_dir, "pathway.vtk"));

    PathwaySet pathways;

    // Initialize 1000 pathways with flux weighting
    pathways.Initialize(1000, PathwaySet::Weighting::FluxWeighted, &g);

    // Track all pathways
    pathways.trackAllPathways(&g, 0.01);  // dx_step = 0.01

    // Get statistics
    std::cout << "Mean path length: " << pathways.meanPathLength() << std::endl;
    std::cout << "Mean travel time: " << pathways.meanTravelTime() << std::endl;
    std::cout << "Completed pathways: " << pathways.countActive() << std::endl;

    // Write output
    pathways.writeToFile(joinPath(output_dir,"pathway_summary.txt"));
    pathways.writeCombinedVTK(joinPath(output_dir,"all_pathways.vtk"));

    // ----- Report -----
    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Assembly time: " << (t_asm1 - t_asm0) << " s\n";
        std::cout << "Solve time:    " << (t_solve1 - t_solve0) << " s\n";
        std::cout << "Total time:    " << (t_total1 - t_total0) << " s\n";
        std::cout << "All outputs saved to: " << output_dir << std::endl;
    }

    // ====================================================================
    // STEP 1: Extract velocity CDF from the flow field
    // ====================================================================
    std::cout << "\n=== Computing velocity CDF from qx field ===" << std::endl;

    TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100);
    qx_inverse_cdf = qx_inverse_cdf.make_uniform(0.01);
    qx_inverse_cdf.writefile(joinPath(output_dir, "qx_inverse_cdf.txt"));

    std::cout << "qx inverse CDF extracted with " << qx_inverse_cdf.size() << " points" << std::endl;
    std::cout << "Velocity range: [" << qx_inverse_cdf.interpol(0.0)
              << ", " << qx_inverse_cdf.interpol(1.0) << "]" << std::endl;

    // ====================================================================
    // STEP 2: Create grid for mixing PDF in (x,u) space
    // ====================================================================
    std::cout << "\n=== Creating Mixing PDF Grid ===" << std::endl;

    Grid2D g_u(nx, ny, Lx, Ly);

    std::cout << "Grid created: " << nx << " x " << nu
              << " (x is spatial, y is u-space)" << std::endl;

    // ====================================================================
    // STEP 3: Initialize flux arrays (IMPORTANT!)
    // ====================================================================
    std::cout << "\n=== Initializing flux arrays ===" << std::endl;

    // Calculate expected sizes for flux arrays
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
        double u = static_cast<double>(j) / (nu - 1);
        double v_at_u = qx_inverse_cdf.interpol(u);

        // Set qx for all x-positions at this u-level
        for (int i = 0; i < nx + 1; ++i) {
            int idx = j * (nx + 1) + i;
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

    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(output_dir, "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(output_dir, "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(output_dir, "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(output_dir, "qy_u_initial.txt"));

    std::cout << "Initial velocity fields written" << std::endl;

    // ====================================================================
    // STEP 6: Set mixing parameters
    // ====================================================================
    std::cout << "\n=== Setting mixing parameters ===" << std::endl;

    double lc = 0.05;        // Correlation length scale
    double lambda_x = 1; // Transverse dispersion length in x
    double lambda_y = 0.1; // Transverse dispersion length in y

    g_u.setMixingParams(lc, lambda_x, lambda_y);
    g_u.SetVal("diffusion", 0.01);  // Molecular diffusion
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);       // Uniform PDF at inlet

    // ====================================================================
    // STEP 7: Solve mixing PDF equation
    // ====================================================================
    std::cout << "\n=== Solving Mixing PDF Equation ===" << std::endl;

    double t_end_pdf = 10;
    double dt_pdf = dt_optimal;
    int output_interval_pdf = 50;
    std::cout << "\n=== Computing effective D_y ===" << std::endl;
    g_u.computeMixingDiffusionCoefficient();
    std::cout << "\n=== Writing D_y into vti ===" << std::endl;
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(output_dir, "D_y.vti"));
    std::cout << "\n=== Solving ===" << std::endl;
    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
    TimeSeriesSet<double> BTCs_Upscaled;
    g_u.setBTCLocations(xLocations);
    g_u.SolveTransport(t_end_pdf, dt_pdf, "transport_", output_interval_pdf, output_dir, "Cu",&BTCs_Upscaled);
    BTCs_Upscaled.write(joinPath(output_dir, "BTC_Upscaled.csv"));

    // ====================================================================
    // STEP 8: Save final PDF distribution
    // ====================================================================
    std::cout << "\n=== Saving final results ===" << std::endl;

    g_u.writeNamedVTI_Auto("C", joinPath(output_dir, "Cu.vti"));
    g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell,
                         joinPath(output_dir, "Cu.txt"));

    std::cout << "\nMixing PDF simulation complete!" << std::endl;

    PetscTime(&t_total1);
    return 0;
}
