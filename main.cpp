// main.cpp
#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"
#include <cmath>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <ctime>
#include <mpi.h>

#include "grid.h"
#include "TimeSeries.h"

static inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx) { return j*nx + i; }
static inline bool on_bc(PetscInt i, PetscInt j, PetscInt nx, PetscInt ny, double x, double y, double Lx, double Ly) {
    const double eps = 1e-14;
    return (i==0) || (j==0) || (std::abs(x - Lx) < eps) || (std::abs(y - Ly) < eps);
}

// ensure directory exists
static void ensureDirectory(const std::string& path) {
    if (mkdir(path.c_str(), 0755) != 0) {
        if (errno == EEXIST) {
            return; // already exists
        } else {
            perror(("Error creating directory " + path).c_str());
            exit(EXIT_FAILURE);
        }
    }
}

// timestamp for run folder
static std::string makeTimestamp() {
    std::time_t now = std::time(nullptr);
    std::tm tm_now;
    localtime_r(&now, &tm_now);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d_%H-%M-%S", &tm_now);
    return std::string(buf);
}

int main(int argc, char** argv) {
    PETScInit petsc(argc, argv);

    // --- results folder under project path ---
    std::string baseDir = "/mnt/3rd900/Projects/PETSc_test/Results/";
    ensureDirectory(baseDir);
    std::string timestamp = makeTimestamp();
    std::string resultsDir = baseDir + "run_" + timestamp + "/";
    ensureDirectory(resultsDir);

    Grid2D g(300,100,3,1);
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

    g.makeGaussianFieldSGS("K_normal_score",0.5,0.1,10);
    PetscTime(&t_asm0);
    PetscTime(&t_total0);
    g.writeNamedMatrix("K_normal_score",Grid2D::ArrayKind::Cell, resultsDir + "K_normal_score.txt");
    g.writeNamedVTI("K_normal_score",Grid2D::ArrayKind::Cell, resultsDir + "NormalScore.vti");
    g.createExponentialField("K_normal_score",1,0,"K");
    PetscTime(&t_asm1);

    PetscTime(&t_solve0);
    g.DarcySolve(1,0,"K","K");
    std::cout<<"Darcy solved ... " <<std::endl;
    g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, resultsDir + "K.vti");
    PetscTime(&t_solve1);

    g.computeMassBalanceError("MassBalanceError");
    g.writeNamedMatrix("MassBalanceError",Grid2D::ArrayKind::Cell ,resultsDir + "error.txt");
    g.writeNamedVTI("head",Grid2D::ArrayKind::Cell, resultsDir + "Head.vti");
    g.writeNamedMatrix("head",Grid2D::ArrayKind::Cell, resultsDir + "Head.txt");
    g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, resultsDir + "qx.vti");
    g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, resultsDir + "qy.vti");

    TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx",Grid2D::ArrayKind::Fx);
    TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
    g.assignFromTimeSeries(QxNormalScores,"qx_normal_score",Grid2D::ArrayKind::Fx);

    std::cout<<"Sampling points for derivative ..."<<std::endl;
    g.writeNamedVTI("qx_normal_score",Grid2D::ArrayKind::Fx,resultsDir + "qx_normal_scores.vti","qx normal score");

    TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score",Grid2D::ArrayKind::Fx, Grid2D::DerivDir::X, 10000, 0.05);
    curvture.writefile(resultsDir + "2nd_deriv.txt");

    TimeSeries<double> x_correlations;
    for (double logdx = -5; logdx<=-1; logdx+=0.1) {
        TimeSeries<double> diff_purturbed = g.sampleGaussianPerturbation("qx_normal_score",Grid2D::ArrayKind::Fx,10000,exp(logdx),0, PerturbDir::XOnly);
        double correlation = diff_purturbed.correlation_tc();
        std::cout << "Diffusion purturbed correlation dx = " << exp(logdx) << "=" << correlation << std::endl;
        x_correlations.append(exp(logdx),correlation);
    }
    x_correlations.writefile(resultsDir + "Correlations_x.csv");

    TimeSeries<double> diff_purturbed2 = g.sampleGaussianPerturbation("qx_normal_score",Grid2D::ArrayKind::Fx,10000,0.005,0, PerturbDir::Radial);
    diff_purturbed2.writefile(resultsDir + "Diff_purturbed2.txt");
    std::cout << "Diffusion purturbed correlation dx = 0.005: " << diff_purturbed2.correlation_tc() << std::endl;

    double dt_optimal = 0.5*g.dx()/g.fieldMinMax("qx",Grid2D::ArrayKind::Fx).second;
    std::cout<<"Optimal Time-Step: " << dt_optimal<<std::endl;

    g.assignConstant("C",Grid2D::ArrayKind::Cell, 0);
    g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, resultsDir + "qx.txt");
    g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, resultsDir + "qy.txt");
    g.writeNamedMatrix("K", Grid2D::ArrayKind::Cell, resultsDir + "K.txt");

    g.SetVal("diffusion", 0.0);
    g.SetVal("porosity", 1);
    g.SetVal("c_left", 1.0);

    g.SolveTransport(10,std::min(dt_optimal,0.5/10.0), (resultsDir + "transport_").c_str(),50);
    g.writeNamedVTI_Auto("C",resultsDir + "C.vti");

    PetscTime(&t_total1);

    // ----- Report -----
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank==0) {
        std::cout << "Assembly time: " << (t_asm1 - t_asm0) << " s\n";
        std::cout << "Solve time: " << (t_solve1 - t_solve0) << " s\n";
        std::cout << "Total time: " << (t_total1 - t_total0) << " s\n";
        std::cout << "Results saved in: " << resultsDir << std::endl;
    }

    return 0;
}
