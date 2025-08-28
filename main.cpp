// main.cpp
#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"
#include <cmath>
#include <iostream>
#include "grid.h"
#include "TimeSeries.h"

static inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx) { return j*nx + i; }
static inline bool on_bc(PetscInt i, PetscInt j, PetscInt nx, PetscInt ny, double x, double y, double Lx, double Ly) {
    const double eps = 1e-14;
    return (i==0) || (j==0) || (std::abs(x - Lx) < eps) || (std::abs(y - Ly) < eps);
}

int main(int argc, char** argv) {
    PETScInit petsc(argc, argv);

    Grid2D g(600,200,3,1);
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;
    g.makeGaussianFieldSGS("K_normal_score",0.5,0.1,10);
    PetscTime(&t_asm0);
    PetscTime(&t_total0);
    g.writeNamedMatrix("K_normal_score",Grid2D::ArrayKind::Cell, "K_normal_score.txt");
    g.writeNamedVTI("K_normal_score",Grid2D::ArrayKind::Cell, "NormalScore.vti");
    g.createExponentialField("K_normal_score",1,0,"K");
    PetscTime(&t_asm1);
    PetscTime(&t_solve0);
    g.DarcySolve(1,0,"K","K");
    std::cout<<"Darcy solved ... " <<std::endl;
    g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, "K.vti");
    PetscTime(&t_solve1);
    g.computeMassBalanceError("MassBalanceError");
    g.writeNamedMatrix("MassBalanceError",Grid2D::ArrayKind::Cell ,"error.txt");
    g.writeNamedVTI("head",Grid2D::ArrayKind::Cell, "Head.vti");
    g.writeNamedMatrix("head",Grid2D::ArrayKind::Cell, "Head.txt");
    g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, "qx.vti");
    g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, "qy.vti");

    TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx",Grid2D::ArrayKind::Fx);
    TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
    g.assignFromTimeSeries(QxNormalScores,"qx_normal_score",Grid2D::ArrayKind::Fx);
    TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score",Grid2D::ArrayKind::Fx, Grid2D::DerivDir::X, 10000, 0.05);
    curvture.writefile("2nd_deriv.txt");
    std::cout<<"Sampling points for derivative ..."<<std::endl;
    double dt_optimal = 0.5*g.dx()/g.fieldMinMax("qx",Grid2D::ArrayKind::Fx).second;
    std::cout<<"Optimal Time-Step: " << dt_optimal<<std::endl;
    g.assignConstant("C",Grid2D::ArrayKind::Cell, 0);
    g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, "qx.txt");
    g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, "qy.txt");
    g.writeNamedMatrix("K", Grid2D::ArrayKind::Cell, "K.txt");

    g.SetVal("diffusion", 0.0);
    g.SetVal("porosity", 1);
    g.SetVal("c_left", 1.0);
    g.SolveTransport(10,std::min(dt_optimal,0.5/10.0), "transport_",50);
    g.writeNamedVTI_Auto("C","C.vti");
    PetscTime(&t_total1);


    // ----- Report -----
    int rank=0; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank==0) {
        std::cout << "Assembly time: " << (t_asm1 - t_asm0) << " s\n";
        std::cout << "Solve time:    " << (t_solve1 - t_solve0) << " s\n";
        std::cout << "Total time:    " << (t_total1 - t_total0) << " s\n";
    }
    return 0;
}
