// main.cpp
#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"
#include <cmath>
#include <iostream>
#include "grid.h"

static inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx) { return j*nx + i; }
static inline bool on_bc(PetscInt i, PetscInt j, PetscInt nx, PetscInt ny, double x, double y, double Lx, double Ly) {
    const double eps = 1e-14;
    return (i==0) || (j==0) || (std::abs(x - Lx) < eps) || (std::abs(y - Ly) < eps);
}

int main(int argc, char** argv) {
    PETScInit petsc(argc, argv);

    Grid2D g(3000,1000,3,1);
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;
    g.makeGaussianFieldSGS("K_normal_score",0.5,0.1,10);
    PetscTime(&t_asm0);
    PetscTime(&t_total0);
    g.writeNamedMatrix("K_normal_score",Grid2D::ArrayKind::Cell, "K_normal_score.txt");
    g.writeNamedVTI_Auto("K_normal_score","NormalScore.vti");
    g.createExponentialField("K_normal_score",1,0,"K");
    PetscTime(&t_asm1);
    PetscTime(&t_solve0);
    g.DarcySolve(1,0,"K","K");
    g.writeNamedVTI_Auto("K", "K.vti");
    PetscTime(&t_solve1);
    g.writeNamedVTI_Auto("head", "Head.vti");
    g.writeNamedVTI_Auto("qx", "qx.vti");
    g.writeNamedVTI_Auto("qy", "qy.vti");
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
