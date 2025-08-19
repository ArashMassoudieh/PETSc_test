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

    Grid2D g(1000,1000,1,1);

    g.makeGaussianFieldSGS("K_normal_score",0.5,0.01,10);
    g.writeNamedMatrix("K_normal_score",Grid2D::ArrayKind::Cell, "K_normal_score.txt");
    g.writeNamedVTI_Auto("K_normal_score","NormalScore.vti");
    // ----- Problem size & coefficients -----
    PetscInt  nx = 1000, ny = 1000;
    PetscReal Lx = 1.0,  Ly = 1.0;
    PetscReal Kx = 1.0,  Ky = 1.0;

    (void)PetscOptionsGetInt (nullptr,nullptr,"-nx",&nx,nullptr);
    (void)PetscOptionsGetInt (nullptr,nullptr,"-ny",&ny,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Lx",&Lx,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Ly",&Ly,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Kx",&Kx,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Ky",&Ky,nullptr);

    const PetscInt N = nx * ny;
    const double dx = Lx / (nx - 1);
    const double dy = Ly / (ny - 1);
    const double ax = Kx / (dx*dx);
    const double ay = Ky / (dy*dy);

    // ----- Timers -----
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;
    PetscTime(&t_total0);

    // ----- Allocate matrix and vectors -----
    PETScMatrix A(N, N, /*nzPerRow=*/5);
    PETScVector b(N);  b.set(0.0);

    // ----- Assembly (eliminate BCs during assembly) -----
    PetscTime(&t_asm0);

    PetscInt Istart=0, Iend=0;
    A.ownershipRange(Istart, Iend);

    for (PetscInt I = Istart; I < Iend; ++I) {
        const PetscInt j = I / nx;
        const PetscInt i = I - j*nx;
        const double x = i * dx;
        const double y = j * dy;

        // Dirichlet: h=1 at x=0; h=0 on x=Lx, y=0, y=Ly
        const bool isBC = on_bc(i,j,nx,ny,x,y,Lx,Ly);
        if (isBC) {
            const double hbc = (i==0 ? 1.0 : 0.0);
            A.setValue(I, I, 1.0, INSERT_VALUES);
            b.setValue(I, hbc, INSERT_VALUES);
            continue;
        }

        // Interior: build stencil, but eliminate neighbors that are BC by:
        //  - skipping their off-diagonal entry
        //  - adding that contribution to b (i.e., move to RHS)
        PetscInt    cols[5];
        PetscScalar vals[5];
        int n = 0;

        // West neighbor
        {
            const PetscInt in = i-1, jn = j;
            const double   xn = in * dx, yn = jn * dy;
            const PetscScalar c = -ax;
            if (on_bc(in,jn,nx,ny,xn,yn,Lx,Ly)) {
                const double hbc = (in==0 ? 1.0 : 0.0);
                b.setValue(I, -c * hbc, ADD_VALUES); // move A(I,nb)*hbc to RHS
            } else {
                cols[n] = idx(in,jn,nx); vals[n] = c; ++n;
            }
        }

        // East neighbor
        {
            const PetscInt in = i+1, jn = j;
            const double   xn = in * dx, yn = jn * dy;
            const PetscScalar c = -ax;
            if (on_bc(in,jn,nx,ny,xn,yn,Lx,Ly)) {
                const double hbc = (in==0 ? 1.0 : 0.0);
                b.setValue(I, -c * hbc, ADD_VALUES);
            } else {
                cols[n] = idx(in,jn,nx); vals[n] = c; ++n;
            }
        }

        // South neighbor
        {
            const PetscInt in = i, jn = j-1;
            const double   xn = in * dx, yn = jn * dy;
            const PetscScalar c = -ay;
            if (on_bc(in,jn,nx,ny,xn,yn,Lx,Ly)) {
                const double hbc = (in==0 ? 1.0 : 0.0);
                b.setValue(I, -c * hbc, ADD_VALUES);
            } else {
                cols[n] = idx(in,jn,nx); vals[n] = c; ++n;
            }
        }

        // North neighbor
        {
            const PetscInt in = i, jn = j+1;
            const double   xn = in * dx, yn = jn * dy;
            const PetscScalar c = -ay;
            if (on_bc(in,jn,nx,ny,xn,yn,Lx,Ly)) {
                const double hbc = (in==0 ? 1.0 : 0.0);
                b.setValue(I, -c * hbc, ADD_VALUES);
            } else {
                cols[n] = idx(in,jn,nx); vals[n] = c; ++n;
            }
        }

        // Center diagonal
        cols[n] = I; vals[n] = 2.0*ax + 2.0*ay; ++n;

        A.setValues(1, &I, n, cols, vals, INSERT_VALUES);
    }

    A.assemble();
    b.assemble();

    // Assert properties (now truly symmetric/SPD)
    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(A.raw(), MAT_SPD, PETSC_TRUE));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(A.raw(), MAT_SYMMETRIC, PETSC_TRUE));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(A.raw(), MAT_SYMMETRY_ETERNAL, PETSC_TRUE));

    PetscTime(&t_asm1);

    // ----- Solve -----
    PetscTime(&t_solve0);
    PETScVector sol = A / b;  // our implementation zeros the initial guess internally
    PetscTime(&t_solve1);

    // ----- Residual -----
    PETScVector r(N);
    A.multiply(sol, r);                                           // r = A*sol
    PetscCallAbort(PETSC_COMM_WORLD, VecAYPX(r.raw(), -1.0, b.raw())); // r = r - b
    PetscReal rnorm = r.norm(NORM_2);

    // Stop total timing before host reads
    PetscTime(&t_total1);

    // ----- Center value (host read) -----
    PetscScalar h_center = 0.0;
    if ((nx % 2 == 1) && (ny % 2 == 1)) {
        const PetscInt ic = nx / 2, jc = ny / 2;
        const PetscInt I  = jc * nx + ic;
        PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(sol.raw(), 1, &I, &h_center));
    } else {
        const PetscInt i0 = nx/2 - 1, j0 = ny/2 - 1;
        const PetscInt ids[4] = { j0*nx + i0, j0*nx + (i0+1), (j0+1)*nx + i0, (j0+1)*nx + (i0+1) };
        PetscScalar vals[4];
        PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(sol.raw(), 4, ids, vals));
        h_center = 0.25*(vals[0] + vals[1] + vals[2] + vals[3]);
    }
    const double x_center = (nx % 2 == 1) ? (nx/2) * (Lx/(nx-1)) : 0.5 * Lx;
    const double y_center = (ny % 2 == 1) ? (ny/2) * (Ly/(ny-1)) : 0.5 * Ly;

    // ----- Report -----
    int rank=0; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank==0) {
        std::cout << "Solved Darcy/Laplace: N=" << N
                  << "  residual ||Ax-b||_2 = " << rnorm << "\n";
        std::cout << "Center value h(" << x_center << "," << y_center << ") = "
                  << h_center << "\n";
        std::cout << "Assembly time: " << (t_asm1 - t_asm0) << " s\n";
        std::cout << "Solve time:    " << (t_solve1 - t_solve0) << " s\n";
        std::cout << "Total time:    " << (t_total1 - t_total0) << " s\n";
    }
    return 0;
}
