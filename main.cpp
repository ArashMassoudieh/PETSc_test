#include <petscksp.h>
#include <petscdmda.h>
#include <vector>

static inline PetscScalar g_dirichlet(PetscReal x, PetscReal, PetscReal, PetscReal) {
    // h = 1 on x=0, 0 on other sides
    const PetscReal eps = 1e-14;
    return (x < eps) ? 1.0 : 0.0;
}

int main(int argc, char** argv) {
    PetscCall(PetscInitialize(&argc,&argv,nullptr,nullptr));

    // --- timers ---
    PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;
    PetscCall(PetscTime(&t_total0));

    // Problem size and coefficients (overridable via -nx -ny -Kx -Ky -Lx -Ly)
    PetscInt  nx = 1000, ny = 1000;
    PetscReal Lx = 1.0, Ly = 1.0, Kx = 1.0, Ky = 1.0;
    (void)PetscOptionsGetInt (nullptr,nullptr,"-nx",&nx,nullptr);
    (void)PetscOptionsGetInt (nullptr,nullptr,"-ny",&ny,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Lx",&Lx,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Ly",&Ly,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Kx",&Kx,nullptr);
    (void)PetscOptionsGetReal(nullptr,nullptr,"-Ky",&Ky,nullptr);

    const PetscReal dx = Lx/(nx-1), dy = Ly/(ny-1);

    // DMDA: 2D, 5-point STAR stencil, dof=1
    DM da;
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                           DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DMDA_STENCIL_STAR,
                           nx, ny, PETSC_DECIDE, PETSC_DECIDE,
                           1, 1, nullptr, nullptr, &da));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));

    // Optional: uncomment to force GPU defaults without flags
    // PetscCall(DMSetMatType(da, MATAIJCUSPARSE));
    // PetscCall(DMSetVecType(da, VECCUDA));

    Mat A; Vec x, b;
    PetscCall(DMCreateMatrix(da, &A));
    PetscCall(DMCreateGlobalVector(da, &x));
    PetscCall(VecDuplicate(x, &b));
    PetscCall(MatZeroEntries(A));
    PetscCall(VecSet(b, 0.0));

    // --- ASSEMBLY timer start ---
    PetscCall(PetscTime(&t_asm0));

    // Assemble interior stencil; collect boundary nodes (we'll enforce BCs properly later)
    PetscInt xs, ys, xm, ym;
    PetscCall(DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr));

    std::vector<MatStencil> bnd; bnd.reserve(2*(nx+ny)); // rough upper bound

    for (PetscInt j = ys; j < ys+ym; ++j) {
        for (PetscInt i = xs; i < xs+xm; ++i) {
            const PetscBool is_bnd =
                ((i==0) || (i==nx-1) || (j==0) || (j==ny-1)) ? PETSC_TRUE : PETSC_FALSE;

            if (is_bnd) {
                MatStencil row; row.i=i; row.j=j; row.k=0; row.c=0;
                bnd.push_back(row);            // don’t insert identity rows here
                continue;                      // BC handled after assembly
            }

            // interior 5-point FD with anisotropic Kx,Ky
            MatStencil row; row.i=i; row.j=j; row.k=0; row.c=0;

            const PetscScalar w = -Kx/(dx*dx);
            const PetscScalar e = -Kx/(dx*dx);
            const PetscScalar s = -Ky/(dy*dy);
            const PetscScalar n = -Ky/(dy*dy);
            const PetscScalar c = -(w+e+s+n);

            MatStencil cols[5];
            PetscScalar vals[5];

            cols[0].i=i-1; cols[0].j=j;   cols[0].k=0; cols[0].c=0; vals[0]=w; // west
            cols[1].i=i+1; cols[1].j=j;   cols[1].k=0; cols[1].c=0; vals[1]=e; // east
            cols[2].i=i;   cols[2].j=j-1; cols[2].k=0; cols[2].c=0; vals[2]=s; // south
            cols[3].i=i;   cols[3].j=j+1; cols[3].k=0; cols[3].c=0; vals[3]=n; // north
            cols[4] = row;                                         vals[4]=c; // center

            PetscCall(MatSetValuesStencil(A, 1, &row, 5, cols, vals, INSERT_VALUES));
        }
    }

    // Finish initial assembly (interior operator only)
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,   MAT_FINAL_ASSEMBLY));
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));

    // Build a vector with Dirichlet values at boundary nodes
    Vec xbc;
    PetscCall(VecDuplicate(x, &xbc));
    PetscCall(VecSet(xbc, 0.0));
    for (const auto& r : bnd) {
        const PetscReal xcoord = r.i * dx, ycoord = r.j * dy;
        const PetscScalar g = g_dirichlet(xcoord, ycoord, Lx, Ly);
        const PetscInt gi = r.j*nx + r.i; // dof=1, row-major global index
        PetscCall(VecSetValue(xbc, gi, g, INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(xbc));
    PetscCall(VecAssemblyEnd(xbc));

    // Properly enforce Dirichlet: zero rows AND columns, set diag=1, and adjust RHS with xbc
    if (!bnd.empty()) {
        PetscCall(MatZeroRowsColumnsStencil(A, (PetscInt)bnd.size(), bnd.data(),
                                            1.0, xbc, b));
    }

    // Reassemble after modification
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,   MAT_FINAL_ASSEMBLY));
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));

    // --- ASSEMBLY timer end ---
    PetscCall(PetscTime(&t_asm1));

    // Solve
    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetFromOptions(ksp));

    // --- SOLVE timer start ---
    PetscCall(PetscTime(&t_solve0));
    PetscCall(KSPSolve(ksp, b, x));
    // --- SOLVE timer end ---
    PetscCall(PetscTime(&t_solve1));

    // Center value
    const PetscInt ic = (nx-1)/2, jc = (ny-1)/2;
    PetscInt centerIdx = jc*nx + ic;
    PetscScalar xc;
    PetscCall(VecGetValues(x, 1, &centerIdx, &xc));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "h(center) = %g\n", (double)xc));

    // Report backend types
    const char *vtype=nullptr,*mtype=nullptr;
    PetscCall(VecGetType(x,&vtype));
    PetscCall(MatGetType(A,&mtype));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Vec type = %s, Mat type = %s\n", vtype, mtype));

    // Timers
    PetscCall(PetscTime(&t_total1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "Timing: assemble = %.6f s, solve = %.6f s, total = %.6f s\n",
                          (double)(t_asm1 - t_asm0), (double)(t_solve1 - t_solve0), (double)(t_total1 - t_total0)));

    // Cleanup
    PetscCall(VecDestroy(&xbc));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));
    PetscCall(DMDestroy(&da));
    PetscCall(PetscFinalize());
    return 0;
}
