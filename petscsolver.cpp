// File overview: petscsolver.cpp is part of the PETSc_test simulation/analysis workflow.
#include "petscsolver.h"
#include "petscmatrix.h"
#include "petscvector.h"

// Build a KSP (Krylov Subspace) solver object that will be configured and reused.
// PetscCallAbort is used consistently in this wrapper so failures are surfaced
// immediately with a clear PETSc error instead of being silently ignored.
PETScSolver::PETScSolver() { PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp_)); }

// Release PETSc resources deterministically.
// KSPDestroy() accepts a pointer-to-handle and nulls it out internally.
PETScSolver::~PETScSolver() { if (ksp_) PetscCallAbort(PETSC_COMM_WORLD, KSPDestroy(&ksp_)); }

void PETScSolver::setOptionsPrefix(const char* prefix) {
    // Prefix allows independent solver blocks to be configured via command-line
    // options (e.g., -darcy_ksp_type gmres, -darcy_pc_type ilu).
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetOptionsPrefix(ksp_, prefix));
}

void PETScSolver::setOperator(const PETScMatrix& A, const PETScMatrix* P) {
    // Set the linear operator used by KSPSolve.
    // If no preconditioning matrix is provided, reuse A for both purposes,
    // which is PETSc's standard "single matrix" configuration.
    Mat Am = A.raw();
    Mat Pm = P ? P->raw() : Am;
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp_, Am, Pm));
}

void PETScSolver::setFromOptions() {
    // Pull KSP/PC/tolerance choices from the PETSc options database.
    // Typical call order is:
    //   setOptionsPrefix(...)
    //   setOperator(...)
    //   setFromOptions()
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp_));
}

void PETScSolver::setUp() {
    // Force explicit setup/allocation now instead of on first solve.
    // Useful for catching setup problems earlier in a simulation workflow.
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetUp(ksp_));
}

void PETScSolver::setReusePreconditioner(bool reuse) {
    // PCSetReusePreconditioner controls whether PETSc is allowed to keep the
    // current preconditioner across solves. This can reduce setup cost when
    // matrix sparsity and conditioning change slowly between timesteps.
    PC pc = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp_, &pc));
    PetscCallAbort(PETSC_COMM_WORLD, PCSetReusePreconditioner(pc, reuse ? PETSC_TRUE : PETSC_FALSE));
}

void PETScSolver::setTolerances(PetscReal rtol, PetscReal atol, PetscReal dtol, PetscInt maxits) {
    // Configure convergence criteria:
    //   rtol   -> relative residual reduction target
    //   atol   -> absolute residual threshold
    //   dtol   -> divergence threshold
    //   maxits -> iteration cap
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetTolerances(ksp_, rtol, atol, dtol, maxits));
}

void PETScSolver::solve(const PETScVector& b, PETScVector& x) const {
    // Solve A x = b using the currently bound operators and options.
    // x must have a layout compatible with b/A (typically VecDuplicate(b)).
    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp_, b.raw(), x.raw()));
}

PETScVector PETScSolver::solveNew(const PETScVector& b) const {
    // Convenience path for single-call solve sites:
    // allocate an output vector with b's parallel layout, solve, and transfer
    // ownership into the PETScVector RAII wrapper.
    Vec xraw = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b.raw(), &xraw));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp_, b.raw(), xraw));
    return PETScVector::adopt(xraw);
}
