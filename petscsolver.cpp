#include "petscsolver.h"
#include "petscmatrix.h"
#include "petscvector.h"

PETScSolver::PETScSolver() { PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp_)); }
PETScSolver::~PETScSolver() { if (ksp_) PetscCallAbort(PETSC_COMM_WORLD, KSPDestroy(&ksp_)); }

void PETScSolver::setOptionsPrefix(const char* prefix) {
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetOptionsPrefix(ksp_, prefix));
}

void PETScSolver::setOperator(const PETScMatrix& A, const PETScMatrix* P) {
    Mat Am = A.raw();
    Mat Pm = P ? P->raw() : Am;
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp_, Am, Pm));
}

void PETScSolver::setFromOptions() {
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp_));
}

void PETScSolver::setUp() {
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetUp(ksp_));
}

void PETScSolver::setReusePreconditioner(bool reuse) {
    PC pc = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp_, &pc));
    PetscCallAbort(PETSC_COMM_WORLD, PCSetReusePreconditioner(pc, reuse ? PETSC_TRUE : PETSC_FALSE));
}

void PETScSolver::setTolerances(PetscReal rtol, PetscReal atol, PetscReal dtol, PetscInt maxits) {
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetTolerances(ksp_, rtol, atol, dtol, maxits));
}

void PETScSolver::solve(const PETScVector& b, PETScVector& x) const {
    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp_, b.raw(), x.raw()));
}

PETScVector PETScSolver::solveNew(const PETScVector& b) const {
    Vec xraw = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b.raw(), &xraw));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp_, b.raw(), xraw));
    return PETScVector::adopt(xraw);
}
