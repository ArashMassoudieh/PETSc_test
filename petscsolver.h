#pragma once
/**
 * @file petsc_solver.h
 * @brief Reusable PETSc KSP wrapper for high-throughput solves.
 *
 * Create once, set operators/options, then call solve() many times.
 * Works with PETScMatrix / PETScVector wrappers.
 */

#include <petscksp.h>
#include <string>

class PETScMatrix;
class PETScVector;

class PETScSolver {
public:
    PETScSolver();
    ~PETScSolver();

    PETScSolver(const PETScSolver&) = delete;
    PETScSolver& operator=(const PETScSolver&) = delete;
    PETScSolver(PETScSolver&&) = delete;
    PETScSolver& operator=(PETScSolver&&) = delete;

    /** Optional: set an options database prefix (read -<prefix>ksp_type, etc.). */
    void setOptionsPrefix(const char* prefix);     // e.g., "darcy_"

    /** Bind operator matrices (A for operator and (optionally) P for preconditioning). */
    void setOperator(const PETScMatrix& A, const PETScMatrix* P = nullptr);

    /** Read options (-ksp_type, -pc_type, tolerances, etc.). */
    void setFromOptions();

    /** Optional: explicitly allocate/set up the solver now. */
    void setUp();

    /** Reuse preconditioner between solves (helpful when A’s values change but pattern doesn’t). */
    void setReusePreconditioner(bool reuse);

    /** Set tolerances (relative, absolute, divergence, max its). */
    void setTolerances(PetscReal rtol, PetscReal atol = PETSC_DEFAULT,
                       PetscReal dtol = PETSC_DEFAULT, PetscInt maxits = PETSC_DEFAULT);

    /** Solve A x = b. x must be created with compatible layout (or duplicate b). */
    void solve(const PETScVector& b, PETScVector& x) const;

    /** Quick one-liner: returns a new vector with solution (duplicates b’s layout). */
    PETScVector solveNew(const PETScVector& b) const;

    /** Access raw KSP. */
    KSP raw() const { return ksp_; }

private:
    KSP ksp_ = nullptr;
};
