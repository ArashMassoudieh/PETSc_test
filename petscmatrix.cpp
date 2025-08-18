#include "petscmatrix.h"
#include "petscvector.h"   // for setDiagonal/multiply
#include <sstream>

PETScMatrix::PETScMatrix() = default;

PETScMatrix::PETScMatrix(PetscInt m, PetscInt n, PetscInt nzPerRow) {
    create(m, n, nzPerRow);
}

PETScMatrix::PETScMatrix(PETScMatrix&& other) noexcept : A_(other.A_) {
    other.A_ = nullptr;
}

PETScMatrix& PETScMatrix::operator=(PETScMatrix&& other) noexcept {
    if (this != &other) {
        destroy();
        A_ = other.A_;
        other.A_ = nullptr;
    }
    return *this;
}

PETScMatrix::~PETScMatrix() { destroy(); }

void PETScMatrix::create(PetscInt m, PetscInt n, PetscInt nzPerRow) {
    destroy();
    PetscCallAbort(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &A_));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetSizes(A_, PETSC_DECIDE, PETSC_DECIDE, m, n));

    // Force GPU AIJ (MPI vs SEQ decided by communicator size)
    int size=1; MPI_Comm_size(PETSC_COMM_WORLD, &size);
#if defined(PETSC_HAVE_CUDA)
    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(A_, size>1 ? MATMPIAIJCUSPARSE : MATSEQAIJCUSPARSE));
#elif defined(PETSC_HAVE_HIP)
    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(A_, size>1 ? MATMPIAIJHIPSPARSE : MATSEQAIJHIPSPARSE));
#elif defined(PETSC_HAVE_KOKKOS_KERNELS)
    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(A_, size>1 ? MATMPIAIJKOKKOS : MATSEQAIJKOKKOS));
#else
    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(A_, size>1 ? MATMPIAIJ : MATSEQAIJ));
#endif

    // Let options still override if you pass -mat_type later
    PetscCallAbort(PETSC_COMM_WORLD, MatSetFromOptions(A_));

    // Preallocate AIJ
    PetscBool isSeqAIJ = PETSC_FALSE, isMPIAIJ = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscObjectTypeCompare((PetscObject)A_, MATSEQAIJ, &isSeqAIJ));
    PetscCallAbort(PETSC_COMM_WORLD, PetscObjectTypeCompare((PetscObject)A_, MATMPIAIJ, &isMPIAIJ));
    if (isSeqAIJ || isMPIAIJ) {
        PetscCallAbort(PETSC_COMM_WORLD, MatSeqAIJSetPreallocation(A_, nzPerRow, nullptr));
        PetscCallAbort(PETSC_COMM_WORLD, MatMPIAIJSetPreallocation(A_, nzPerRow, nullptr, nzPerRow, nullptr));
    }
    PetscCallAbort(PETSC_COMM_WORLD, MatSetUp(A_));
}

void PETScMatrix::setValue(PetscInt i, PetscInt j, PetscScalar a, InsertMode mode) {
    PetscCallAbort(PETSC_COMM_WORLD, MatSetValue(A_, i, j, a, mode));
}

void PETScMatrix::setValues(PetscInt ni, const PetscInt rows[], PetscInt nj,
                            const PetscInt cols[], const PetscScalar *vals,
                            InsertMode mode) {
    PetscCallAbort(PETSC_COMM_WORLD, MatSetValues(A_, ni, rows, nj, cols, vals, mode));
}

void PETScMatrix::assemble() {
    PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY));
}

void PETScMatrix::getGlobalSize(PetscInt& m, PetscInt& n) const {
    PetscCallAbort(PETSC_COMM_WORLD, MatGetSize(A_, &m, &n));
}

void PETScMatrix::ownershipRange(PetscInt& rStart, PetscInt& rEnd) const {
    PetscCallAbort(PETSC_COMM_WORLD, MatGetOwnershipRange(A_, &rStart, &rEnd));
}

PetscInt PETScMatrix::nnz() const {
    MatInfo info;
    PetscCallAbort(PETSC_COMM_WORLD, MatGetInfo(A_, MAT_GLOBAL_SUM, &info));
    // info.nz_used is a double; round to PetscInt for a human-friendly nnz.
    return static_cast<PetscInt>(info.nz_used + 0.5);
}

void PETScMatrix::zeroEntries() {
    PetscCallAbort(PETSC_COMM_WORLD, MatZeroEntries(A_));
}

void PETScMatrix::scale(PetscScalar a) {
    PetscCallAbort(PETSC_COMM_WORLD, MatScale(A_, a));
}

void PETScMatrix::shift(PetscScalar a) {
    PetscCallAbort(PETSC_COMM_WORLD, MatShift(A_, a));
}

void PETScMatrix::setDiagonal(const PETScVector& d) {
    PetscCallAbort(PETSC_COMM_WORLD, MatDiagonalSet(A_, d.raw(), INSERT_VALUES));
}

void PETScMatrix::zeroRowsColumns(const std::vector<PetscInt>& rows, PetscScalar diag) {
    if (rows.empty()) return;
    // const_cast is safe: PETSc API takes non-const pointer but won't modify values.
    auto* rptr = const_cast<PetscInt*>(rows.data());
    PetscCallAbort(PETSC_COMM_WORLD, MatZeroRowsColumns(A_, static_cast<PetscInt>(rows.size()), rptr, diag, nullptr, nullptr));
}

void PETScMatrix::multiply(const PETScVector& x, PETScVector& y) const {
    PetscCallAbort(PETSC_COMM_WORLD, MatMult(A_, x.raw(), y.raw()));
}

PETScMatrix PETScMatrix::duplicate(bool copyValues) const {
    Mat B = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, MatDuplicate(A_, copyValues ? MAT_COPY_VALUES : MAT_DO_NOT_COPY_VALUES, &B));
    PETScMatrix out;
    out.A_ = B; // take ownership
    return out;
}

void PETScMatrix::CopyFrom(const PETScMatrix& other, bool samePattern) {
    PetscCallAbort(PETSC_COMM_WORLD, MatCopy(other.A_, A_,
                                             samePattern ? SAME_NONZERO_PATTERN : DIFFERENT_NONZERO_PATTERN));
}

std::string PETScMatrix::toStringSummary() const {
    PetscInt m=0, n=0;
    getGlobalSize(m, n);
    const char* t = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, MatGetType(A_, &t));
    std::ostringstream oss;
    oss << "Mat(type=" << (t ? t : "unknown") << ", " << m << "x" << n
        << ", nnz=" << nnz() << ")";
    return oss.str();
}

void PETScMatrix::destroy() {
    if (A_) {
        PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&A_));
        A_ = nullptr;
    }
}

bool PETScMatrix::solve(const PETScVector& b, PETScVector& x, const char* optionsPrefix) const {
    // Ensure x exists with compatible layout
    if (!x.raw()) {
        Vec xr = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b.raw(), &xr));
        x = PETScVector::adopt(xr);
    }
    // Zero initial guess (important on some GPU backends)
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(x.raw(), 0.0));

    KSP ksp = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp));
    if (optionsPrefix && *optionsPrefix) {
        PetscCallAbort(PETSC_COMM_WORLD, KSPSetOptionsPrefix(ksp, optionsPrefix));
    }
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp, A_, A_));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp, b.raw(), x.raw()));

    KSPConvergedReason reason;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetConvergedReason(ksp, &reason));
    PetscCallAbort(PETSC_COMM_WORLD, KSPDestroy(&ksp));
    return reason >= 0;  // true if converged, false if diverged
}

PETScVector PETScMatrix::solveNew(const PETScVector& b, const char* optionsPrefix) const {
    PETScVector x;  // empty
    bool ok = solve(b, x, optionsPrefix);
    if (!ok) {
        // Either throw, or return x as-is and let caller check;
        // choose one. Throwing provides a clear failure mode:
        throw std::runtime_error("KSP did not converge");
    }
    return x; // NRVO/move
}


PETScVector PETScMatrix::operator/(const PETScVector& b) const {
    return solveNew(b, nullptr);
}
