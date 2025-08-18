#include "petscvector.h"
#include <vector>
#include <sstream>
#include <iomanip>

PETScVector::PETScVector() = default;

PETScVector::PETScVector(PetscInt n) { create(n); }

PETScVector::PETScVector(PETScVector&& other) noexcept : v_(other.v_) {
    other.v_ = nullptr;
}

PETScVector& PETScVector::operator=(PETScVector&& other) noexcept {
    if (this != &other) {
        destroy();
        v_ = other.v_;
        other.v_ = nullptr;
    }
    return *this;
}

PETScVector::~PETScVector() { destroy(); }

void PETScVector::create(PetscInt n) {
    destroy();
    PetscCallAbort(PETSC_COMM_WORLD, VecCreate(PETSC_COMM_WORLD, &v_));
    PetscCallAbort(PETSC_COMM_WORLD, VecSetSizes(v_, PETSC_DECIDE, n));
    PetscCallAbort(PETSC_COMM_WORLD, VecSetFromOptions(v_));
}

void PETScVector::set(PetscScalar a) {
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(v_, a));
}

void PETScVector::setValue(PetscInt i, PetscScalar a, InsertMode mode) {
    PetscCallAbort(PETSC_COMM_WORLD, VecSetValue(v_, i, a, mode));
}

void PETScVector::assemble() {
    PetscCallAbort(PETSC_COMM_WORLD, VecAssemblyBegin(v_));
    PetscCallAbort(PETSC_COMM_WORLD, VecAssemblyEnd(v_));
}

PetscReal PETScVector::norm(NormType type) const {
    PetscReal n{};
    PetscCallAbort(PETSC_COMM_WORLD, VecNorm(v_, type, &n));
    return n;
}

void PETScVector::destroy() {
    if (v_) {
        PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&v_));
        v_ = nullptr;
    }
}

std::string PETScVector::toString(int precision) const {
    PetscInt n = 0;
    PetscCallAbort(PETSC_COMM_WORLD, VecGetSize(v_, &n));

    // Build index list [0,1,...,n-1]
    std::vector<PetscInt> idx(n);
    for (PetscInt i = 0; i < n; ++i) idx[i] = i;

    std::vector<PetscScalar> vals(n);
    PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(v_, n, idx.data(), vals.data()));

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << "[";
    for (PetscInt i = 0; i < n; ++i) {
        oss << vals[i];
        if (i + 1 < n) oss << ",";
    }
    oss << "]";
    return oss.str();
}

PETScVector PETScVector::duplicate() const {
    Vec w = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(v_, &w));
    PetscCallAbort(PETSC_COMM_WORLD, VecCopy(v_, w));

    PETScVector out;      // empty wrapper
    out.v_ = w;           // take ownership of duplicated Vec
    return out;           // NRVO/move
}

void PETScVector::CopyFrom(const PETScVector& other) {
    // PETSc will check compatibility (sizes/layout). If mismatch, it errors out.
    PetscCallAbort(PETSC_COMM_WORLD, VecCopy(other.v_, v_));
}

PETScVector PETScVector::adopt(Vec v) {
    PETScVector out;
    out.v_ = v; // take ownership; destructor will VecDestroy(v)
    return out;
}
