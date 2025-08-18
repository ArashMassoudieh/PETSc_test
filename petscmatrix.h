#pragma once
/**
 * @file petsc_matrix.h
 * @brief RAII C++ wrapper around a PETSc sparse Mat with convenience helpers.
 *
 * Design:
 *  - Unique ownership of underlying Mat (no copy; move-enabled).
 *  - Minimal sugar over PETSc semantics; raw() exposed for advanced use.
 *  - AIJ preallocation hint via nzPerRow (skipped if not AIJ).
 *
 * New utilities:
 *  - getGlobalSize(), ownershipRange(), nnz()
 *  - zeroEntries(), scale(), shift(), setDiagonal()
 *  - zeroRowsColumns()  // useful for Dirichlet via row/col elimination
 *  - multiply()         // y = A*x
 *  - duplicate(copyValues) // deep duplicate new matrix
 *  - CopyFrom(other, samePattern) // deep copy into *this*
 *  - toStringSummary()   // quick "type/m×n/nnz" string for logging
 */

#include <petscksp.h>
#include <string>
#include <vector>

class PETScVector;

class PETScMatrix {
public:
    /// Construct an empty (null) wrapper. Call create() before use.
    PETScMatrix();

    /**
   * @brief Construct and allocate a distributed Mat of size m x n.
   * @param m        Global rows
   * @param n        Global cols
   * @param nzPerRow Estimated nonzeros per row (AIJ preallocation hint)
   */
    PETScMatrix(PetscInt m, PetscInt n, PetscInt nzPerRow = 5);

    // Non-copyable (unique ownership)
    PETScMatrix(const PETScMatrix&) = delete;
    PETScMatrix& operator=(const PETScMatrix&) = delete;

    // Movable
    PETScMatrix(PETScMatrix&& other) noexcept;
    PETScMatrix& operator=(PETScMatrix&& other) noexcept;

    /// Destroy and free any owned PETSc Mat.
    ~PETScMatrix();

    // ---- Core lifetime / assembly ----

    /**
   * @brief Create/replace the underlying Mat with global size m x n.
   *        Any previously held Mat is destroyed first.
   */
    void create(PetscInt m, PetscInt n, PetscInt nzPerRow = 5);

    /** Set a single entry (INSERT_VALUES default, or ADD_VALUES). */
    void setValue(PetscInt i, PetscInt j, PetscScalar a,
                  InsertMode mode = INSERT_VALUES);

    /** Set a dense block of entries (row-major vals). */
    void setValues(PetscInt ni, const PetscInt rows[], PetscInt nj,
                   const PetscInt cols[], const PetscScalar *vals,
                   InsertMode mode = INSERT_VALUES);

    /** Finalize enqueued modifications. Safe to call multiple times. */
    void assemble();

    // ---- Info / convenience ----

    /** Get global size (m, n). */
    void getGlobalSize(PetscInt& m, PetscInt& n) const;

    /**
   * @brief Get owned row range [rStart, rEnd) for this rank.
   *        Use when assembling in parallel to only touch owned rows.
   */
    void ownershipRange(PetscInt& rStart, PetscInt& rEnd) const;

    /**
   * @brief Get approximate number of used nonzeros (sum over rows).
   * @return nnz (rounded from MatInfo.nz_used)
   */
    PetscInt nnz() const;

    /** Set all numeric entries to zero (structure preserved). */
    void zeroEntries();

    /** Scale all numeric entries: A := a * A. */
    void scale(PetscScalar a);

    /** Shift diagonal: A := A + a * I. */
    void shift(PetscScalar a);

    /** Set diagonal from a vector (length must match global size). */
    void setDiagonal(const class PETScVector& d);

    /**
   * @brief Zero selected rows and columns, placing `diag` on the diagonal of zeroed rows.
   * @param rows  Global row indices to eliminate
   * @param diag  Diagonal value to set on each eliminated row (e.g., 1.0 for Dirichlet)
   *
   * This is often used to enforce Dirichlet BCs: zero rows/cols, set b[i]=h_bc, and put diag=1.
   */
    void zeroRowsColumns(const std::vector<PetscInt>& rows, PetscScalar diag);

    /** y := A * x (sizes/layout must be compatible). */
    void multiply(const class PETScVector& x, class PETScVector& y) const;

    /**
   * @brief Deep duplicate into a new matrix with same layout, optionally values.
   * @param copyValues  if true, copies numeric values; otherwise only structure
   */
    PETScMatrix duplicate(bool copyValues = true) const;

    /**
   * @brief Deep copy from another matrix into *this*.
   * @param other       Source matrix
   * @param samePattern If true, assume same nonzero pattern (faster).
   *                    If false, allow different pattern (may reallocate).
   */
    void CopyFrom(const PETScMatrix& other, bool samePattern = true);

    /** Quick summary string: "Mat(type=..., m×n, nnz=...)" */
    std::string toStringSummary() const;

    /// Non-owning access to the underlying PETSc Mat.
    Mat raw() const { return A_; }

    /**
   * @brief Solve A x = b using PETSc KSP and return x.
   * @param b              Right-hand side vector
   * @param optionsPrefix  Optional options database prefix (e.g., "mypfx_").
   *                       If provided, the solver reads -mypfx_ksp_type, -mypfx_pc_type, etc.
   *                       If null, reads the global -ksp_-pc_* options.
            */

    // Reuse: solve into existing x (allocates if x is empty)
    bool solve(const PETScVector& b, PETScVector& x, const char* optionsPrefix = nullptr) const;

    // One-shot: allocate new x, solve, return it
    PETScVector solveNew(const PETScVector& b, const char* optionsPrefix = nullptr) const;

    /// Sugar: x = A / b; equivalent to solve(b).
    PETScVector operator/(const PETScVector& b) const;



private:
    /// Helper: destroy any owned Mat and reset handle to null.
    void destroy();

    Mat A_ = nullptr;  ///< Owned PETSc matrix handle (or nullptr if empty)
};
