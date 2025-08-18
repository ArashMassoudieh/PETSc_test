#pragma once
/**
 * @file petsc_vector.h
 * @brief RAII C++ wrapper around a PETSc Vec.
 *
 * Design goals:
 *  - Keep a minimal 1:1 mapping to PETSc semantics.
 *  - Own and lifetime-manage the underlying Vec handle.
 *  - Disallow copying (unique ownership), allow moves.
 *  - Be usable in both serial and MPI-parallel PETSc runs.
 *
 * Typical usage:
 *   PETScVector x(N);       // create distributed vector of size N
 *   x.set(0.0);             // fill with zeros
 *   x.setValue(i, v);       // enqueue set(s)
 *   x.assemble();           // finalize assembly
 *   auto n2 = x.norm(NORM_2);
 *
 * Notes:
 *  - All PETSc calls are wrapped with PetscCallAbort(PETSC_COMM_WORLD, ...),
 *    which aborts on error across MPI ranks.
 *  - You can fetch the raw PETSc Vec with raw() for advanced operations.
 */

#include <petscksp.h>

class PETScVector {
public:
    /// Construct an empty (null) wrapper. Call create() before use.
    PETScVector();

    /// Construct and allocate a distributed PETSc Vec of global size n.
    explicit PETScVector(PetscInt n);

    // Non-copyable (unique ownership of the PETSc handle)
    PETScVector(const PETScVector&) = delete;
    PETScVector& operator=(const PETScVector&) = delete;

    // Moveable
    PETScVector(PETScVector&& other) noexcept;
    PETScVector& operator=(PETScVector&& other) noexcept;

    /// Destroy and free any owned PETSc Vec.
    ~PETScVector();

    /**
   * @brief Create/replace the underlying Vec with global size n.
   *        Any previously held Vec is destroyed first.
   * @param n Global vector length
   */
    void create(PetscInt n);

    /**
   * @brief Set all entries to a scalar.
   * @param a Value to assign to each entry
   */
    void set(PetscScalar a);

    /**
   * @brief Set one entry (insert or add).
   * @param i    Global index
   * @param a    Value
   * @param mode INSERT_VALUES (default) or ADD_VALUES
   *
   * After a sequence of setValue calls you must call assemble().
   */
    void setValue(PetscInt i, PetscScalar a, InsertMode mode = INSERT_VALUES);

    /**
   * @brief Finalize any enqueued modifications (must be called after setValue).
   *        Safe to call multiple times; subsequent calls are cheap no-ops.
   */
    void assemble();

    /**
   * @brief Compute a vector norm.
   * @param type Norm type (e.g., NORM_2, NORM_INFINITY, ...)
   * @return     The norm value
   */
    PetscReal norm(NormType type = NORM_2) const;

    /// Non-owning access to the underlying PETSc Vec.
    Vec raw() const { return v_; }

    /// Convert vector contents to a string "[x0,x1,x2,...]"
    std::string toString(int precision = 6) const;

    /**
   * @brief Deep copy creator: returns a new PETScVector with the same layout and values.
   * @return A new PETScVector owning a duplicated PETSc Vec.
   *
   * Implementation detail:
   *  - Uses VecDuplicate(raw(), &w) and VecCopy(raw(), w).
   */
    PETScVector duplicate() const;

    /**
   * @brief Deep copy into *this* from another vector with compatible layout.
   * @param other Source vector
   *
   * Preconditions:
   *  - *this must already be created (create()) and have the same distribution
   *    and global size as @p other. PETSc will error out if layouts mismatch.
   */
    void CopyFrom(const PETScVector& other);

    /** Wrap and take ownership of an existing PETSc Vec (created elsewhere). */
    static PETScVector adopt(Vec v);
private:
    /// Helper: destroy any owned Vec and reset handle to null.
    void destroy();

    Vec v_ = nullptr;  ///< Owned PETSc vector handle (or nullptr if empty)
};
