// petsc_init.h
#pragma once
#include <petscsys.h>

/**
 * @brief RAII wrapper for PETSc initialization/finalization.
 * Construct at the start of main(); it calls PetscInitialize.
 * When it goes out of scope, it calls PetscFinalize (if appropriate).
 */
class PETScInit {
public:
    PETScInit(int &argc, char **&argv) {
        PetscCallAbort(PETSC_COMM_WORLD, PetscInitialize(&argc, &argv, nullptr, nullptr));
    }

    // non-copyable
    PETScInit(const PETScInit&) = delete;
    PETScInit& operator=(const PETScInit&) = delete;

    ~PETScInit() noexcept {
        // Only finalize if PETSc is initialized and not yet finalized
        PetscBool inited = PETSC_FALSE, fin = PETSC_FALSE;
        (void)PetscInitialized(&inited);  // safe even if not initialized
        if (inited) {
            (void)PetscFinalized(&fin);
            if (!fin) {
                (void)PetscFinalize();        // avoid throwing from destructor
            }
        }
    }
};
