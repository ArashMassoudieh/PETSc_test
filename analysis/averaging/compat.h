#include <petscsys.h>
#ifndef PetscCallAbort
#define PetscCallAbort(comm, ...) do { PetscErrorCode review_ierr = (__VA_ARGS__); CHKERRABORT(comm, review_ierr); } while (0)
#endif
