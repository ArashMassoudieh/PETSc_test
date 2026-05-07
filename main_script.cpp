// main.cpp
// Script-driven entry point for the upscaling simulation framework.
//
// Usage:
//   mpirun -np 4 ./PETSc_test --script=my_run.uscript
//
// All simulation parameters, run modes, output paths, and loops are
// specified in the .uscript file.  See example.uscript for reference.
//
// The old hard-coded main is preserved as main_legacy.cpp if you need
// to revert.

#include "petsc_init.h"
#include <iostream>
#include <string>
#include <mpi.h>

#include "script_executor.h"

int main(int argc, char** argv)
{
    PETScInit petsc(argc, argv);

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // ---- Find --script= argument ----
    std::string script_path;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--script=", 0) == 0) {
            script_path = a.substr(9);
            break;
        }
    }

    if (script_path.empty()) {
        if (rank == 0) {
            std::cerr << "Usage: PETSc_test --script=<file.uscript>\n\n"
                      << "No script file specified.  Please provide a .uscript file.\n"
                      << "See example.uscript for the script syntax reference.\n";
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        return 1;
    }

    // ---- Run the script ----
    const bool ok = runScriptFile(script_path, rank);

    MPI_Barrier(PETSC_COMM_WORLD);
    return ok ? 0 : 1;
}
