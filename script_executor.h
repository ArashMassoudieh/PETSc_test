// script_executor.h
#pragma once
#include "script_types.h"
#include "sim_runner.h"
#include <string>

// Run a parsed script program.
// rank: MPI rank (rank 0 does I/O, all ranks run solver)
// Returns false on fatal error.
bool executeScript(const ScriptProgram& program, int rank);

// Convenience: parse + execute a script file
bool runScriptFile(const std::string& path, int rank);
