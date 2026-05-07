// script2_executor.h
#pragma once
#include "script2_types.h"
#include <string>

// Execute a parsed S2Program directly against Grid2D objects.
// rank: MPI rank (all ranks participate; rank 0 does I/O).
bool executeScript2(const S2Program& program, int rank);

// Convenience: parse + execute
bool runScript2File(const std::string& path, int rank);
