// script2_parser.h
#pragma once
#include "script2_types.h"
#include <string>

// Parse a .uscript2 file into an S2Program.
// Throws std::runtime_error with file:lineno message on syntax errors.
S2Program parseScript2File  (const std::string& path);
S2Program parseScript2String(const std::string& text,
                              const std::string& source = "<string>");
