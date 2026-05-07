// script_parser.h
#pragma once
#include "script_types.h"
#include <string>

// Parse a .uscript file into a ScriptProgram.
// Throws std::runtime_error with file:lineno message on syntax errors.
ScriptProgram parseScriptFile(const std::string& path);

// Parse from a string (for testing)
ScriptProgram parseScriptString(const std::string& text, const std::string& source_name = "<string>");
