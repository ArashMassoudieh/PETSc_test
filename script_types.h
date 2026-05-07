// script_types.h
// Token definitions and runtime context for the upscaling script engine.
#pragma once

#include <string>
#include <vector>
#include <map>
#include <variant>
#include <stdexcept>

// -----------------------------------------------------------------------
// Statement types (one per line kind in the script)
// -----------------------------------------------------------------------
enum class StmtKind {
    Set,              // set key = value
    Run,              // run fine | upscale | copula_transport
    Write,            // write kind = filename
    Score,            // score btc_compare against path
    Plot,             // plot kind = filename
    SweepBegin,       // sweep var in [a,b,c] or start:end
    SweepEnd,         // end sweep
    Comment,          // # ...
    Blank             // empty line
};

// -----------------------------------------------------------------------
// A single parsed statement
// -----------------------------------------------------------------------
struct Statement {
    StmtKind    kind   = StmtKind::Blank;
    int         lineno = 0;

    // Set:  key / value
    // Run:  key = stage name ("fine", "upscale", "copula_transport")
    // Write: key = output kind, value = filename
    // Score: key = reference path
    // Plot:  key = plot kind, value = filename
    std::string key;
    std::string value;

    // Sweep: variable name, list of values (all as strings; numerics resolved at runtime)
    std::string            sweep_var;
    std::vector<std::string> sweep_values;  // discrete list OR expanded integer range

    // Body of sweep block (child statements)
    std::vector<Statement> body;
};

// -----------------------------------------------------------------------
// A parsed program = flat list of top-level statements
// (sweep bodies are nested inside Statement::body)
// -----------------------------------------------------------------------
using ScriptProgram = std::vector<Statement>;

// -----------------------------------------------------------------------
// Runtime variable map
// -----------------------------------------------------------------------
class ScriptContext
{
public:
    // Set / get string variables
    void        set(const std::string& key, const std::string& value);
    std::string get(const std::string& key) const;
    bool        has(const std::string& key) const;

    // Convenience typed getters (throw on missing / bad parse)
    double      getDouble (const std::string& key) const;
    int         getInt    (const std::string& key) const;
    bool        getBool   (const std::string& key) const;

    // Expand $variable references in a string
    std::string expand(const std::string& s) const;

    // Dump all vars (for debugging)
    void dump() const;

    // Raw map access (for iterating)
    const std::map<std::string, std::string>& vars() const { return vars_; }

private:
    std::map<std::string, std::string> vars_;
};
