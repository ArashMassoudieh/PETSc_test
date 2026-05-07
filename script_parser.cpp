// script_parser.cpp
#include "script_parser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

// -----------------------------------------------------------------------
// Internal helpers
// -----------------------------------------------------------------------
namespace {

// Trim leading/trailing whitespace
std::string trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last  = s.find_last_not_of (" \t\r\n");
    return s.substr(first, last - first + 1);
}

// Lower-case copy
std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return s;
}

// Strip inline comment (everything from # onward, unless inside quotes)
std::string stripComment(const std::string& line)
{
    bool in_quotes = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') in_quotes = !in_quotes;
        if (!in_quotes && line[i] == '#') return line.substr(0, i);
    }
    return line;
}

// Split on first occurrence of delim
std::pair<std::string,std::string> splitFirst(const std::string& s, char delim)
{
    const auto pos = s.find(delim);
    if (pos == std::string::npos)
        return {s, ""};
    return {trim(s.substr(0, pos)), trim(s.substr(pos + 1))};
}

// Split comma-separated list:  "a, b, c" -> ["a","b","c"]
std::vector<std::string> splitList(const std::string& s)
{
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        const auto t = trim(token);
        if (!t.empty()) result.push_back(t);
    }
    return result;
}

// Expand integer range "start:end" (inclusive, step 1) to string list
// e.g. "1:5" -> ["1","2","3","4","5"]
std::vector<std::string> expandRange(const std::string& s,
                                     const std::string& src, int lineno)
{
    const auto pos = s.find(':');
    if (pos == std::string::npos) {
        // single value
        return {trim(s)};
    }
    std::string lo_s = trim(s.substr(0, pos));
    std::string hi_s = trim(s.substr(pos + 1));
    int lo, hi;
    try { lo = std::stoi(lo_s); hi = std::stoi(hi_s); }
    catch (...) {
        throw std::runtime_error(src + ":" + std::to_string(lineno)
                                 + ": bad range '" + s + "' (must be int:int)");
    }
    if (lo > hi)
        throw std::runtime_error(src + ":" + std::to_string(lineno)
                                 + ": range start > end in '" + s + "'");
    std::vector<std::string> result;
    for (int v = lo; v <= hi; ++v)
        result.push_back(std::to_string(v));
    return result;
}

// Parse a sweep values spec: either [a, b, c] (discrete list)
// or start:end (integer range), or a single value
std::vector<std::string> parseSweepValues(const std::string& spec,
                                          const std::string& src, int lineno)
{
    const std::string s = trim(spec);
    if (s.empty())
        throw std::runtime_error(src + ":" + std::to_string(lineno)
                                 + ": empty sweep values");

    // Bracketed list: [a, b, c]
    if (s.front() == '[' && s.back() == ']') {
        const std::string inner = s.substr(1, s.size() - 2);
        return splitList(inner);
    }

    // Integer range: start:end
    if (s.find(':') != std::string::npos)
        return expandRange(s, src, lineno);

    // Single value
    return {s};
}

// -----------------------------------------------------------------------
// Recursive parser:  parse lines[pos..end) until we hit "end sweep"
// or run out of lines.  Returns the index past the last consumed line.
// -----------------------------------------------------------------------
std::size_t parseBlock(
    const std::vector<std::string>& lines,
    const std::vector<int>&         linenums,
    std::size_t                     pos,
    std::size_t                     end,
    const std::string&              src,
    std::vector<Statement>&         out,
    bool                            inside_sweep)
{
    while (pos < end) {
        const int      lineno  = linenums[pos];
        const std::string raw  = lines[pos];
        const std::string line = trim(stripComment(raw));
        ++pos;

        // ---- Blank / pure comment ----
        if (line.empty()) {
            Statement s;
            s.kind   = StmtKind::Blank;
            s.lineno = lineno;
            out.push_back(s);
            continue;
        }

        const std::string low = lower(line);

        // ---- end sweep ----
        if (low == "end sweep" || low == "end_sweep" || low == "endsweep") {
            if (!inside_sweep)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'end sweep' without matching 'sweep'");
            return pos;   // caller consumes this token
        }

        // ---- set key = value ----
        if (low.substr(0, 4) == "set ") {
            const std::string rest = trim(line.substr(4));
            auto [key, value] = splitFirst(rest, '=');
            if (key.empty())
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'set' missing key");
            Statement s;
            s.kind   = StmtKind::Set;
            s.lineno = lineno;
            s.key    = key;
            s.value  = value;   // kept as raw string; $vars expanded at runtime
            out.push_back(s);
            continue;
        }

        // ---- run stage ----
        if (low.substr(0, 4) == "run ") {
            Statement s;
            s.kind   = StmtKind::Run;
            s.lineno = lineno;
            s.key    = trim(lower(line.substr(4)));  // e.g. "fine", "upscale", "copula_transport"
            if (s.key.empty())
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'run' missing stage name");
            out.push_back(s);
            continue;
        }

        // ---- write kind = filename ----
        if (low.substr(0, 6) == "write ") {
            const std::string rest = trim(line.substr(6));
            auto [kind, filename]  = splitFirst(rest, '=');
            if (kind.empty() || filename.empty())
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'write' syntax: write kind = filename");
            Statement s;
            s.kind   = StmtKind::Write;
            s.lineno = lineno;
            s.key    = trim(lower(kind));
            s.value  = trim(filename);
            out.push_back(s);
            continue;
        }

        // ---- score btc_compare against path ----
        if (low.substr(0, 6) == "score ") {
            // "score btc_compare against /some/path"
            const std::string rest = trim(line.substr(6));
            const auto apos = lower(rest).find(" against ");
            if (apos == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'score' syntax: score <kind> against <path>");
            Statement s;
            s.kind   = StmtKind::Score;
            s.lineno = lineno;
            s.key    = trim(rest.substr(apos + 9));   // reference path
            s.value  = trim(lower(rest.substr(0, apos)));  // kind e.g. "btc_compare"
            out.push_back(s);
            continue;
        }

        // ---- plot kind = filename ----
        if (low.substr(0, 5) == "plot ") {
            const std::string rest = trim(line.substr(5));
            auto [kind, filename]  = splitFirst(rest, '=');
            if (kind.empty() || filename.empty())
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'plot' syntax: plot kind = filename");
            Statement s;
            s.kind   = StmtKind::Plot;
            s.lineno = lineno;
            s.key    = trim(lower(kind));
            s.value  = trim(filename);
            out.push_back(s);
            continue;
        }

        // ---- sweep var in values ----
        // Syntax:  sweep VAR in [a, b, c]
        //      or  sweep VAR in start:end
        if (low.substr(0, 6) == "sweep ") {
            const std::string rest = trim(line.substr(6));
            // Find " in "
            const auto in_pos = lower(rest).find(" in ");
            if (in_pos == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                                         + ": 'sweep' syntax: sweep VAR in VALUES");
            const std::string var    = trim(rest.substr(0, in_pos));
            const std::string values = trim(rest.substr(in_pos + 4));

            Statement s;
            s.kind        = StmtKind::SweepBegin;
            s.lineno      = lineno;
            s.sweep_var   = var;
            s.sweep_values = parseSweepValues(values, src, lineno);

            // Recursively parse body until "end sweep"
            pos = parseBlock(lines, linenums, pos, end, src, s.body, /*inside_sweep=*/true);
            out.push_back(std::move(s));
            continue;
        }

        // ---- Unknown ----
        throw std::runtime_error(src + ":" + std::to_string(lineno)
                                 + ": unrecognised statement: '" + line + "'");
    }

    if (inside_sweep)
        throw std::runtime_error(src + ": reached end of file inside a sweep block "
                                 "(missing 'end sweep')");

    return pos;
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------
ScriptProgram parseScriptString(const std::string& text, const std::string& source_name)
{
    // Split into lines, tracking line numbers
    std::vector<std::string> lines;
    std::vector<int>         linenums;
    std::istringstream ss(text);
    std::string line;
    int n = 0;
    while (std::getline(ss, line)) {
        ++n;
        lines.push_back(line);
        linenums.push_back(n);
    }

    ScriptProgram program;
    parseBlock(lines, linenums, 0, lines.size(), source_name, program, /*inside_sweep=*/false);
    return program;
}

ScriptProgram parseScriptFile(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("parseScriptFile: cannot open '" + path + "'");

    std::ostringstream buf;
    buf << f.rdbuf();
    return parseScriptString(buf.str(), path);
}
