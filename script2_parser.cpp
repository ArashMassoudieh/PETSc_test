// script2_parser.cpp
#include "script2_parser.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

// -----------------------------------------------------------------------
// Internal helpers
// -----------------------------------------------------------------------
namespace {

std::string trim(const std::string& s) {
    auto f = s.find_first_not_of(" \t\r\n");
    if (f == std::string::npos) return "";
    return s.substr(f, s.find_last_not_of(" \t\r\n") - f + 1);
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return s;
}

// Strip inline # comment (outside quotes)
std::string stripComment(const std::string& line) {
    bool inq = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') inq = !inq;
        if (!inq && line[i] == '#') return line.substr(0, i);
    }
    return line;
}

// Split "key=value, key=value, ..." into S2Arg list.
// Values may contain commas only if they don't conflict with arg separators,
// so we split on commas that are NOT inside brackets [] or quotes.
std::vector<S2Arg> parseArgs(const std::string& s,
                              const std::string& src, int lineno)
{
    std::vector<S2Arg> result;
    std::string cur;
    int depth = 0;
    bool inq = false;

    for (std::size_t i = 0; i <= s.size(); ++i) {
        char c = (i < s.size()) ? s[i] : ',';  // sentinel comma at end
        if (c == '"') inq = !inq;
        if (!inq) {
            if (c == '[' || c == '(') ++depth;
            else if (c == ']' || c == ')') --depth;
        }
        if (!inq && depth == 0 && c == ',') {
            const std::string tok = trim(cur);
            cur.clear();
            if (tok.empty()) continue;
            const auto eq = tok.find('=');
            if (eq == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": argument missing '=': '" + tok + "'");
            S2Arg a;
            a.key   = trim(tok.substr(0, eq));
            a.value = trim(tok.substr(eq + 1));
            result.push_back(a);
        } else {
            cur += c;
        }
    }
    return result;
}

// Parse "start:end" integer range → (start, end)
std::pair<int,int> parseIntRange(const std::string& s,
                                  const std::string& src, int lineno)
{
    const auto pos = s.find(':');
    if (pos == std::string::npos)
        throw std::runtime_error(src + ":" + std::to_string(lineno)
                                 + ": expected start:end range, got '" + s + "'");
    try {
        int a = std::stoi(trim(s.substr(0, pos)));
        int b = std::stoi(trim(s.substr(pos + 1)));
        return {a, b};
    } catch (...) {
        throw std::runtime_error(src + ":" + std::to_string(lineno)
                                 + ": non-integer in range '" + s + "'");
    }
}

// -----------------------------------------------------------------------
// Recursive block parser.
// Reads from lines[pos] until it finds a matching "}" or end of input.
// inside_block == true means we expect "}" to close the block.
// Returns index past the last consumed line.
// -----------------------------------------------------------------------
std::size_t parseBlock(
    const std::vector<std::string>& lines,
    const std::vector<int>&         linenums,
    std::size_t                     pos,
    const std::string&              src,
    std::vector<S2Stmt>&            out,
    bool                            inside_block)
{
    const std::size_t end = lines.size();

    while (pos < end) {
        const int         lineno = linenums[pos];
        const std::string raw    = lines[pos];
        const std::string line   = trim(stripComment(raw));
        ++pos;

        // ---- Blank / comment ----
        if (line.empty()) {
            S2Stmt s; s.kind = S2Kind::Blank; s.lineno = lineno;
            out.push_back(s);
            continue;
        }

        const std::string low = lower(line);

        // ---- Block close ----
        if (line == "}") {
            if (!inside_block)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": unexpected '}'");
            return pos;   // caller sees block is done
        }

        // ---- Block open (standalone "{") — should not appear at top level ----
        if (line == "{") {
            throw std::runtime_error(src + ":" + std::to_string(lineno)
                + ": '{' must follow a 'repeat' statement on the NEXT line");
        }

        // ---- set var = value ----
        if (low.substr(0, 4) == "set ") {
            const std::string rest = trim(line.substr(4));
            const auto eq = rest.find('=');
            if (eq == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": 'set' missing '='");
            S2Stmt s;
            s.kind   = S2Kind::Set;
            s.lineno = lineno;
            s.var    = trim(rest.substr(0, eq));
            s.val    = trim(rest.substr(eq + 1));
            out.push_back(s);
            continue;
        }

        // ---- grid NAME = nx:N, ny:M, Lx:L, Ly:L ----
        if (low.substr(0, 5) == "grid ") {
            const std::string rest = trim(line.substr(5));
            const auto eq = rest.find('=');
            if (eq == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": 'grid' syntax: grid NAME = nx:N, ny:M, Lx:L, Ly:L");
            S2Stmt s;
            s.kind   = S2Kind::GridCreate;
            s.lineno = lineno;
            s.obj    = trim(rest.substr(0, eq));
            // Re-parse args as key:value (using colon as separator in values)
            const std::string argtxt = trim(rest.substr(eq + 1));
            // Split on commas first, then on first colon per token
            std::stringstream ss(argtxt);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                tok = trim(tok);
                if (tok.empty()) continue;
                const auto co = tok.find(':');
                if (co == std::string::npos)
                    throw std::runtime_error(src + ":" + std::to_string(lineno)
                        + ": grid arg missing ':': '" + tok + "'");
                S2Arg a;
                a.key   = trim(tok.substr(0, co));
                a.value = trim(tok.substr(co + 1));
                s.args.push_back(a);
            }
            out.push_back(s);
            continue;
        }

        // ---- accumulator NAME = bins:N ----
        if (low.substr(0, 12) == "accumulator ") {
            const std::string rest = trim(line.substr(12));
            const auto eq = rest.find('=');
            if (eq == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": 'accumulator' syntax: accumulator NAME = bins:N");
            S2Stmt s;
            s.kind   = S2Kind::AccumCreate;
            s.lineno = lineno;
            s.obj    = trim(rest.substr(0, eq));
            const std::string argtxt = trim(rest.substr(eq + 1));
            const auto co = argtxt.find(':');
            if (co != std::string::npos) {
                S2Arg a; a.key = trim(argtxt.substr(0, co));
                a.value = trim(argtxt.substr(co + 1));
                s.args.push_back(a);
            }
            out.push_back(s);
            continue;
        }

        // ---- repeat VAR = start:end ----
        // Body must be on subsequent lines inside { }
        if (low.substr(0, 7) == "repeat ") {
            const std::string rest = trim(line.substr(7));
            const auto eq = rest.find('=');
            if (eq == std::string::npos)
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": 'repeat' syntax: repeat VAR = start:end");

            S2Stmt s;
            s.kind     = S2Kind::RepeatBegin;
            s.lineno   = lineno;
            s.loop_var = trim(rest.substr(0, eq));
            const auto [lo, hi] = parseIntRange(trim(rest.substr(eq + 1)), src, lineno);
            s.loop_start = lo;
            s.loop_end   = hi;

            // Next non-blank line must be "{"
            while (pos < end && trim(stripComment(lines[pos])).empty()) ++pos;
            if (pos >= end || trim(stripComment(lines[pos])) != "{")
                throw std::runtime_error(src + ":" + std::to_string(lineno)
                    + ": 'repeat' block must be followed by '{' on the next line");
            ++pos;  // consume "{"

            // Recursively parse body until "}"
            pos = parseBlock(lines, linenums, pos, src, s.body, /*inside_block=*/true);
            out.push_back(std::move(s));
            continue;
        }

        // ---- NAME.method  arg=val, arg=val, ... ----
        // Also handles NAME.method with no args (e.g. adv_copula.save_mean file=...)
        {
            const auto dot = line.find('.');
            if (dot != std::string::npos) {
                const std::string obj    = trim(line.substr(0, dot));
                const std::string rest   = trim(line.substr(dot + 1));

                // Split method name from args (first whitespace)
                const auto sp = rest.find_first_of(" \t");
                const std::string method = (sp == std::string::npos)
                                           ? rest : trim(rest.substr(0, sp));
                const std::string argtxt = (sp == std::string::npos)
                                           ? "" : trim(rest.substr(sp + 1));

                if (obj.empty() || method.empty())
                    throw std::runtime_error(src + ":" + std::to_string(lineno)
                        + ": bad method call: '" + line + "'");

                S2Stmt s;
                s.kind   = S2Kind::MethodCall;
                s.lineno = lineno;
                s.obj    = obj;
                s.method = lower(method);
                if (!argtxt.empty())
                    s.args = parseArgs(argtxt, src, lineno);
                out.push_back(s);
                continue;
            }
        }

        // ---- Unknown ----
        throw std::runtime_error(src + ":" + std::to_string(lineno)
            + ": unrecognised statement: '" + line + "'");
    }

    if (inside_block)
        throw std::runtime_error(src + ": reached end of file inside a block "
                                 "(missing '}')");
    return pos;
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------
S2Program parseScript2String(const std::string& text, const std::string& src)
{
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

    S2Program prog;
    parseBlock(lines, linenums, 0, src, prog, /*inside_block=*/false);
    return prog;
}

S2Program parseScript2File(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("parseScript2File: cannot open '" + path + "'");
    std::ostringstream buf;
    buf << f.rdbuf();
    return parseScript2String(buf.str(), path);
}
