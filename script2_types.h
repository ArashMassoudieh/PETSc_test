// script2_types.h
// Types for the direct Grid2D script engine (no sim_runner dependency).
#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <sstream>

// -----------------------------------------------------------------------
// Statement kinds
// -----------------------------------------------------------------------
enum class S2Kind {
    Blank,          // empty line or comment
    Set,            // set var = value
    GridCreate,     // grid NAME = nx:N, ny:M, Lx:L, Ly:L
    AccumCreate,    // accumulator NAME = bins:N
    MethodCall,     // NAME.method  arg=val, arg=val, ...
    RepeatBegin,    // repeat VAR = start:end
    BlockOpen,      // {
    BlockClose,     // }
};

// -----------------------------------------------------------------------
// A key=value argument on a method call line
// -----------------------------------------------------------------------
struct S2Arg {
    std::string key;
    std::string value;   // raw string; $vars expanded at runtime
};

// -----------------------------------------------------------------------
// A single parsed statement
// -----------------------------------------------------------------------
struct S2Stmt {
    S2Kind               kind    = S2Kind::Blank;
    int                  lineno  = 0;

    // Set
    std::string          var;       // variable name
    std::string          val;       // raw value (may contain $refs)

    // GridCreate / AccumCreate: object name + args
    std::string          obj;       // object name  (e.g. "g", "adv_copula")
    std::vector<S2Arg>   args;      // parsed key=value pairs

    // MethodCall: object name + method name + args
    std::string          method;    // method name  (e.g. "create_field")

    // RepeatBegin: loop variable + range
    std::string          loop_var;
    int                  loop_start = 1;
    int                  loop_end   = 1;

    // Body of repeat block (nested statements)
    std::vector<S2Stmt>  body;
};

using S2Program = std::vector<S2Stmt>;

// -----------------------------------------------------------------------
// Runtime variable map
// -----------------------------------------------------------------------
class S2Context
{
public:
    void        set(const std::string& k, const std::string& v) { vars_[k] = v; }
    std::string get(const std::string& k) const {
        auto it = vars_.find(k);
        if (it == vars_.end())
            throw std::runtime_error("S2Context: undefined variable '" + k + "'");
        return it->second;
    }
    bool has(const std::string& k) const { return vars_.count(k) > 0; }

    double getDouble(const std::string& k) const { return std::stod(get(k)); }
    int    getInt   (const std::string& k) const { return std::stoi(get(k)); }
    bool   getBool  (const std::string& k) const {
        std::string s = get(k);
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return (char)std::tolower(c); });
        return s == "true" || s == "yes" || s == "1";
    }

    // Expand $var and ${var} references
    std::string expand(const std::string& s) const {
        std::string out;
        out.reserve(s.size());
        std::size_t i = 0;
        while (i < s.size()) {
            if (s[i] != '$') { out += s[i++]; continue; }
            ++i;
            bool braced = (i < s.size() && s[i] == '{');
            if (braced) ++i;
            std::size_t j = i;
            while (j < s.size() && (std::isalnum((unsigned char)s[j]) || s[j] == '_')) ++j;
            std::string name = s.substr(i, j - i);
            if (braced && j < s.size() && s[j] == '}') ++j;
            if (has(name)) out += get(name);
            else { std::cerr << "[warn] expand: unknown $" << name << "\n"; out += '$'; out += name; }
            i = j;
        }
        return out;
    }

    // Lookup arg by key from a list, return expanded value or default
    std::string argVal(const std::vector<S2Arg>& args,
                       const std::string& key,
                       const std::string& def = "") const
    {
        for (auto& a : args)
            if (a.key == key) return expand(a.value);
        return def;
    }

    bool argHas(const std::vector<S2Arg>& args, const std::string& key) const {
        for (auto& a : args) if (a.key == key) return true;
        return false;
    }

    double argDouble(const std::vector<S2Arg>& args,
                     const std::string& key, double def = 0.0) const {
        const std::string v = argVal(args, key, "");
        return v.empty() ? def : std::stod(v);
    }

    int argInt(const std::vector<S2Arg>& args,
               const std::string& key, int def = 0) const {
        const std::string v = argVal(args, key, "");
        return v.empty() ? def : std::stoi(v);
    }

    void dump() const {
        for (auto& [k, v] : vars_) std::cout << "  " << k << " = " << v << "\n";
    }

private:
    std::map<std::string, std::string> vars_;
};

// -----------------------------------------------------------------------
// Self-contained copula accumulator
// Accumulates binned copula matrices across realizations and computes
// the ensemble mean.
//
// Each call to add() contributes one (n_bins x n_bins) matrix.
// The mean is the element-wise average, re-normalised so rows sum to 1/n_bins.
// -----------------------------------------------------------------------
class CopulaAccumulator
{
public:
    explicit CopulaAccumulator(int bins = 20) : bins_(bins) {}

    // Add one realization's normalized binned copula (row-major, size bins x bins,
    // values = theta * du so each row sums to du = 1/bins).
    void add(const std::vector<double>& mat)
    {
        if ((int)mat.size() != bins_ * bins_)
            throw std::runtime_error("CopulaAccumulator::add: size mismatch");
        if (sum_.empty()) sum_.assign(bins_ * bins_, 0.0);
        for (int k = 0; k < bins_ * bins_; ++k) sum_[k] += mat[k];
        ++count_;
    }

    // Mean matrix (rows sum to 1/bins)
    std::vector<double> mean() const {
        if (count_ == 0 || sum_.empty())
            throw std::runtime_error("CopulaAccumulator::mean: no data");
        std::vector<double> m(sum_.size());
        for (std::size_t k = 0; k < sum_.size(); ++k) m[k] = sum_[k] / count_;
        return m;
    }

    int bins()  const { return bins_;  }
    int count() const { return count_; }
    bool empty() const { return count_ == 0; }

    // Write mean as CSV (bins x bins, comma-separated, rows sum to 1/bins)
    void saveMean(const std::string& path) const;

private:
    int                  bins_  = 20;
    int                  count_ = 0;
    std::vector<double>  sum_;
};
