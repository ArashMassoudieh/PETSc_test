// script_types.cpp
#include "script_types.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

// -----------------------------------------------------------------------
// ScriptContext
// -----------------------------------------------------------------------

void ScriptContext::set(const std::string& key, const std::string& value)
{
    vars_[key] = value;
}

std::string ScriptContext::get(const std::string& key) const
{
    auto it = vars_.find(key);
    if (it == vars_.end())
        throw std::runtime_error("ScriptContext: undefined variable '" + key + "'");
    return it->second;
}

bool ScriptContext::has(const std::string& key) const
{
    return vars_.count(key) > 0;
}

double ScriptContext::getDouble(const std::string& key) const
{
    const std::string s = get(key);
    try { return std::stod(s); }
    catch (...) {
        throw std::runtime_error("ScriptContext: cannot convert '" + key
                                 + "' = '" + s + "' to double");
    }
}

int ScriptContext::getInt(const std::string& key) const
{
    const std::string s = get(key);
    try { return std::stoi(s); }
    catch (...) {
        throw std::runtime_error("ScriptContext: cannot convert '" + key
                                 + "' = '" + s + "' to int");
    }
}

bool ScriptContext::getBool(const std::string& key) const
{
    std::string s = get(key);
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    if (s == "true"  || s == "yes" || s == "1") return true;
    if (s == "false" || s == "no"  || s == "0") return false;
    throw std::runtime_error("ScriptContext: cannot convert '" + key
                             + "' = '" + s + "' to bool");
}

// Expand $varname or ${varname} references in a string.
// Unknown variables are left as-is with a warning.
std::string ScriptContext::expand(const std::string& s) const
{
    std::string out;
    out.reserve(s.size());
    std::size_t i = 0;
    while (i < s.size()) {
        if (s[i] == '$') {
            ++i;
            bool braced = (i < s.size() && s[i] == '{');
            if (braced) ++i;

            // Collect variable name: alphanumeric + underscore
            std::size_t j = i;
            while (j < s.size() && (std::isalnum((unsigned char)s[j]) || s[j] == '_'))
                ++j;

            if (braced && j < s.size() && s[j] == '}') ++j;

            std::string varname = s.substr(i, j - i - (braced ? 1 : 0));
            if (braced) varname = s.substr(i, j - i - 1);
            else        varname = s.substr(i, j - i);

            if (has(varname))
                out += get(varname);
            else {
                std::cerr << "[warn] expand: unknown variable '$" << varname << "'\n";
                out += '$';
                out += varname;
            }
            i = j;
        } else {
            out += s[i++];
        }
    }
    return out;
}

void ScriptContext::dump() const
{
    std::cout << "--- ScriptContext variables ---\n";
    for (auto& [k, v] : vars_)
        std::cout << "  " << k << " = " << v << "\n";
    std::cout << "------------------------------\n";
}
