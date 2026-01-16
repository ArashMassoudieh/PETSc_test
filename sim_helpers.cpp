#include "sim_helpers.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <cmath>        // NEW: std::isfinite

// POSIX directory scan (Linux)
#include <dirent.h>

// --------------------
// FS / path helpers
// --------------------
bool createDirectory(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        if (info.st_mode & S_IFDIR) return true;
        std::cerr << "Error: Path exists but is not a directory: " << path << "\n";
        return false;
    }
#ifdef _WIN32
    if (mkdir(path.c_str()) != 0) {
#else
    if (mkdir(path.c_str(), 0755) != 0) {
#endif
        if (errno != EEXIST) {
            std::cerr << "Error creating directory " << path << ": " << std::strerror(errno) << "\n";
            return false;
        }
    }
    std::cout << "Created output directory: " << path << "\n";
    return true;
}

bool fileExists(const std::string& path) {
    struct stat st;
    return (stat(path.c_str(), &st) == 0) && (st.st_mode & S_IFREG);
}

bool dirExists(const std::string& path) {
    struct stat st;
    return (stat(path.c_str(), &st) == 0) && (st.st_mode & S_IFDIR);
}

std::string joinPath(const std::string& dir, const std::string& filename) {
    if (dir.empty()) return filename;
    if (dir.back() == '/' || dir.back() == '\\') return dir + filename;
    return dir + "/" + filename;
}

// --------------------
// naming helpers
// --------------------
std::string makeTimestamp() {
    std::time_t now = std::time(nullptr);
    std::tm tm_now;
    localtime_r(&now, &tm_now);
    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

std::string makeRealLabel(int r1) {
    std::ostringstream oss;
    oss << "r" << std::setw(4) << std::setfill('0') << r1;
    return oss.str();
}

std::string makeFineFolder(int r1) {
    return "fine_" + makeRealLabel(r1);
}

double mean_of(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / (double)v.size();
}

void accumulate_sum_count(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& count
) {
    if (sum.empty())   sum.assign(x.size(), 0.0);
    if (count.empty()) count.assign(x.size(), 0);

    // If sizes don't match, fail fast (or you could throw/log).
    if (sum.size() != x.size() || count.size() != x.size()) return;

    for (size_t i = 0; i < x.size(); ++i) {
        const double v = x[i];
        if (is_finite_number(v)) {
            sum[i]   += v;
            count[i] += 1;
        }
    }
}

std::vector<double> finalize_mean_vec(
    const std::vector<double>& sum,
    const std::vector<int>& count
) {
    std::vector<double> m(sum.size(), std::numeric_limits<double>::quiet_NaN());
    if (sum.size() != count.size()) return m;

    for (size_t i = 0; i < sum.size(); ++i) {
        if (count[i] > 0) m[i] = sum[i] / static_cast<double>(count[i]);
    }
    return m;
}

// --------------------
// delimiter-robust parsing
// --------------------
char detect_delim(const std::string& header) {
    if (header.find(',')  != std::string::npos) return ',';
    if (header.find('\t') != std::string::npos) return '\t';
    if (header.find(';')  != std::string::npos) return ';';
    return ' ';
}

std::vector<std::string> split_line_delim(const std::string& line, char delim) {
    std::vector<std::string> out;

    if (delim == ' ') {
        std::istringstream iss(line);
        std::string tok;
        while (iss >> tok) out.push_back(tok);
        return out;
    }

    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') { in_quotes = !in_quotes; continue; }
        if (c == delim && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

bool try_parse_double(const std::string& s, double& v) {
    char* endp = nullptr;
    v = std::strtod(s.c_str(), &endp);
    return endp != s.c_str() && *endp == '\0';
}

// --------------------
// CSV utilities
// --------------------
bool read_time_series_table_csv(
    const std::string& path,
    std::vector<double>& times,
    std::vector<std::string>& colnames,
    std::vector<std::vector<double>>& cols
) {
    std::ifstream f(path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    char delim = detect_delim(header);
    auto h = split_line_delim(header, delim);
    if (h.size() < 2) return false;

    colnames.assign(h.begin() + 1, h.end());
    cols.assign(colnames.size(), std::vector<double>{});
    times.clear();

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto parts = split_line_delim(line, delim);
        if (parts.size() < 1) continue;

        double t;
        if (!try_parse_double(parts[0], t)) continue;
        times.push_back(t);

        for (size_t j = 0; j < colnames.size(); ++j) {
            double val = std::numeric_limits<double>::quiet_NaN();
            if (j + 1 < parts.size()) {
                double tmp;
                if (try_parse_double(parts[j + 1], tmp)) val = tmp;
            }
            cols[j].push_back(val);
        }
    }

    for (auto& c : cols) if (c.size() != times.size()) return false;
    return true;
}

bool write_comparison_csv(
    const std::string& out_path,
    const std::vector<double>& base_time,
    const std::vector<std::string>& out_colnames,
    const std::vector<std::vector<double>>& out_cols
) {
    if (out_colnames.size() != out_cols.size()) return false;
    for (auto& c : out_cols) if (c.size() != base_time.size()) return false;

    std::ofstream f(out_path);
    if (!f) return false;

    f << "t";
    for (auto& n : out_colnames) f << "," << n;
    f << "\n";

    for (size_t i = 0; i < base_time.size(); ++i) {
        f << std::setprecision(12) << base_time[i];
        for (size_t j = 0; j < out_cols.size(); ++j) {
            f << "," << std::setprecision(12) << out_cols[j][i];
        }
        f << "\n";
    }
    return true;
}

// --------------------
// resampling
// --------------------
double lerp(double x0, double y0, double x1, double y1, double x) {
    if (x1 == x0) return y0;
    const double a = (x - x0) / (x1 - x0);
    return y0 + a * (y1 - y0);
}

std::vector<double> resample_col_linear(
    const std::vector<double>& t_src,
    const std::vector<double>& y_src,
    const std::vector<double>& t_dst
) {
    std::vector<double> y_dst;
    y_dst.reserve(t_dst.size());

    if (t_src.empty() || y_src.size() != t_src.size()) {
        y_dst.assign(t_dst.size(), std::numeric_limits<double>::quiet_NaN());
        return y_dst;
    }

    size_t k = 0;
    for (double td : t_dst) {
        if (td <= t_src.front()) { y_dst.push_back(y_src.front()); continue; }
        if (td >= t_src.back())  { y_dst.push_back(y_src.back());  continue; }

        while (k + 1 < t_src.size() && t_src[k + 1] < td) ++k;
        y_dst.push_back(lerp(t_src[k], y_src[k], t_src[k+1], y_src[k+1], td));
    }
    return y_dst;
}

void resample_table_linear(
    const std::vector<double>& t_src,
    const std::vector<std::vector<double>>& cols_src,
    const std::vector<double>& t_dst,
    std::vector<std::vector<double>>& cols_dst
) {
    cols_dst.clear();
    cols_dst.reserve(cols_src.size());
    for (const auto& c : cols_src) cols_dst.push_back(resample_col_linear(t_src, c, t_dst));
}

// --------------------
// FineMean utilities (mean of what is on disk; no std)
// --------------------
void accumulate_sum(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& n
) {
    if (sum.empty()) sum.assign(x.size(), 0.0);
    if (n.empty())   n.assign(x.size(), 0);

    // If sizes mismatch, reset (caller should keep consistent sizes)
    if (sum.size() != x.size()) sum.assign(x.size(), 0.0);
    if (n.size()   != x.size()) n.assign(x.size(), 0);

    for (size_t i = 0; i < x.size(); ++i) {
        if (std::isfinite(x[i])) {
            sum[i] += x[i];
            n[i]   += 1;
        }
    }
}

std::vector<double> finalize_mean(
    const std::vector<double>& sum,
    const std::vector<int>& n
) {
    std::vector<double> out;
    out.resize(sum.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < sum.size(); ++i) {
        if (i < n.size() && n[i] > 0) out[i] = sum[i] / (double)n[i];
    }
    return out;
}

// --------------------
// resume/upscale-only readers
// --------------------
bool parse_keyval_file(const std::string& path, std::map<std::string, std::string>& kv) {
    std::ifstream f(path);
    if (!f) return false;

    kv.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        kv[line.substr(0, pos)] = line.substr(pos + 1);
    }
    return !kv.empty();
}

bool read_mean_params_txt(
    const std::string& path,
    double& lc_mean, double& lx_mean, double& ly_mean, double& dt_mean
) {
    std::map<std::string, std::string> kv;
    if (!parse_keyval_file(path, kv)) return false;

    auto getd = [&](const std::string& k, double& out)->bool{
        auto it = kv.find(k);
        if (it == kv.end()) return false;
        out = std::atof(it->second.c_str());
        return true;
    };

    bool ok = true;
    ok = ok && getd("lc_mean", lc_mean);
    ok = ok && getd("lambda_x_mean", lx_mean);
    ok = ok && getd("lambda_y_mean", ly_mean);
    ok = ok && getd("dt_mean", dt_mean);
    return ok;
}

bool read_mean_inverse_cdf_csv(const std::string& path,
                               std::vector<double>& u, std::vector<double>& v)
{
    std::ifstream f(path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    char delim = detect_delim(header);

    u.clear(); v.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto parts = split_line_delim(line, delim);
        if (parts.size() < 2) continue;

        double uu, vv;
        if (!try_parse_double(parts[0], uu)) continue;
        if (!try_parse_double(parts[1], vv)) continue;
        u.push_back(uu);
        v.push_back(vv);
    }
    return !u.empty() && u.size() == v.size();
}

bool read_xy_table(const std::string& path, std::vector<double>& x, std::vector<double>& y) {
    std::ifstream f(path);
    if (!f) return false;

    x.clear(); y.clear();
    std::string line;

    // detect delimiter from first non-empty line
    char delim = ' ';
    std::streampos pos0 = f.tellg();
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        delim = detect_delim(line);
        break;
    }
    f.clear();
    f.seekg(pos0);

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto parts = split_line_delim(line, delim);
        if (parts.size() < 2) continue;

        double xx, yy;
        if (!try_parse_double(parts[0], xx)) continue;
        if (!try_parse_double(parts[1], yy)) continue;

        x.push_back(xx);
        y.push_back(yy);
    }
    return !x.empty() && x.size() == y.size();
}

double interp1_linear(const std::vector<double>& x, const std::vector<double>& y, double xv) {
    if (x.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (x.size() == 1) return y[0];

    if (xv <= x.front()) return y.front();
    if (xv >= x.back())  return y.back();

    size_t k = 0;
    while (k + 1 < x.size() && x[k + 1] < xv) ++k;
    return lerp(x[k], y[k], x[k+1], y[k+1], xv);
}

bool read_fine_params_all_csv(
    const std::string& path,
    std::vector<double>& lc,
    std::vector<double>& lx,
    std::vector<double>& ly,
    std::vector<double>& dt
) {
    std::ifstream f(path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    char delim = detect_delim(header);

    lc.clear(); lx.clear(); ly.clear(); dt.clear();

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto p = split_line_delim(line, delim);
        if (p.size() < 5) continue;

        double v_lc, v_lx, v_ly, v_dt;
        if (!try_parse_double(p[1], v_lc)) continue;
        if (!try_parse_double(p[2], v_lx)) continue;
        if (!try_parse_double(p[3], v_ly)) continue;
        if (!try_parse_double(p[4], v_dt)) continue;

        lc.push_back(v_lc);
        lx.push_back(v_lx);
        ly.push_back(v_ly);
        dt.push_back(v_dt);
    }
    return !lc.empty();
}

std::vector<std::pair<int,std::string>> list_fine_folders(const std::string& run_dir) {
    std::vector<std::pair<int,std::string>> out;

    DIR* d = opendir(run_dir.c_str());
    if (!d) return out;

    struct dirent* ent;
    while ((ent = readdir(d)) != nullptr) {
        std::string name = ent->d_name;
        if (name == "." || name == "..") continue;
        if (name.rfind("fine_r", 0) != 0) continue;

        std::string s = name.substr(std::string("fine_r").size());
        if (s.empty()) continue;

        int rid = std::atoi(s.c_str());
        if (rid <= 0) continue;

        std::string full = joinPath(run_dir, name);
        if (dirExists(full)) out.push_back({rid, full});
    }
    closedir(d);

    std::sort(out.begin(), out.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });
    return out;
}

std::string fine_qx_cdf_path(const std::string& fine_dir, int r) {
    const std::string pfx = makeRealLabel(r) + "_";
    return joinPath(fine_dir, pfx + "qx_inverse_cdf.txt");
}

bool accumulate_inverse_cdf_on_grid(
    const std::string& fine_dir, int r,
    double du, int nU,
    std::vector<double>& invcdf_sum
) {
    std::string path = fine_qx_cdf_path(fine_dir, r);
    std::vector<double> u, v;
    if (!read_xy_table(path, u, v)) return false;

    for (int k = 0; k < nU; ++k) {
        double uk = k * du;
        invcdf_sum[k] += interp1_linear(u, v, uk);
    }
    return true;
}
