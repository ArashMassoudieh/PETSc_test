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

#include <set>        // <-- REQUIRED for std::set

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

/*
// python plot
bool write_btc_compare_plot_py(const std::string& py_path,
                              const std::string& csv_path,
                              const std::string& fig_prefix,
                              const std::string& y_label)
{
    std::ofstream f(py_path);
    if (!f) {
        std::cerr << "ERROR: cannot write python plot script: " << py_path << "\n";
        return false;
    }

    // Notes:
    // - your CSV layout is: first column "t", then many columns named like
    //   Fine_r0001_series_0, ..., FineMean_series_0, Upscaled_mean_series_0, etc.
    // - script groups by suffix after the last underscore: series_0, series_1, series_2, ...
    // - for each group it draws: all Fine_r#### in gray, FineMean in black, Upscaled_mean in red dashed
    f <<
R"PY(
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

csv_path   = r")PY" << csv_path << R"PY("
fig_prefix = r")PY" << fig_prefix << R"PY("
y_label    = r")PY" << y_label << R"PY("

df = pd.read_csv(csv_path)
if "t" not in df.columns:
    raise RuntimeError("CSV must have a 't' column as the first column.")

t = df["t"].to_numpy()

def group_key(col):
    # everything after last '_' (e.g., 'series_0')
    if col == "t":
        return None
    return col.split("_")[-2] + "_" + col.split("_")[-1] if col.count("_") >= 1 else col

# Build groups (series_0, series_1, series_2, ...)
groups = {}
for c in df.columns:
    if c == "t":
        continue
    k = group_key(c)
    groups.setdefault(k, []).append(c)

out_dir = os.path.dirname(csv_path)
if out_dir == "":
    out_dir = "."

for k, cols in sorted(groups.items()):
    fine_cols = [c for c in cols if c.startswith("Fine_r")]
    fine_mean = [c for c in cols if c.startswith("FineMean_")]
    up_mean   = [c for c in cols if c.startswith("Upscaled_mean") or c.startswith("UpscaledDeriv_mean") or c.startswith("Upscaled_mean_") or c.startswith("UpscaledDeriv_mean_")]

    # If your Upscaled columns are named like "Upscaled_mean_series_0" (as you requested earlier),
    # this will match "Upscaled_mean_..." already.
    # For safety, also accept "Upscaled_meanseries_0" etc. by using a fallback:
    if len(up_mean) == 0:
        up_mean = [c for c in cols if c.startswith("Upscaled")]

    plt.figure()
    # Realizations in gray
    for c in fine_cols:
        plt.plot(t, df[c].to_numpy(), alpha=0.25, linewidth=1.0)

    # FineMean (black)
    for c in fine_mean:
        plt.plot(t, df[c].to_numpy(), linewidth=2.5, label="FineMean")

    # Upscaled mean (red dashed)
    for c in up_mean:
        plt.plot(t, df[c].to_numpy(), linestyle="--", linewidth=2.5, label="Upscaled mean")

    plt.xlabel("Time")
    plt.ylabel(y_label)
    plt.title(k)
    plt.grid(True, alpha=0.25)

    # avoid duplicate legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    uniq = []
    seen = set()
    for h, lab in zip(handles, labels):
        if lab not in seen:
            uniq.append((h, lab))
            seen.add(lab)
    if uniq:
        plt.legend([u[0] for u in uniq], [u[1] for u in uniq], loc="best")

    out_png = os.path.join(out_dir, f"{fig_prefix}_{k}.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

print("Saved figures with prefix:", fig_prefix, "in", out_dir)
)PY";

    return true;
}

int run_python_script(const std::string& py_path)
{
    // best-effort: try python3 then python
    std::string cmd1 = "python3 \"" + py_path + "\"";
    int rc = std::system(cmd1.c_str());
    if (rc == 0) return 0;

    std::string cmd2 = "python \"" + py_path + "\"";
    rc = std::system(cmd2.c_str());
    return rc;
}

// gnuplot
// ------------------------------------------------------------
// Write gnuplot script for Fine vs Upscaled BTC comparison
// ------------------------------------------------------------
bool write_btc_compare_plot_gnuplot(const std::string& gp_path,
                                    const std::string& csv_path,
                                    const std::string& fig_prefix,
                                    const std::string& y_label)
{
    std::ofstream f(gp_path);
    if (!f) {
        std::cerr << "ERROR: cannot write gnuplot script: " << gp_path << "\n";
        return false;
    }

    f <<
R"GP(
# -------------------------------
# Fine vs Upscaled BTC comparison
# -------------------------------

set datafile separator ","
set key outside
set grid
set term pngcairo size 1200,800 enhanced font "Helvetica,14"

csvfile = ")GP" << csv_path << R"GP("
figpref = ")GP" << fig_prefix << R"GP("
ylabel  = ")GP" << y_label << R"GP("

# Read header
stats csvfile using 1 nooutput
ncols = STATS_columns

# First column is time
tcol = 1

# Identify unique series names (series_0, series_1, ...)
# We rely on column headers
series_names = ""

do for [c=2:ncols] {
    header = system(sprintf("awk -F, 'NR==1{print $%d}' %s", c, csvfile))
    # extract suffix after last underscore
    if (strlen(header) > 0) {
        n = words(header)
    }
}

# Manually loop series indices until no match found
# (robust for series_0, series_1, series_2 ...)
do for [k=0:20] {

    found = 0
    set output sprintf("%s_series_%d.png", figpref, k)
    set xlabel "Time"
    set ylabel ylabel
    set title sprintf("series_%d", k)

    plot \
)GP";

    // We generate plotting commands dynamically using awk
    // Gray lines: Fine_r####
    // Black: FineMean
    // Red dashed: Upscaled_mean
    f <<
R"GP(
    for [c=2:ncols] \
        ( system(sprintf("awk -F, 'NR==1 && $%d ~ /^Fine_r.*series_%d$/{print 1}' %s", c, k, csvfile)) \
          ? csvfile using tcol:c with lines lc rgb "#aaaaaa" lw 1 notitle : 1/0 ), \
    for [c=2:ncols] \
        ( system(sprintf("awk -F, 'NR==1 && $%d ~ /^FineMean.*series_%d$/{print 1}' %s", c, k, csvfile)) \
          ? csvfile using tcol:c with lines lc rgb "black" lw 3 title "Fine mean" : 1/0 ), \
    for [c=2:ncols] \
        ( system(sprintf("awk -F, 'NR==1 && $%d ~ /^Upscaled.*series_%d$/{print 1}' %s", c, k, csvfile)) \
          ? csvfile using tcol:c with lines lc rgb "red" lw 3 dt 2 title "Upscaled mean" : 1/0 )

    if (GPVAL_ERRNO == 0) {
        found = 1
    }

    if (!found) {
        unset output
        break
    }
}

print "Gnuplot BTC comparison figures written with prefix:", figpref
)GP";

    return true;
}

// ------------------------------------------------------------
// Run gnuplot
// ------------------------------------------------------------
int run_gnuplot_script(const std::string& gp_path)
{
    std::string cmd = "gnuplot \"" + gp_path + "\"";
    return std::system(cmd.c_str());
}

// simple gnuplot
static inline bool read_csv_header(const std::string& csv_path, std::vector<std::string>& headers)
{
    std::ifstream f(csv_path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;

    headers.clear();
    std::string cur;
    std::stringstream ss(line);
    while (std::getline(ss, cur, ',')) {
        // trim quotes/spaces
        cur.erase(cur.begin(), std::find_if(cur.begin(), cur.end(), [](unsigned char ch){ return !std::isspace(ch); }));
        cur.erase(std::find_if(cur.rbegin(), cur.rend(), [](unsigned char ch){ return !std::isspace(ch); }).base(), cur.end());
        if (!cur.empty() && cur.front()=='"' && cur.back()=='"' && cur.size()>=2) cur = cur.substr(1, cur.size()-2);
        headers.push_back(cur);
    }
    return true;
}

static inline bool ends_with(const std::string& s, const std::string& suf)
{
    return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}

int run_gnuplot_script(const std::string& gp_path)
{
    std::string cmd = "gnuplot \"" + gp_path + "\"";
    return std::system(cmd.c_str());
}

bool write_btc_compare_plot_gnuplot_simple(const std::string& gp_path,
                                           const std::string& csv_path,
                                           const std::string& fig_prefix,
                                           const std::string& y_label,
                                           int max_series)
{
    std::vector<std::string> hdr;
    if (!read_csv_header(csv_path, hdr) || hdr.size() < 2) {
        std::cerr << "ERROR: cannot read CSV header: " << csv_path << "\n";
        return false;
    }
    if (hdr[0] != "t") {
        std::cerr << "WARNING: expected first column 't' but got '" << hdr[0] << "'. Using col 1 as x.\n";
    }

    std::ofstream gp(gp_path);
    if (!gp) {
        std::cerr << "ERROR: cannot write gnuplot script: " << gp_path << "\n";
        return false;
    }

    gp << "set datafile separator \",\"\n";
    gp << "set term pngcairo size 1200,800 enhanced font \"Helvetica,14\"\n";
    gp << "set grid\n";
    gp << "set key outside\n";
    // Extract directory of csv_path
    std::string csv_dir = ".";
    {
        auto pos = csv_path.find_last_of("/\\");
        if (pos != std::string::npos) csv_dir = csv_path.substr(0, pos);
    }

    gp << "csvfile = \"" << csv_path << "\"\n";
    gp << "outdir  = \"" << csv_dir << "\"\n";
    gp << "ylabel  = \"" << y_label << "\"\n";
    gp << "figpref = \"" << fig_prefix << "\"\n\n";

    // For each series_k, gather columns:
    //  - Fine_r????_series_k  (many)
    //  - FineMean_series_k    (one)
    //  - Upscaled_mean_series_k (one)
    int made = 0;

    for (int k = 0; k <= max_series; ++k) {
        const std::string suf = "series_" + std::to_string(k);

        std::vector<int> fine_cols;
        int fineMean_col = -1;
        int upMean_col   = -1;

        for (int c = 1; c < (int)hdr.size(); ++c) { // c is 0-based index, but gnuplot uses 1-based
            const std::string& name = hdr[c];
            if (!ends_with(name, suf)) continue;

            if (name.rfind("Fine_r", 0) == 0) {
                fine_cols.push_back(c+1); // to 1-based
            } else if (name.rfind("FineMean_", 0) == 0) {
                fineMean_col = c+1;
            } else if (name.rfind("Upscaled_mean", 0) == 0 || name.rfind("UpscaledDeriv_mean", 0) == 0 || name.rfind("Upscaled", 0) == 0) {
                upMean_col = c+1;
            }
        }

        if (fine_cols.empty() && fineMean_col < 0 && upMean_col < 0) {
            // no such series => stop after we have made at least one
            if (made > 0) break;
            else continue;
        }

        gp << "set output sprintf(\"%s/%s_" << base_file << ".png\", outdir, figpref)\n";
        gp << "set xlabel \"Time\"\n";
        gp << "set ylabel ylabel\n";
        gp << "set title \"" << suf << "\"\n";

        gp << "plot \\\n";

        bool first = true;

        // fine realizations (gray, no title)
        for (int col : fine_cols) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << col << " with lines lw 1 lc rgb \"#aaaaaa\" notitle";
        }

        // fine mean (black)
        if (fineMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << fineMean_col << " with lines lw 3 lc rgb \"black\" title \"Fine mean\"";
        }

        // upscaled mean (red dashed)
        if (upMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << upMean_col << " with lines lw 3 dt 2 lc rgb \"red\" title \"Upscaled mean\"";
        }

        gp << "\n\n";
        made++;
    }

    gp << "unset output\n";
    gp << "print \"Wrote " << fig_prefix << "_series_*.png\"\n";

    if (made == 0) {
        std::cerr << "WARNING: no series plots generated; check column names in: " << csv_path << "\n";
    }

    return true;
}

// gnuplot by name
static inline bool read_csv_header(const std::string& csv_path, std::vector<std::string>& headers)
{
    std::ifstream f(csv_path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;

    headers.clear();
    std::string cur;
    std::stringstream ss(line);
    while (std::getline(ss, cur, ',')) {
        cur.erase(cur.begin(), std::find_if(cur.begin(), cur.end(),
            [](unsigned char ch){ return !std::isspace(ch); }));
        cur.erase(std::find_if(cur.rbegin(), cur.rend(),
            [](unsigned char ch){ return !std::isspace(ch); }).base(), cur.end());
        if (!cur.empty() && cur.front()=='"' && cur.back()=='"' && cur.size()>=2)
            cur = cur.substr(1, cur.size()-2);
        headers.push_back(cur);
    }
    return true;
}

static inline bool starts_with(const std::string& s, const std::string& p)
{
    return s.rfind(p, 0) == 0;
}

// Return safe filename token
static inline std::string sanitize_token(std::string s)
{
    for (char& c : s) {
        if (!(std::isalnum((unsigned char)c) || c=='_' || c=='-' )) c = '_';
    }
    // compress multiple underscores
    std::string out;
    out.reserve(s.size());
    bool prev_us = false;
    for (char c : s) {
        if (c=='_') {
            if (!prev_us) out.push_back(c);
            prev_us = true;
        } else {
            out.push_back(c);
            prev_us = false;
        }
    }
    return out;
}

// Extract base name after known prefixes
static inline std::string base_name_from_header(const std::string& h)
{
    // Fine realiz: Fine_r0001_<base>
    if (starts_with(h, "Fine_r")) {
        auto pos = h.find('_');                // after Fine
        pos = h.find('_', pos+1);              // after r0001
        if (pos != std::string::npos && pos+1 < h.size()) return h.substr(pos+1);
        return "";
    }
    if (starts_with(h, "FineMean_")) return h.substr(std::string("FineMean_").size());
    if (starts_with(h, "Upscaled_mean_")) return h.substr(std::string("Upscaled_mean_").size());
    if (starts_with(h, "UpscaledDeriv_mean_")) return h.substr(std::string("UpscaledDeriv_mean_").size());
    if (starts_with(h, "Upscaled_")) {
        // fallback
        auto pos = h.find('_');
        if (pos != std::string::npos && pos+1 < h.size()) return h.substr(pos+1);
        return "";
    }
    return "";
}

int run_gnuplot_script(const std::string& gp_path)
{
    std::string cmd = "gnuplot \"" + gp_path + "\"";
    return std::system(cmd.c_str());
}

bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label)
{
    std::vector<std::string> hdr;
    if (!read_csv_header(csv_path, hdr) || hdr.size() < 2) {
        std::cerr << "ERROR: cannot read CSV header: " << csv_path << "\n";
        return false;
    }

    // Collect unique base names
    std::set<std::string> bases;
    for (size_t i = 1; i < hdr.size(); ++i) {
        std::string b = base_name_from_header(hdr[i]);
        if (!b.empty()) bases.insert(b);
    }

    if (bases.empty()) {
        std::cerr << "WARNING: no recognizable columns (Fine_r*//*_____________FineMean/Upscaled*) in: " << csv_path << "\n";
        std::cerr << "First few headers:\n";
        for (size_t i = 0; i < std::min<size_t>(hdr.size(), 12); ++i)
            std::cerr << "  [" << i << "] " << hdr[i] << "\n";
        return false;
    }

    std::ofstream gp(gp_path);
    if (!gp) {
        std::cerr << "ERROR: cannot write gnuplot script: " << gp_path << "\n";
        return false;
    }

    gp << "set datafile separator \",\"\n";
    gp << "set term pngcairo size 1200,800 enhanced font \"Helvetica,14\"\n";
    gp << "set grid\n";
    gp << "set key outside\n";
    // Extract directory of csv_path
    std::string csv_dir = ".";
    {
        auto pos = csv_path.find_last_of("/\\");
        if (pos != std::string::npos) csv_dir = csv_path.substr(0, pos);
    }

    gp << "csvfile = \"" << csv_path << "\"\n";
    gp << "outdir  = \"" << csv_dir << "\"\n";
    gp << "ylabel  = \"" << y_label << "\"\n";
    gp << "figpref = \"" << fig_prefix << "\"\n\n";

    int made = 0;

    for (const auto& base : bases) {
        std::vector<int> fine_cols;
        int fineMean_col = -1;
        int upMean_col   = -1;

        for (int c0 = 1; c0 < (int)hdr.size(); ++c0) { // 0-based index
            const std::string& name = hdr[c0];
            const int col = c0 + 1; // gnuplot 1-based

            if (starts_with(name, "Fine_r") && base_name_from_header(name) == base) fine_cols.push_back(col);
            else if (starts_with(name, "FineMean_") && base_name_from_header(name) == base) fineMean_col = col;
            else if ((starts_with(name, "Upscaled_mean_") || starts_with(name, "UpscaledDeriv_mean_") || starts_with(name, "Upscaled_"))
                      && base_name_from_header(name) == base) upMean_col = col;
        }

        if (fine_cols.empty() && fineMean_col < 0 && upMean_col < 0) continue;

        std::string base_file = sanitize_token(base);

        gp << "set output sprintf(\"%s/%s_" << base_file << ".png\", outdir, figpref)\n";
        gp << "set xlabel \"Time\"\n";
        gp << "set ylabel ylabel\n";
        gp << "set title \"" << base << "\"\n";
        gp << "plot \\\n";

        bool first = true;
        for (int col : fine_cols) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << col << " with lines lw 1 lc rgb \"#aaaaaa\" notitle";
        }
        if (fineMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << fineMean_col << " with lines lw 3 lc rgb \"black\" title \"Fine mean\"";
        }
        if (upMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << upMean_col << " with lines lw 3 dt 2 lc rgb \"red\" title \"Upscaled mean\"";
        }
        gp << "\n\n";

        made++;
    }

    gp << "unset output\n";
    gp << "print \"Wrote " << fig_prefix << "_*.png\"\n";

    if (made == 0) {
        std::cerr << "WARNING: no plots generated (bases existed but no matching columns?) for: " << csv_path << "\n";
    }
    return true;
}
*/

// ============================================================
// GNUPlot helpers + plot writer (Fine vs Upscaled)
// ============================================================

// ---- forward declarations (gnuplot helpers) ----
static inline bool read_csv_header(const std::string& csv_path,
                                   std::vector<std::string>& headers);

static inline bool starts_with(const std::string& s,
                               const std::string& p);

static inline std::string sanitize_token(std::string s);

static inline std::string base_name_from_header(const std::string& h);

// ------------------------------------------------------------
// Helper implementations
// ------------------------------------------------------------
static inline bool read_csv_header(const std::string& csv_path, std::vector<std::string>& headers)
{
    std::ifstream f(csv_path);
    if (!f) return false;

    std::string line;
    if (!std::getline(f, line)) return false;

    headers.clear();
    std::string cur;
    std::stringstream ss(line);

    while (std::getline(ss, cur, ',')) {
        // trim
        cur.erase(cur.begin(), std::find_if(cur.begin(), cur.end(),
            [](unsigned char ch){ return !std::isspace(ch); }));
        cur.erase(std::find_if(cur.rbegin(), cur.rend(),
            [](unsigned char ch){ return !std::isspace(ch); }).base(), cur.end());

        // strip surrounding quotes
        if (!cur.empty() && cur.front()=='"' && cur.back()=='"' && cur.size()>=2)
            cur = cur.substr(1, cur.size()-2);

        headers.push_back(cur);
    }
    return true;
}

static inline bool starts_with(const std::string& s, const std::string& p)
{
    return s.rfind(p, 0) == 0;
}

// Return safe filename token
static inline std::string sanitize_token(std::string s)
{
    for (char& c : s) {
        if (!(std::isalnum((unsigned char)c) || c=='_' || c=='-')) c = '_';
    }

    // compress multiple underscores
    std::string out;
    out.reserve(s.size());
    bool prev_us = false;

    for (char c : s) {
        if (c == '_') {
            if (!prev_us) out.push_back(c);
            prev_us = true;
        } else {
            out.push_back(c);
            prev_us = false;
        }
    }
    return out;
}

// Extract base name after known prefixes
static inline std::string base_name_from_header(const std::string& h)
{
    // Fine realiz: Fine_r0001_<base>
    if (starts_with(h, "Fine_r")) {
        auto pos = h.find('_');            // after Fine
        pos = h.find('_', pos + 1);        // after r0001
        if (pos != std::string::npos && pos + 1 < h.size()) return h.substr(pos + 1);
        return "";
    }
    if (starts_with(h, "FineMean_")) return h.substr(std::string("FineMean_").size());
    if (starts_with(h, "Upscaled_mean_")) return h.substr(std::string("Upscaled_mean_").size());
    if (starts_with(h, "UpscaledDeriv_mean_")) return h.substr(std::string("UpscaledDeriv_mean_").size());

    // fallback
    if (starts_with(h, "Upscaled_")) {
        auto pos = h.find('_');
        if (pos != std::string::npos && pos + 1 < h.size()) return h.substr(pos + 1);
        return "";
    }
    return "";
}

// ------------------------------------------------------------
// Run gnuplot (linker needs this definition in some .cpp)
// ------------------------------------------------------------
int run_gnuplot_script(const std::string& gp_path)
{
    std::string cmd = "gnuplot \"" + gp_path + "\"";
    return std::system(cmd.c_str());
}

// ------------------------------------------------------------
// gnuplot by basename
//   - Realizations: gray thin
//   - FineMean: thick black
//   - Upscaled mean: thick red dashed
//   - Output PNGs saved next to CSV (run folder)
// ------------------------------------------------------------
bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label,
                                                bool skip_base_t)
{
    std::vector<std::string> hdr;
    if (!read_csv_header(csv_path, hdr) || hdr.size() < 2) {
        std::cerr << "ERROR: cannot read CSV header: " << csv_path << "\n";
        return false;
    }

    // Collect unique base names
    std::set<std::string> bases;
    for (size_t i = 1; i < hdr.size(); ++i) {
        std::string b = base_name_from_header(hdr[i]);
        if (!b.empty()) bases.insert(b);
    }

    if (bases.empty()) {
        std::cerr << "WARNING: no recognizable columns (Fine_r*/FineMean/Upscaled*) in: " << csv_path << "\n";
        std::cerr << "First few headers:\n";
        for (size_t i = 0; i < std::min<size_t>(hdr.size(), 12); ++i)
            std::cerr << "  [" << i << "] " << hdr[i] << "\n";
        return false;
    }

    std::ofstream gp(gp_path);
    if (!gp) {
        std::cerr << "ERROR: cannot write gnuplot script: " << gp_path << "\n";
        return false;
    }

    // Extract directory of csv_path (where we want PNGs)
    std::string csv_dir = ".";
    {
        auto pos = csv_path.find_last_of("/\\");
        if (pos != std::string::npos) csv_dir = csv_path.substr(0, pos);
    }

    gp << "set datafile separator \",\"\n";
    gp << "set term pngcairo size 1400,900 enhanced font \"Helvetica,16\"\n";
    gp << "set grid\n";
    gp << "set border linewidth 1.2\n";
    gp << "set tics out\n";
    gp << "set key outside\n";
    gp << "set key samplen 2 spacing 1.2\n";
    gp << "csvfile = \"" << csv_path << "\"\n";
    gp << "outdir  = \"" << csv_dir << "\"\n";
    gp << "ylabel  = \"" << y_label << "\"\n";
    gp << "figpref = \"" << fig_prefix << "\"\n\n";

    int made = 0;

    for (const auto& base : bases) {

        // Skip plotting time-like bases if requested
        if (skip_base_t && (base == "t" || base == "time" || base == "Time")) continue;

        std::vector<int> fine_cols;
        int fineMean_col = -1;
        int upMean_col   = -1;

        for (int c0 = 1; c0 < (int)hdr.size(); ++c0) { // 0-based index
            const std::string& name = hdr[c0];
            const int col = c0 + 1; // gnuplot 1-based

            if (starts_with(name, "Fine_r") && base_name_from_header(name) == base) {
                fine_cols.push_back(col);
            }
            else if (starts_with(name, "FineMean_") && base_name_from_header(name) == base) {
                fineMean_col = col;
            }
            else if ((starts_with(name, "Upscaled_mean_") || starts_with(name, "UpscaledDeriv_mean_") || starts_with(name, "Upscaled_"))
                      && base_name_from_header(name) == base) {
                upMean_col = col;
            }
        }

        if (fine_cols.empty() && fineMean_col < 0 && upMean_col < 0) continue;

        std::string base_file = sanitize_token(base);

        gp << "set output sprintf(\"%s/%s_" << base_file << ".png\", outdir, figpref)\n";
        gp << "set xlabel \"Time\"\n";
        gp << "set ylabel ylabel\n";
        gp << "set title \"" << base << "\"\n";
        gp << "plot \\\n";

        bool first = true;

        // If only 1 realization and FineMean exists, skip gray curve to avoid overlap/double thickness
        bool draw_gray = !(fine_cols.size() == 1 && fineMean_col > 0);

        // Fine realizations: gray thin
        if (draw_gray) {
            for (int col : fine_cols) {
                if (!first) gp << ", \\\n";
                first = false;
                gp << "  csvfile using 1:" << col
                   << " with lines lw 1 lc rgb \"#b0b0b0\" notitle";
            }
        }

        // FineMean: thick black
        if (fineMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << fineMean_col
               << " with lines lw 4 lc rgb \"black\" title \"FineMean\"";
        }

        // Upscaled mean: thick red dashed
        if (upMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << upMean_col
               << " with lines lw 4 dt 2 lc rgb \"red\" title \"Upscaled mean\"";
        }

        gp << "\n\n";
        made++;
    }

    gp << "unset output\n";
    gp << "print \"Wrote " << fig_prefix << "_*.png\"\n";

    if (made == 0) {
        std::cerr << "WARNING: no plots generated for: " << csv_path << "\n";
    }

    return true;
}
