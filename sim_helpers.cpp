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
#include <cmath>
#include <cctype>
#include <set>

#include <dirent.h>

// ============================================================
// Internal helpers (only in this .cpp)
// ============================================================

static inline bool starts_with(const std::string& s, const std::string& p)
{
    return s.rfind(p, 0) == 0;
}

static inline std::string trim_copy(std::string s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
        [](unsigned char ch){ return !std::isspace(ch); }));
    s.erase(std::find_if(s.rbegin(), s.rend(),
        [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
    return s;
}

// IMPORTANT: remove ALL whitespace inside header tokens.
// This fixes: "FineDeriv_r0001_ x=0.50" -> "FineDeriv_r0001_x=0.50"
static inline std::string remove_all_whitespace(std::string s)
{
    s.erase(std::remove_if(s.begin(), s.end(),
                           [](unsigned char ch){ return std::isspace(ch); }),
            s.end());
    return s;
}

// strip surrounding quotes if present
static inline std::string strip_quotes(std::string s)
{
    if (!s.empty() && s.front()=='"' && s.back()=='"' && s.size() >= 2)
        return s.substr(1, s.size()-2);
    return s;
}

static inline std::string normalize_header_token(std::string s)
{
    s = trim_copy(std::move(s));
    s = strip_quotes(std::move(s));
    s = trim_copy(std::move(s));
    s = remove_all_whitespace(std::move(s)); // <-- the big fix
    return s;
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
static inline std::string base_name_from_header(const std::string& h_in)
{
    // normalize first (removes embedded spaces!)
    const std::string h = normalize_header_token(h_in);

    // Fine BTC realizations: Fine_r0001_<base>
    if (starts_with(h, "Fine_r")) {
        auto pos = h.find('_');            // after Fine
        pos = h.find('_', pos + 1);        // after r0001
        if (pos != std::string::npos && pos + 1 < h.size()) return h.substr(pos + 1);
        return "";
    }

    // Fine derivative realizations: FineDeriv_r0001_<base>
    if (starts_with(h, "FineDeriv_r")) {
        auto pos = h.find('_');            // after FineDeriv
        pos = h.find('_', pos + 1);        // after r0001
        if (pos != std::string::npos && pos + 1 < h.size()) return h.substr(pos + 1);
        return "";
    }

    // Fine mean
    if (starts_with(h, "FineMean_"))
        return h.substr(std::string("FineMean_").size());

    // Fine derivative mean
    if (starts_with(h, "FineDerivMean_"))
        return h.substr(std::string("FineDerivMean_").size());

    // Upscaled BTC mean
    if (starts_with(h, "Upscaled_mean_"))
        return h.substr(std::string("Upscaled_mean_").size());

    // Upscaled derivative mean
    if (starts_with(h, "UpscaledDeriv_mean_"))
        return h.substr(std::string("UpscaledDeriv_mean_").size());

    // fallback
    if (starts_with(h, "Upscaled_")) {
        auto pos = h.find('_');
        if (pos != std::string::npos && pos + 1 < h.size()) return h.substr(pos + 1);
        return "";
    }

    if (starts_with(h, "UpscaledDeriv_")) {
        auto pos = h.find('_');
        if (pos != std::string::npos && pos + 1 < h.size()) return h.substr(pos + 1);
        return "";
    }

    return "";
}

static inline bool read_csv_header(const std::string& csv_path,
                                  std::vector<std::string>& headers)
{
    std::ifstream f(csv_path);
    if (!f) return false;

    std::string line;
    if (!std::getline(f, line)) return false;

    headers.clear();
    std::string cur;
    std::stringstream ss(line);

    while (std::getline(ss, cur, ',')) {
        headers.push_back(normalize_header_token(cur));
    }
    return true;
}

// ============================================================
// FS / path helpers
// ============================================================

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

// ============================================================
// naming helpers
// ============================================================

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

// ============================================================
// NaN-safe accumulation helpers
// ============================================================

void accumulate_sum_count(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& count
) {
    if (sum.empty())   sum.assign(x.size(), 0.0);
    if (count.empty()) count.assign(x.size(), 0);

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

// ============================================================
// delimiter-robust parsing
// ============================================================

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

// ============================================================
// CSV utilities
// ============================================================

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

    // normalize headers (removes embedded spaces)
    for (auto& s : h) s = normalize_header_token(s);

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

// ============================================================
// resampling
// ============================================================

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

// ============================================================
// FineMean utilities
// ============================================================

void accumulate_sum(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& n
) {
    if (sum.empty()) sum.assign(x.size(), 0.0);
    if (n.empty())   n.assign(x.size(), 0);

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
    std::vector<double> out(sum.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < sum.size(); ++i) {
        if (i < n.size() && n[i] > 0) out[i] = sum[i] / (double)n[i];
    }
    return out;
}

// ============================================================
// resume/upscale-only readers
// ============================================================

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

// ============================================================
// gnuplot runner
// ============================================================

int run_gnuplot_script(const std::string& gp_path)
{
    std::string cmd = "gnuplot \"" + gp_path + "\"";
    return std::system(cmd.c_str());
}

// ============================================================
// gnuplot by basename
// ============================================================

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
        std::cerr << "WARNING: no recognizable columns (Fine*/FineMean*/Upscaled*) in: " << csv_path << "\n";
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

    // Decide which bases are "locations" for this CSV
    bool has_x = false, has_series = false;
    for (const auto& b : bases) {
        if (starts_with(b, "x=")) has_x = true;
        if (starts_with(b, "series_")) has_series = true;
    }

    // Prefer x= otherwise series_
    std::string must_prefix = "";
    if (has_x) must_prefix = "x=";
    else if (has_series) must_prefix = "series_";

    for (const auto& base : bases) {

        if (skip_base_t && (base == "t" || base == "time" || base == "Time")) continue;
        if (!must_prefix.empty() && !starts_with(base, must_prefix)) continue;

        std::vector<int> fine_cols;
        int fineMean_col = -1;
        int upMean_col   = -1;

        for (int c0 = 1; c0 < (int)hdr.size(); ++c0) {
            const std::string& name = hdr[c0]; // already normalized by read_csv_header
            const int col = c0 + 1;

            if ((starts_with(name, "Fine_r") || starts_with(name, "FineDeriv_r")) &&
                base_name_from_header(name) == base) {
                fine_cols.push_back(col);
            }
            else if ((starts_with(name, "FineMean_") || starts_with(name, "FineDerivMean_")) &&
                     base_name_from_header(name) == base) {
                fineMean_col = col;
            }
            else if ((starts_with(name, "Upscaled_mean_") ||
                      starts_with(name, "UpscaledDeriv_mean_") ||
                      starts_with(name, "UpscaledDeriv_") ||
                      starts_with(name, "Upscaled_")) &&
                     base_name_from_header(name) == base) {
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

        bool draw_gray = !(fine_cols.size() == 1 && fineMean_col > 0);

        if (draw_gray) {
            for (int col : fine_cols) {
                if (!first) gp << ", \\\n";
                first = false;
                gp << "  csvfile using 1:" << col
                   << " with lines lw 1 lc rgb \"#b0b0b0\" notitle";
            }
        }

        if (fineMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << fineMean_col
               << " with lines lw 4 lc rgb \"black\" title \"FineMean\"";
        }

        if (upMean_col > 0) {
            if (!first) gp << ", \\\n";
            first = false;
            gp << "  csvfile using 1:" << upMean_col
               << " with lines lw 4 lc rgb \"red\" title \"Upscaled mean\"";
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
