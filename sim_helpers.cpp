// sim_helpers.cpp
#include "sim_helpers.h"
#include "sim_runner.h" // for SimParams (used by run-tag helpers)

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <dirent.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <regex>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <limits>

// ============================================================
// FS / path helpers
// ============================================================
bool createDirectory(const std::string& path)
{
#if defined(_WIN32)
    int status = _mkdir(path.c_str());
#else
    int status = mkdir(path.c_str(), 0755);
#endif
    if (status == 0) return true;
    if (errno == EEXIST) return true;
    return false;
}

bool fileExists(const std::string& path)
{
    std::ifstream f(path.c_str());
    return (bool)f;
}

bool dirExists(const std::string& path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0) return false;
    return (info.st_mode & S_IFDIR) != 0;
}

std::string joinPath(const std::string& dir, const std::string& filename)
{
    if (dir.empty()) return filename;
    char last = dir.back();
    if (last == '/' || last == '\\') return dir + filename;
    return dir + "/" + filename;
}

// ============================================================
// naming helpers
// ============================================================
std::string makeTimestamp()
{
    std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream ss;
    ss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return ss.str();
}

std::string makeRealLabel(int r1)
{
    std::ostringstream ss;
    ss << "r" << std::setw(4) << std::setfill('0') << r1;
    return ss.str();
}

std::string makeFineFolder(int r1)
{
    return std::string("fine_") + makeRealLabel(r1);
}

// --------------------------------------------------
// Resume folder parsing helpers (robust)
// --------------------------------------------------

static inline std::string lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (unsigned char)std::tolower(c); });
    return s;
}

static bool contains_word_token_ci(const std::string& s_in, const std::string& word_in)
{
    const std::string s = lower_copy(s_in);
    const std::string word = lower_copy(word_in);

    const auto is_alpha = [](unsigned char c){ return std::isalpha(c) != 0; };

    size_t pos = 0;
    while ((pos = s.find(word, pos)) != std::string::npos) {
        const bool left_ok  = (pos == 0) || !is_alpha((unsigned char)s[pos - 1]);
        const bool right_ok = (pos + word.size() >= s.size()) ||
                              !is_alpha((unsigned char)s[pos + word.size()]);
        if (left_ok && right_ok) return true;
        pos += word.size();
    }
    return false;
}

bool parse_resume_aniso(const std::string& s, bool& is_aniso)
{
    if (contains_word_token_ci(s, "aniso")) { is_aniso = true;  return true; }
    if (contains_word_token_ci(s, "iso"))   { is_aniso = false; return true; }
    return false;
}

static inline bool is_key_boundary_left_ci(const std::string& s_lower, size_t pos)
{
    if (pos == 0) return true;
    unsigned char c = (unsigned char)s_lower[pos - 1];
    return !std::isalnum(c);
}

static inline bool is_number_start(unsigned char c)
{
    return (std::isdigit(c) || c == '.' || c == '+' || c == '-');
}

static bool extract_last_number_after_key_ci(
    const std::string& s_original,
    const std::string& key_lower,
    double& value_out
)
{
    value_out = 0.0;

    const std::string s = lower_copy(s_original);
    const size_t key_len = key_lower.size();

    bool found = false;
    size_t pos = 0;

    while (true) {
        pos = s.find(key_lower, pos);
        if (pos == std::string::npos) break;

        if (!is_key_boundary_left_ci(s, pos)) { pos += 1; continue; }

        size_t i = pos + key_len;
        if (i >= s.size()) { pos += 1; continue; }

        while (i < s.size() && std::isspace((unsigned char)s[i])) ++i;
        if (i < s.size() && s[i] == '=') ++i;
        while (i < s.size() && std::isspace((unsigned char)s[i])) ++i;
        if (i >= s.size()) { pos += 1; continue; }

        if (!is_number_start((unsigned char)s[i])) { pos += 1; continue; }

        const char* start = s_original.c_str() + i;
        char* endp = nullptr;
        double v = std::strtod(start, &endp);
        if (endp == start) { pos += 1; continue; }

        value_out = v;
        found = true;

        pos += 1;
    }

    return found;
}

bool parse_resume_std(const std::string& s, int& std_val)
{
    double v = 0.0;
    if (!extract_last_number_after_key_ci(s, "std", v)) return false;

    const double vr = std::round(v);
    if (std::abs(v - vr) > 1e-12) return false;

    std_val = (int)vr;
    return true;
}

bool parse_resume_D(const std::string& s, double& D_val)
{
    double v = 0.0;
    if (!extract_last_number_after_key_ci(s, "d", v)) return false;

    D_val = v;
    return true;
}

bool validate_resume_run_dir(const std::string& resume_run_dir,
                             const SimParams& P,
                             std::string& err_msg)
{
    err_msg.clear();

    int std_folder = 0;
    double D_folder = 0.0;
    bool aniso_folder = false;

    bool ok_std   = parse_resume_std(resume_run_dir, std_folder);
    bool ok_D     = parse_resume_D(resume_run_dir, D_folder);
    bool ok_aniso = parse_resume_aniso(resume_run_dir, aniso_folder);

    bool any_error = false;

    if (!ok_std || !ok_D || !ok_aniso) {
        any_error = true;
        err_msg += "Could not parse resume_run_dir tokens:\n";
        if (!ok_std)   err_msg += "  - missing/invalid token: std=<int>\n";
        if (!ok_D)     err_msg += "  - missing/invalid token: D=<number>\n";
        if (!ok_aniso) err_msg += "  - missing token: iso or aniso\n";
        err_msg += "Expected styles like: \"std=2, D=0.1, aniso\" OR \"std2_D0.1_aniso\" (spaces ok)\n";
    }

    const int std_param = static_cast<int>(std::round(P.stdev));
    if (ok_std && (std_folder != std_param)) {
        any_error = true;
        err_msg += "stdev mismatch:\n";
        err_msg += "  Folder std = " + std::to_string(std_folder) + "\n";
        err_msg += "  Param  std = " + std::to_string(std_param) + "\n";
    }

    const double tol = 1e-12;
    if (ok_D && (std::abs(D_folder - P.Diffusion_coefficient) > tol)) {
        any_error = true;
        std::ostringstream os;
        os << "Diffusion coefficient mismatch:\n"
           << "  Folder D = " << D_folder << "\n"
           << "  Param  D = " << P.Diffusion_coefficient << "\n";
        err_msg += os.str();
    }

    const bool aniso_param = (std::abs(P.correlation_ls_x - P.correlation_ls_y) > 1e-12);
    if (ok_aniso && (aniso_folder != aniso_param)) {
        any_error = true;
        err_msg += "iso/aniso mismatch:\n";
        err_msg += std::string("  Folder = ") + (aniso_folder ? "aniso" : "iso") + "\n";
        err_msg += std::string("  Param  = ") + (aniso_param  ? "aniso" : "iso") + "\n";
        err_msg += "  (Rule: correlation_ls_x == correlation_ls_y => iso, else aniso)\n";
    }

    if (any_error) {
        err_msg += "resume_run_dir = " + resume_run_dir + "\n";
    }

    return !any_error;
}

// ============================================================
// run folder tag helpers
// ============================================================
std::string fmt_compact_double_tag(double v)
{
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os << std::setprecision(6) << v;

    std::string s = os.str();

    if (s.find('.') != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    if (s.empty()) s = "0";
    return s;
}

std::string make_run_tag_std_D_aniso_df(const SimParams& P)
{
    const int std_int = static_cast<int>(std::round(P.stdev));

    const bool aniso = (std::abs(P.correlation_ls_x - P.correlation_ls_y) > 1e-12);
    const std::string aniso_tag = aniso ? "aniso" : "iso";

    const std::string Dtag  = fmt_compact_double_tag(P.Diffusion_coefficient);
    const std::string DFtag = fmt_compact_double_tag(P.diffusion_factor);

    return "std" + std::to_string(std_int) + "_D" + Dtag + "_" + aniso_tag + "_df" + DFtag;
}

double mean_of(const std::vector<double>& v)
{
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / (double)v.size();
}

// ============================================================
// NaN-safe accumulation helpers
// ============================================================
void accumulate_sum_count(const std::vector<double>& x,
                          std::vector<double>& sum,
                          std::vector<int>& count)
{
    if (sum.empty()) {
        sum.assign(x.size(), 0.0);
        count.assign(x.size(), 0);
    }
    if (x.size() != sum.size()) return;

    for (size_t i = 0; i < x.size(); ++i) {
        double v = x[i];
        if (std::isfinite(v)) {
            sum[i] += v;
            count[i] += 1;
        }
    }
}

std::vector<double> finalize_mean_vec(const std::vector<double>& sum,
                                      const std::vector<int>& count)
{
    std::vector<double> out(sum.size(), std::numeric_limits<double>::quiet_NaN());
    if (sum.size() != count.size()) return out;

    for (size_t i = 0; i < sum.size(); ++i) {
        if (count[i] > 0) out[i] = sum[i] / (double)count[i];
    }
    return out;
}

bool parse_range3(const std::string& s, double& a, double& b, double& c)
{
    // expects "min:max:step"
    auto p1 = s.find(':');
    if (p1 == std::string::npos) return false;

    auto p2 = s.find(':', p1 + 1);
    if (p2 == std::string::npos) return false;

    std::string s1 = s.substr(0, p1);
    std::string s2 = s.substr(p1 + 1, p2 - (p1 + 1));
    std::string s3 = s.substr(p2 + 1);

    return try_parse_double(s1, a)
        && try_parse_double(s2, b)
        && try_parse_double(s3, c);
}

// ============================================================
// delimiter-robust parsing
// ============================================================
static inline std::string trim_copy(std::string s)
{
    auto notsp = [](unsigned char c){ return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notsp));
    s.erase(std::find_if(s.rbegin(), s.rend(), notsp).base(), s.end());
    return s;
}

char detect_delim(const std::string& header)
{
    size_t c1 = std::count(header.begin(), header.end(), ',');
    size_t c2 = std::count(header.begin(), header.end(), '\t');
    size_t c3 = std::count(header.begin(), header.end(), ';');

    if (c2 >= c1 && c2 >= c3) return '\t';
    if (c3 >= c1 && c3 >= c2) return ';';
    return ',';
}

std::vector<std::string> split_line_delim(const std::string& line, char delim)
{
    std::vector<std::string> out;
    std::string cur;
    std::stringstream ss(line);
    while (std::getline(ss, cur, delim)) out.push_back(trim_copy(cur));
    return out;
}

bool try_parse_double(const std::string& s, double& v)
{
    std::string t = trim_copy(s);
    if (t.empty()) return false;
    char* endp = nullptr;
    v = std::strtod(t.c_str(), &endp);
    if (endp == t.c_str()) return false;
    return true;
}

// ============================================================
// CSV utilities: read a "t + columns" table
// ============================================================
bool read_time_series_table_csv(const std::string& path,
                                std::vector<double>& times,
                                std::vector<std::string>& colnames,
                                std::vector<std::vector<double>>& cols)
{
    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    char delim = detect_delim(header);
    auto h = split_line_delim(header, delim);
    if (h.size() < 2) return false;

    colnames.clear();
    cols.clear();

    for (size_t j = 1; j < h.size(); ++j) colnames.push_back(h[j]);
    cols.assign(colnames.size(), std::vector<double>{});
    times.clear();

    std::string line;
    while (std::getline(f, line)) {
        if (trim_copy(line).empty()) continue;
        auto a = split_line_delim(line, delim);
        if (a.size() < 2) continue;

        double t;
        if (!try_parse_double(a[0], t)) continue;

        times.push_back(t);

        for (size_t j = 1; j < a.size() && (j-1) < cols.size(); ++j) {
            double v = std::numeric_limits<double>::quiet_NaN();
            try_parse_double(a[j], v);
            cols[j-1].push_back(v);
        }

        for (size_t j = a.size(); j < h.size(); ++j) {
            cols[j-1].push_back(std::numeric_limits<double>::quiet_NaN());
        }
    }

    for (auto& c : cols) {
        if (c.size() != times.size()) return false;
    }
    return true;
}

bool write_comparison_csv(const std::string& out_path,
                          const std::vector<double>& base_time,
                          const std::vector<std::string>& out_colnames,
                          const std::vector<std::vector<double>>& out_cols)
{
    if (out_colnames.size() != out_cols.size()) return false;
    for (auto& c : out_cols) if (c.size() != base_time.size()) return false;

    std::ofstream f(out_path.c_str());
    if (!f) return false;

    f << "t";
    for (auto& n : out_colnames) f << "," << n;
    f << "\n";

    f << std::setprecision(15);
    for (size_t i = 0; i < base_time.size(); ++i) {
        f << base_time[i];
        for (size_t j = 0; j < out_cols.size(); ++j) {
            double v = out_cols[j][i];
            if (std::isfinite(v)) f << "," << v;
            else f << ",";
        }
        f << "\n";
    }
    return true;
}

// ============================================================
// resampling
// ============================================================
double lerp(double x0, double y0, double x1, double y1, double x)
{
    if (x1 == x0) return y0;
    double a = (x - x0) / (x1 - x0);
    return y0 + a * (y1 - y0);
}

std::vector<double> resample_col_linear(const std::vector<double>& t_src,
                                        const std::vector<double>& y_src,
                                        const std::vector<double>& t_dst)
{
    std::vector<double> out(t_dst.size(), std::numeric_limits<double>::quiet_NaN());
    if (t_src.empty() || y_src.empty() || t_src.size() != y_src.size()) return out;

    size_t k = 0;
    for (size_t i = 0; i < t_dst.size(); ++i) {
        double td = t_dst[i];

        while (k + 1 < t_src.size() && t_src[k+1] < td) k++;

        if (td < t_src.front() || td > t_src.back()) {
            out[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        if (k + 1 >= t_src.size()) {
            out[i] = y_src.back();
            continue;
        }

        double t0 = t_src[k], t1 = t_src[k+1];
        double y0 = y_src[k], y1 = y_src[k+1];

        if (!std::isfinite(y0) || !std::isfinite(y1)) {
            out[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        out[i] = lerp(t0, y0, t1, y1, td);
    }
    return out;
}

void resample_table_linear(const std::vector<double>& t_src,
                           const std::vector<std::vector<double>>& cols_src,
                           const std::vector<double>& t_dst,
                           std::vector<std::vector<double>>& cols_dst)
{
    cols_dst.clear();
    cols_dst.reserve(cols_src.size());
    for (auto& c : cols_src) cols_dst.push_back(resample_col_linear(t_src, c, t_dst));
}

// ============================================================
// FineMean utilities (legacy)
// ============================================================
void accumulate_sum(const std::vector<double>& x,
                    std::vector<double>& sum,
                    std::vector<int>& n)
{
    accumulate_sum_count(x, sum, n);
}

std::vector<double> finalize_mean(const std::vector<double>& sum,
                                  const std::vector<int>& n)
{
    return finalize_mean_vec(sum, n);
}

// ============================================================
// resume/upscale-only readers
// ============================================================
bool parse_keyval_file(const std::string& path,
                       std::map<std::string, std::string>& kv)
{
    kv.clear();
    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;
        auto p = line.find('=');
        if (p == std::string::npos) continue;
        std::string k = trim_copy(line.substr(0, p));
        std::string v = trim_copy(line.substr(p + 1));
        if (!k.empty()) kv[k] = v;
    }
    return !kv.empty();
}

bool read_mean_params_txt(const std::string& path,
                          double& lc_mean, double& lx_mean, double& ly_mean, double& dt_mean,
                          double& nu_x_mean, double& nu_y_mean)
{
    std::map<std::string,std::string> kv;
    if (!parse_keyval_file(path, kv)) return false;

    auto getd = [&](const std::string& k, double& out)->bool{
        auto it = kv.find(k);
        if (it == kv.end()) return false;
        out = std::atof(it->second.c_str());
        return true;
    };

    bool ok =
        getd("lc_mean", lc_mean) &&
        getd("lambda_x_mean", lx_mean) &&
        getd("lambda_y_mean", ly_mean) &&
        getd("dt_mean", dt_mean);

    if (!getd("matern_nu_x_mean", nu_x_mean)) nu_x_mean = 1.5;
    if (!getd("matern_nu_y_mean", nu_y_mean)) nu_y_mean = 1.5;

    return ok;
}

bool read_mean_inverse_cdf_csv(const std::string& path,
                               TimeSeries<double>& mean_cdf)
{
    mean_cdf.readfile(path);
    return !mean_cdf.empty();
}

bool read_xy_table(const std::string& path, std::vector<double>& x, std::vector<double>& y)
{
    x.clear(); y.clear();
    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string line;
    while (std::getline(f, line)) {
        if (trim_copy(line).empty()) continue;

        char delim = ',';
        if (line.find('\t') != std::string::npos) delim = '\t';
        else if (line.find(',') != std::string::npos) delim = ',';
        else delim = ' ';

        auto a = split_line_delim(line, delim);
        if (a.size() < 2) continue;

        double xx, yy;
        if (!try_parse_double(a[0], xx)) continue;
        if (!try_parse_double(a[1], yy)) continue;
        x.push_back(xx);
        y.push_back(yy);
    }
    return !x.empty();
}

double interp1_linear(const std::vector<double>& x,
                      const std::vector<double>& y,
                      double xv)
{
    if (x.empty() || y.empty() || x.size() != y.size()) return 0.0;
    if (xv <= x.front()) return y.front();
    if (xv >= x.back()) return y.back();

    size_t k = 0;
    while (k + 1 < x.size() && x[k+1] < xv) k++;
    if (k + 1 >= x.size()) return y.back();

    return lerp(x[k], y[k], x[k+1], y[k+1], xv);
}

bool read_fine_params_all_csv(const std::string& path,
                              std::vector<double>& lc,
                              std::vector<double>& lx,
                              std::vector<double>& ly,
                              std::vector<double>& dt)
{
    lc.clear(); lx.clear(); ly.clear(); dt.clear();
    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    std::string line;
    while (std::getline(f, line)) {
        if (trim_copy(line).empty()) continue;
        auto a = split_line_delim(line, ',');
        if (a.size() < 5) continue;

        double v1,v2,v3,v4;
        if (!try_parse_double(a[1], v1)) continue;
        if (!try_parse_double(a[2], v2)) continue;
        if (!try_parse_double(a[3], v3)) continue;
        if (!try_parse_double(a[4], v4)) continue;

        lc.push_back(v1);
        lx.push_back(v2);
        ly.push_back(v3);
        dt.push_back(v4);
    }

    return !lc.empty();
}

std::vector<std::pair<int,std::string>> list_fine_folders(const std::string& run_dir)
{
    std::vector<std::pair<int,std::string>> out;

    DIR* d = opendir(run_dir.c_str());
    if (!d) return out;

    struct dirent* ent = nullptr;
    while ((ent = readdir(d)) != nullptr) {
        std::string name = ent->d_name;
        if (name == "." || name == "..") continue;

        if (name.rfind("fine_r", 0) != 0) continue;

        if (name.size() < 10) continue;
        std::string digits = name.substr(std::string("fine_r").size(), 5);
        digits.erase(std::remove_if(digits.begin(), digits.end(),
                                    [](unsigned char c){ return !std::isdigit(c); }),
                     digits.end());
        if (digits.empty()) continue;

        int r = std::atoi(digits.c_str());
        std::string full = joinPath(run_dir, name);
        if (dirExists(full)) out.push_back({r, full});
    }
    closedir(d);

    std::sort(out.begin(), out.end(), [](auto& a, auto& b){ return a.first < b.first; });
    return out;
}

std::string fine_qx_cdf_path(const std::string& fine_dir, int r)
{
    std::string fd = fine_dir;
    if (!fd.empty() && fd.back() != '/' && fd.back() != '\\') fd += "/";
    std::string pfx = makeRealLabel(r) + "_";
    return joinPath(fd, pfx + "qx_inverse_cdf.txt");
}

bool accumulate_inverse_cdf_on_grid(const std::string& fine_dir, int r,
                                    double du, int nU,
                                    std::vector<double>& invcdf_sum)
{
    std::vector<double> x, y;
    std::string path = fine_qx_cdf_path(fine_dir, r);
    if (!fileExists(path)) return false;
    if (!read_xy_table(path, x, y)) return false;

    if ((int)invcdf_sum.size() != nU) invcdf_sum.assign(nU, 0.0);

    for (int k = 0; k < nU; ++k) {
        double u = k * du;
        invcdf_sum[k] += interp1_linear(x, y, u);
    }
    return true;
}

// ============================================================
// BTC calibration helpers
// ============================================================

// Find column index by case-insensitive exact match
static int find_col_ci_exact(const std::vector<std::string>& cols, const std::string& target)
{
    auto low = [](std::string s){
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return (unsigned char)std::tolower(c); });
        return s;
    };
    const std::string t = low(trim_copy(target));

    for (int i = 0; i < (int)cols.size(); ++i) {
        if (low(trim_copy(cols[i])) == t) return i;
    }
    return -1;
}

// Find the "time" column (prefers t/time/sec/seconds)
static int find_time_col_ci(const std::vector<std::string>& cols)
{
    int idx = -1;
    idx = find_col_ci_exact(cols, "t");       if (idx >= 0) return idx;
    idx = find_col_ci_exact(cols, "time");    if (idx >= 0) return idx;
    idx = find_col_ci_exact(cols, "sec");     if (idx >= 0) return idx;
    idx = find_col_ci_exact(cols, "seconds"); if (idx >= 0) return idx;
    return 0; // fallback: first column
}

// Reads RED curve (Upscaled column) from x=0.50BTC_Compare.csv etc.
// Picks first column containing "Upscaled" but NOT "Deriv".
bool read_upscaled_from_btc_compare(
    const std::string& path,
    std::vector<double>& t_out,
    std::vector<double>& c_out
)
{
    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    char delim = detect_delim(header);
    auto cols = split_line_delim(header, delim);
    if (cols.size() < 2) return false;

    int t_col = find_time_col_ci(cols);

    int c_col = -1;
    for (int i = 0; i < (int)cols.size(); ++i) {
        std::string s = trim_copy(cols[i]);
        std::string low = s;
        std::transform(low.begin(), low.end(), low.begin(),
                       [](unsigned char c){ return (unsigned char)std::tolower(c); });

        if (low.find("upscaled") == std::string::npos) continue;
        if (low.find("deriv") != std::string::npos) continue;
        c_col = i;
        break;
    }

    if (t_col < 0 || c_col < 0) return false;

    t_out.clear();
    c_out.clear();

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        auto a = split_line_delim(line, delim);
        if ((int)a.size() <= std::max(t_col, c_col)) continue;

        double t = std::numeric_limits<double>::quiet_NaN();
        double c = std::numeric_limits<double>::quiet_NaN();
        if (!try_parse_double(a[t_col], t)) continue;
        if (!try_parse_double(a[c_col], c)) continue;

        if (std::isfinite(t) && std::isfinite(c)) {
            t_out.push_back(t);
            c_out.push_back(c);
        }
    }

    return t_out.size() >= 2;
}

double rmse_ignore_nan(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size() || a.empty()) return std::numeric_limits<double>::infinity();

    double sse = 0.0;
    int n = 0;

    for (size_t i = 0; i < a.size(); ++i) {
        double x = a[i];
        double y = b[i];
        if (!std::isfinite(x) || !std::isfinite(y)) continue;
        double d = x - y;
        sse += d * d;
        n++;
    }

    if (n == 0) return std::numeric_limits<double>::infinity();
    return std::sqrt(sse / (double)n);
}

// ------------------------------------------------------------
// Read (t, last_column) from CSV (skip header), like gnuplot ncol
// ------------------------------------------------------------
static bool read_lastcol_series_csv(const std::string& path,
                                   std::vector<double>& t_out,
                                   std::vector<double>& y_out)
{
    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    char delim = detect_delim(header);
    auto h = split_line_delim(header, delim);
    if (h.size() < 2) return false;

    const int t_col = 0;
    const int y_col = (int)h.size() - 1; // LAST COLUMN

    t_out.clear();
    y_out.clear();

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        auto a = split_line_delim(line, delim);
        if ((int)a.size() <= y_col) continue;

        double t = std::numeric_limits<double>::quiet_NaN();
        double y = std::numeric_limits<double>::quiet_NaN();
        if (!try_parse_double(a[t_col], t)) continue;
        if (!try_parse_double(a[y_col], y)) continue;

        if (std::isfinite(t) && std::isfinite(y)) {
            t_out.push_back(t);
            y_out.push_back(y);
        }
    }

    return t_out.size() >= 2;
}

// New scoring API: mean log-RMSE over explicit xLocations
// BLACK from BTC_mean.csv (paired columns), RED from per-x BTC_Compare.csv (LAST column, matches gnuplot)
double score_upscaled_vs_black_mean_from_compare(
    const std::string& black_mean_csv,
    const std::string& run_dir,
    const std::vector<double>& xLocations)
{
    using Series = std::pair<std::vector<double>, std::vector<double>>; // (t, y)

    auto trim2 = [](const std::string& s)->std::string{
        auto isspace_ = [](unsigned char c){ return std::isspace(c); };
        size_t b = 0;
        while (b < s.size() && isspace_((unsigned char)s[b])) ++b;
        size_t e = s.size();
        while (e > b && isspace_((unsigned char)s[e-1])) --e;
        return s.substr(b, e - b);
    };

    auto to_lower = [](std::string s)->std::string{
        for (char& c : s) c = (char)std::tolower((unsigned char)c);
        return s;
    };

    auto split_any = [&](const std::string& line, char delim)->std::vector<std::string>{
        std::vector<std::string> out;
        std::string cur;
        std::stringstream ss(line);
        while (std::getline(ss, cur, delim)) out.push_back(trim2(cur));
        return out;
    };

    auto detect_delim_local = [&](const std::string& header)->char{
        if (header.find(',') != std::string::npos) return ',';
        if (header.find('\t') != std::string::npos) return '\t';
        return ' ';
    };

    auto read_any_series_csv = [&](const std::string& path,
                                  std::map<std::string, Series>& series_out,
                                  std::vector<double>& shared_t)->bool
    {
        std::ifstream f(path);
        if (!f) return false;

        series_out.clear();
        shared_t.clear();

        std::string header;
        if (!std::getline(f, header)) return false;
        char delim = detect_delim_local(header);

        std::vector<std::string> cols = split_any(header, delim);
        if (cols.size() < 2) return false;

        int tcount = 0;
        for (auto& c : cols) if (to_lower(trim2(c)) == "t") ++tcount;
        const bool is_paired = (tcount >= 2) && (cols.size() % 2 == 0);

        if (is_paired) {
            const int nPairs = (int)cols.size() / 2;
            std::vector<std::vector<double>> tcols(nPairs), ycols(nPairs);
            std::string line;
            while (std::getline(f, line)) {
                if (trim2(line).empty()) continue;
                std::vector<std::string> tok = split_any(line, delim);
                if ((int)tok.size() < 2*nPairs) continue;

                for (int k = 0; k < nPairs; ++k) {
                    double t=0.0, y=0.0;
                    if (!try_parse_double(tok[2*k], t)) continue;
                    if (!try_parse_double(tok[2*k+1], y)) continue;
                    tcols[k].push_back(t);
                    ycols[k].push_back(y);
                }
            }
            for (int k = 0; k < nPairs; ++k) {
                const std::string name = trim2(cols[2*k+1]);
                series_out[name] = { std::move(tcols[k]), std::move(ycols[k]) };
            }
            return !series_out.empty();
        } else {
            std::vector<std::string> names(cols.begin()+1, cols.end());
            std::vector<std::vector<double>> ycols(names.size());
            std::string line;
            while (std::getline(f, line)) {
                if (trim2(line).empty()) continue;
                std::vector<std::string> tok = split_any(line, delim);
                if (tok.size() < 1 + names.size()) continue;

                double t=0.0;
                if (!try_parse_double(tok[0], t)) continue;
                shared_t.push_back(t);

                for (size_t j=0; j<names.size(); ++j) {
                    double y=std::numeric_limits<double>::quiet_NaN();
                    try_parse_double(tok[1+j], y);
                    ycols[j].push_back(y);
                }
            }
            if (shared_t.empty()) return false;
            for (size_t j=0; j<names.size(); ++j) {
                series_out[trim2(names[j])] = { shared_t, ycols[j] };
            }
            return !series_out.empty();
        }
    };

    auto find_series_for_x = [&](const std::map<std::string, Series>& M, double x)->const Series*{
        std::ostringstream ss;
        ss << "x=" << std::fixed << std::setprecision(2) << x;
        const std::string want = to_lower(ss.str());
        for (const auto& kv : M) {
            std::string k = to_lower(trim2(kv.first));
            if (k.find(want) != std::string::npos) return &kv.second;
        }
        return nullptr;
    };

    // load BLACK mean series
    std::map<std::string, Series> black_map;
    std::vector<double> black_shared_t;
    if (!read_any_series_csv(black_mean_csv, black_map, black_shared_t)) {
        std::cerr << "WARNING: could not read black mean CSV: " << black_mean_csv << "\n";
        return std::numeric_limits<double>::infinity();
    }

    // score: log-RMSE over all x locations using per-x compare files
    const double eps = 1e-12;
    const double tmin = 2.0;   // tail emphasis; try 3.0 if you want more tail
    const double tmax = 1e100;

    double sum_rmse = 0.0;
    int n_used = 0;

    for (double x : xLocations) {
        const Series* B = find_series_for_x(black_map, x);
        if (!B) {
            std::cerr << "WARNING: missing BLACK series for x=" << x << " in " << black_mean_csv << "\n";
            return std::numeric_limits<double>::infinity();
        }

        const std::string cmp_path = joinPath(run_dir, fmt_x(x) + "BTC_Compare.csv");

        std::vector<double> tu, cu;
        if (!fileExists(cmp_path) || !read_lastcol_series_csv(cmp_path, tu, cu)) {
            std::cerr << "WARNING: could not read compare (upscaled) file: " << cmp_path << "\n";
            return std::numeric_limits<double>::infinity();
        }

        const std::vector<double>& tb = B->first;
        const std::vector<double>& cb = B->second;
        if (tb.size() < 2 || cb.size() != tb.size()) return std::numeric_limits<double>::infinity();

        std::vector<double> cu_i = resample_col_linear(tu, cu, tb);

        double sse = 0.0;
        int n = 0;
        for (size_t i = 0; i < tb.size() && i < cb.size() && i < cu_i.size(); ++i) {
            const double t = tb[i];
            if (!(t >= tmin && t <= tmax)) continue;

            const double b = cb[i];
            const double u = cu_i[i];
            if (!std::isfinite(b) || !std::isfinite(u)) continue;

            const double lb = std::log10(std::max(b, eps));
            const double lu = std::log10(std::max(u, eps));
            const double d = lb - lu;

            sse += d * d;
            ++n;
        }
        if (n < 10) return std::numeric_limits<double>::infinity();

        const double rmse = std::sqrt(sse / (double)n);
        sum_rmse += rmse;
        ++n_used;
    }

    if (n_used == 0) return std::numeric_limits<double>::infinity();
    return sum_rmse / (double)n_used;
}

double score_upscaled_vs_black_mean_from_compare(
    const std::string& black_mean_csv,
    const std::string& run_dir)
{
    return score_upscaled_vs_black_mean_from_compare(black_mean_csv, run_dir, {0.5, 1.5, 2.5});
}

bool list_calibration_runs_with_df(
    const std::string& root_dir,
    std::vector<double>& dfs_out,
    std::vector<std::string>& run_dirs_out)
{
    dfs_out.clear();
    run_dirs_out.clear();

    if (!dirExists(root_dir)) return false;

    std::regex re_df(R"((?:^|[_\-])df(?:=)?([0-9]*\.?[0-9]+))", std::regex::icase);

    for (const auto& ent : std::filesystem::directory_iterator(root_dir)) {
        if (!ent.is_directory()) continue;
        const std::string name = ent.path().filename().string();
        std::smatch m;
        if (!std::regex_search(name, m, re_df)) continue;

        double df = std::numeric_limits<double>::quiet_NaN();
        if (!try_parse_double(m[1].str(), df)) continue;

        dfs_out.push_back(df);
        run_dirs_out.push_back(ent.path().string());
    }

    std::vector<size_t> idxs(dfs_out.size());
    std::iota(idxs.begin(), idxs.end(), 0);
    std::sort(idxs.begin(), idxs.end(), [&](size_t a, size_t b){ return dfs_out[a] < dfs_out[b]; });

    std::vector<double> dfs2;
    std::vector<std::string> dirs2;
    dfs2.reserve(idxs.size());
    dirs2.reserve(idxs.size());
    for (size_t k : idxs) {
        dfs2.push_back(dfs_out[k]);
        dirs2.push_back(run_dirs_out[k]);
    }
    dfs_out.swap(dfs2);
    run_dirs_out.swap(dirs2);
    return !dfs_out.empty();
}

// ============================================================
// gnuplot helper utilities (MISSING BEFORE - NOW IMPLEMENTED)
// ============================================================
bool starts_with(const std::string& s, const std::string& prefix)
{
    if (prefix.size() > s.size()) return false;
    return std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::string sanitize_token(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s) {
        if (std::isalnum(c) || c == '_' || c == '-' || c == '.') out.push_back((char)c);
        else out.push_back('_');
    }
    // collapse repeated underscores
    std::string out2;
    out2.reserve(out.size());
    char prev = '\0';
    for (char c : out) {
        if (c == '_' && prev == '_') continue;
        out2.push_back(c);
        prev = c;
    }
    // trim underscores
    while (!out2.empty() && out2.front() == '_') out2.erase(out2.begin());
    while (!out2.empty() && out2.back()  == '_') out2.pop_back();
    if (out2.empty()) out2 = "plot";
    return out2;
}

// Group columns by the "base" (typically x=0.50) regardless of prefix
std::string base_name_from_header(const std::string& header)
{
    std::string h = trim_copy(header);

    // remove known prefixes
    const char* prefixes[] = {
        "FineDerivMean_", "FineDeriv_", "UpscaledDeriv_",
        "FineMean_", "Fine_", "Upscaled_",
        "UpscaledDeriv", "Upscaled", "FineDerivMean", "FineDeriv", "FineMean", "Fine"
    };
    for (const char* p : prefixes) {
        const std::string P(p);
        if (starts_with(h, P)) {
            h = h.substr(P.size());
            if (!h.empty() && h[0] == '_') h.erase(h.begin());
            break;
        }
    }

    // common patterns: "...x=0.50..." -> keep from x=
    auto pos = h.find("x=");
    if (pos != std::string::npos) return trim_copy(h.substr(pos));

    // else: take last token after underscore (fallback)
    auto us = h.rfind('_');
    if (us != std::string::npos && us + 1 < h.size()) return trim_copy(h.substr(us + 1));

    return h;
}

bool read_csv_header(const std::string& csv_path, std::vector<std::string>& headers)
{
    headers.clear();
    std::ifstream f(csv_path.c_str());
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    if (trim_copy(header).empty()) return false;

    char delim = detect_delim(header);
    headers = split_line_delim(header, delim);
    return !headers.empty();
}

// ============================================================
// gnuplot helpers
// ============================================================
bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label,
                                                bool skip_base_t)
{
    std::vector<std::string> headers;
    if (!read_csv_header(csv_path, headers)) return false;

    std::map<std::string, std::vector<std::pair<int,std::string>>> groups;

    // columns are 1-based in gnuplot using "using 1:col"
    for (int c = 2; c <= (int)headers.size(); ++c) {
        const std::string& h = headers[c-1];

        if (skip_base_t) {
            std::string ht = trim_copy(h);
            if (ht == "t" || ht == "T" || ht == "time") continue;
        }

        std::string base = base_name_from_header(h);
        groups[base].push_back({c, h});
    }

    if (groups.empty()) return false;

    std::ofstream gp(gp_path.c_str());
    if (!gp) return false;

    gp << "set datafile separator ','\n";
    gp << "set grid\n";
    gp << "set xlabel 't'\n";
    gp << "set ylabel '" << y_label << "'\n";
    gp << "set key outside\n";
    gp << "set term pngcairo size 1400,900\n";

    gp << "set style line 1 lc rgb '#9a9a9a' lt 1 lw 1\n"; // fine gray
    gp << "set style line 2 lc rgb '#000000' lt 1 lw 3\n"; // mean bold black
    gp << "set style line 3 lc rgb '#d60000' lt 1 lw 3\n"; // upscaled bold red

    for (auto& kv : groups) {
        const std::string& base = kv.first;
        const auto& cols = kv.second;

        std::string out_png = fig_prefix + "_" + sanitize_token(base) + ".png";
        gp << "set output '" << out_png << "'\n";
        gp << "set title '" << base << "'\n";

        gp << "plot ";
        bool first = true;

        for (auto& pr : cols) {
            int col = pr.first;
            const std::string& name = pr.second;

            int ls = 1;
            if (starts_with(name, "FineMean_") || starts_with(name, "FineDerivMean_")) ls = 2;
            else if (starts_with(name, "Upscaled") || starts_with(name, "UpscaledDeriv")) ls = 3;
            else if (starts_with(name, "Fine_") || starts_with(name, "FineDeriv_")) ls = 1;

            if (!first) gp << ", ";
            first = false;

            gp << "'" << csv_path << "' using 1:" << col
               << " with lines ls " << ls
               << " title '" << name << "'";
        }

        gp << "\n";
    }

    gp << "unset output\n";
    return true;
}

int run_gnuplot_script(const std::string& gp_path)
{
    std::string cmd = "gnuplot \"" + gp_path + "\"";
    return std::system(cmd.c_str());
}
