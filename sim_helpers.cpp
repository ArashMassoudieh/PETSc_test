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
// BTC calibration helpers (NEW)
// ============================================================

bool read_btc_mean_paired_csv(
    const std::string& path,
    std::vector<std::string>& loc_names,
    std::vector<std::vector<double>>& t_cols,
    std::vector<std::vector<double>>& c_cols
)
{
    loc_names.clear();
    t_cols.clear();
    c_cols.clear();

    std::ifstream f(path.c_str());
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    char delim = detect_delim(header);
    auto h = split_line_delim(header, delim);
    if (h.size() < 4) return false; // at least t,x , t,x

    // Expect pairs: t, x=..., t, x=..., ...
    const int nPairs = (int)h.size() / 2;
    loc_names.reserve(nPairs);
    t_cols.assign(nPairs, {});
    c_cols.assign(nPairs, {});

    for (int k = 0; k < nPairs; ++k) {
        std::string name = h[2*k + 1];
        name = trim_copy(name);
        loc_names.push_back(name);
    }

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        auto a = split_line_delim(line, delim);
        if ((int)a.size() < 2) continue;

        for (int k = 0; k < nPairs; ++k) {
            const int it = 2*k;
            const int ic = 2*k + 1;
            if (ic >= (int)a.size()) continue;

            double tt = std::numeric_limits<double>::quiet_NaN();
            double cc = std::numeric_limits<double>::quiet_NaN();
            try_parse_double(a[it], tt);
            try_parse_double(a[ic], cc);

            if (std::isfinite(tt)) t_cols[k].push_back(tt);
            else t_cols[k].push_back(std::numeric_limits<double>::quiet_NaN());

            if (std::isfinite(cc)) c_cols[k].push_back(cc);
            else c_cols[k].push_back(std::numeric_limits<double>::quiet_NaN());
        }
    }

    // basic sanity
    for (int k = 0; k < nPairs; ++k) {
        if (t_cols[k].size() < 2) return false;
        if (t_cols[k].size() != c_cols[k].size()) return false;
    }

    return true;
}

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

    int t_col = -1;
    int c_col = -1;

    for (int i = 0; i < (int)cols.size(); ++i) {
        std::string s = cols[i];
        s = trim_copy(s);
        std::string low = s;
        std::transform(low.begin(), low.end(), low.begin(),
                       [](unsigned char c){ return (unsigned char)std::tolower(c); });

        if (low == "t") t_col = i;

        // We want the Upscaled column, header looks like: Upscaledx=0.50 (no underscore)
        if (low.find("upscaled") != std::string::npos) c_col = i;
    }

    if (t_col < 0 || c_col < 0) return false;

    t_out.clear();
    c_out.clear();

    std::string line;
    while (std::getline(f, line)) {
        auto a = split_line_delim(line, delim);
        if ((int)a.size() <= std::max(t_col, c_col)) continue;

        double t = std::numeric_limits<double>::quiet_NaN();
        double c = std::numeric_limits<double>::quiet_NaN();
        try_parse_double(a[t_col], t);
        try_parse_double(a[c_col], c);

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

double score_upscaled_vs_black_mean_from_compare(
    const std::string& black_btc_mean_csv,
    const std::string& run_dir
)
{
    // read BLACK (BTC_mean.csv)
    std::vector<std::string> names;
    std::vector<std::vector<double>> t_black, c_black;

    if (!read_btc_mean_paired_csv(black_btc_mean_csv, names, t_black, c_black)) {
        std::cerr << "ERROR: cannot read black BTC_mean: " << black_btc_mean_csv << "\n";
        return std::numeric_limits<double>::infinity();
    }

    double total = 0.0;
    int used = 0;

    for (int k = 0; k < (int)names.size(); ++k) {
        const std::string& loc = names[k]; // expected "x=0.50"
        const std::string compare_csv = joinPath(run_dir, loc + "BTC_Compare.csv");

        if (!fileExists(compare_csv)) {
            std::cerr << "WARNING: missing compare file: " << compare_csv << "\n";
            continue;
        }

        std::vector<double> t_red, c_red;
        if (!read_upscaled_from_btc_compare(compare_csv, t_red, c_red)) {
            std::cerr << "WARNING: cannot read upscaled curve from: " << compare_csv << "\n";
            continue;
        }

        // Interpolate RED onto BLACK time grid
        std::vector<double> c_red_i = resample_col_linear(t_red, c_red, t_black[k]);

        const double e = rmse_ignore_nan(c_black[k], c_red_i);
        if (std::isfinite(e)) {
            total += e;
            used++;
        }
    }

    if (used == 0) return std::numeric_limits<double>::infinity();
    return total / (double)used;
}


// ============================================================
// Gnuplot helpers
// ============================================================
static inline bool read_csv_header(const std::string& csv_path,
                                  std::vector<std::string>& headers)
{
    std::ifstream f(csv_path.c_str());
    if (!f) return false;

    std::string line;
    if (!std::getline(f, line)) return false;

    char delim = detect_delim(line);
    headers = split_line_delim(line, delim);
    return headers.size() >= 2;
}

static inline bool starts_with(const std::string& s, const std::string& p)
{
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

static inline std::string sanitize_token(std::string s)
{
    for (char& c : s) {
        if (std::isalnum((unsigned char)c)) continue;
        if (c == '_' || c == '-' ) continue;
        c = '_';
    }
    return s;
}

static inline std::string base_name_from_header(const std::string& h)
{
    // expected: Prefix_"x=0.50"
    // base is everything after last underscore
    auto pos = h.rfind('_');
    if (pos == std::string::npos) return h;
    return h.substr(pos + 1);
}

static inline std::string dirname_of(const std::string& path)
{
    size_t p = path.find_last_of("/\\");
    if (p == std::string::npos) return ".";
    return path.substr(0, p);
}

bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label,
                                                bool skip_base_t)
{
    std::vector<std::string> headers;
    if (!read_csv_header(csv_path, headers)) return false;

    // map base -> list of (colIndex, headerName)
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

    // IMPORTANT: outputs must go to run_dir, not build cwd.
    // We force outputs to fig_prefix (which main should pass as joinPath(run_dir,"..."))
    std::string out_prefix = fig_prefix;

    std::ofstream gp(gp_path.c_str());
    if (!gp) return false;

    gp << "set datafile separator ','\n";
    gp << "set grid\n";
    gp << "set xlabel 't'\n";
    gp << "set ylabel '" << y_label << "'\n";
    gp << "set key outside\n";
    gp << "set term pngcairo size 1400,900\n";

    // Line styles (as you requested):
    // Fine_*        : gray
    // FineMean_*    : bold black
    // Upscaled*     : bold red
    // (Derivative versions also respected)
    gp << "set style line 1 lc rgb '#9a9a9a' lt 1 lw 1\n"; // fine gray
    gp << "set style line 2 lc rgb '#000000' lt 1 lw 3\n"; // mean bold black
    gp << "set style line 3 lc rgb '#d60000' lt 1 lw 3\n"; // upscaled bold red\n";

    for (auto& kv : groups) {
        const std::string& base = kv.first;
        const auto& cols = kv.second;

        std::string out_png = out_prefix + "_" + sanitize_token(base) + ".png";
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
