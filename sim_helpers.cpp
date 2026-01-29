// sim_helpers.cpp
#include "sim_helpers.h"

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
    // choose the most frequent among [',', '\t', ';']
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

    // First column assumed time
    colnames.clear();
    cols.clear();

    // keep only non-time columns
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
            try_parse_double(a[j], v); // if fails, v stays NaN
            cols[j-1].push_back(v);
        }

        // if short row, pad NaNs
        for (size_t j = a.size(); j < h.size(); ++j) {
            cols[j-1].push_back(std::numeric_limits<double>::quiet_NaN());
        }
    }

    // sanity: all columns same length
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

    // assume t_src is increasing (as produced by your solver)
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

        // if endpoints NaN, keep NaN (simple & safe)
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
                          double& lc_mean, double& lx_mean, double& ly_mean, double& dt_mean)
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
        // tolerate commas/tabs/spaces
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

        // a[0] realization (ignored)
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

        // expecting fine_r0001, fine_r0002, ...
        if (name.rfind("fine_r", 0) != 0) continue;

        // parse r#### part
        if (name.size() < 10) continue; // "fine_r" + 4 digits
        std::string digits = name.substr(std::string("fine_r").size(), 5); // like "0001" or "0001/"? (safe)
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
