// sim_helpers.h
#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <cctype>

#include <TimeSeries.h>
#include <TimeSeriesSet.h>

#include <petscsys.h>   // PetscInt

// Forward declaration (avoid including sim_runner.h here)
struct SimParams;

// --------------------
// FS / path helpers
// --------------------
bool createDirectory(const std::string& path);
bool fileExists(const std::string& path);
bool dirExists (const std::string& path);
std::string joinPath(const std::string& dir, const std::string& filename);

// --------------------
// naming helpers
// --------------------
std::string makeTimestamp();            // YYYYMMDD_HHMMSS
std::string makeRealLabel(int r1);      // r0001, r0002, ...
std::string makeFineFolder(int r1);     // fine_r0001, ...

// --------------------
// run folder tag helpers
// --------------------
std::string fmt_compact_double_tag(double v);

// Build run tag used for NEW run_dir folder names:
//   std<stdev>_D<diffusion>_<iso|aniso>_df<diffusion_factor>
std::string make_run_tag_std_D_aniso_df(const SimParams& P);

// --------------------------------------------------
// Resume folder consistency helpers
// --------------------------------------------------
bool parse_resume_std(const std::string& s, int& std_val);
bool parse_resume_D  (const std::string& s, double& D_val);
bool parse_resume_aniso(const std::string& s, bool& is_aniso);

bool validate_resume_run_dir(const std::string& resume_run_dir,
                             const SimParams& P,
                             std::string& err_msg);

double mean_of(const std::vector<double>& v);

// --------------------------------------------------
// Grid indexing helpers
// --------------------------------------------------
inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx)
{
    return j * nx + i;
}

inline bool on_bc(PetscInt i, PetscInt j,
                  PetscInt /*nx*/, PetscInt /*ny*/,
                  double x, double y, double Lx, double Ly)
{
    const double eps = 1e-14;
    return (i == 0) || (j == 0) ||
           (std::abs(x - Lx) < eps) ||
           (std::abs(y - Ly) < eps);
}

// --------------------------------------------------
// NaN-safe accumulation helpers
// --------------------------------------------------
inline bool is_finite_number(double v)
{
    return std::isfinite(v);
}

void accumulate_sum_count(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& count
);

std::vector<double> finalize_mean_vec(
    const std::vector<double>& sum,
    const std::vector<int>& count
);

// --------------------
// delimiter-robust parsing
// --------------------
char detect_delim(const std::string& header);
std::vector<std::string> split_line_delim(const std::string& line, char delim);
bool try_parse_double(const std::string& s, double& v);

// --------------------
// CSV utilities
// --------------------
bool read_time_series_table_csv(
    const std::string& path,
    std::vector<double>& times,
    std::vector<std::string>& colnames,
    std::vector<std::vector<double>>& cols
);

bool write_comparison_csv(
    const std::string& out_path,
    const std::vector<double>& base_time,
    const std::vector<std::string>& out_colnames,
    const std::vector<std::vector<double>>& out_cols
);

inline bool load_csv_time_only(const std::string& path, std::vector<double>& t_out)
{
    std::vector<double> t;
    std::vector<std::string> names;
    std::vector<std::vector<double>> cols;
    if (!read_time_series_table_csv(path, t, names, cols)) return false;
    if (t.empty()) return false;
    t_out = std::move(t);
    return true;
}

// Accepts either:
//   u,v
//   0.0,1.2
// or no header:
//   0.0,1.2
// Also accepts whitespace separated "0.0 1.2"
static inline bool read_inverse_cdf_any_format(const std::string& path, TimeSeries<double>& out)
{
    std::ifstream f(path);
    if (!f) return false;

    out.clear();

    std::string line;
    bool header_checked = false;

    auto trim = [](std::string& s) {
        auto isspace_ = [](unsigned char c){ return std::isspace(c); };
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](unsigned char c){ return !isspace_(c); }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [&](unsigned char c){ return !isspace_(c); }).base(), s.end());
    };

    while (std::getline(f, line)) {
        trim(line);
        if (line.empty()) continue;

        if (!header_checked) {
            header_checked = true;
            std::string low = line;
            std::transform(low.begin(), low.end(), low.begin(),
                           [](unsigned char c){ return (unsigned char)std::tolower(c); });
            if ((low.find('u') != std::string::npos) &&
                (low.find('v') != std::string::npos) &&
                (low.find(',') != std::string::npos || low.find(' ') != std::string::npos))
            {
                continue;
            }
        }

        for (char& c : line) if (c == ',') c = ' ';

        std::stringstream ss(line);
        double u = 0.0, v = 0.0;
        if (!(ss >> u >> v)) continue;

        out.append(u, v);
    }

    return out.size() >= 2;
}

// --------------------
// resampling
// --------------------
double lerp(double x0, double y0, double x1, double y1, double x);

std::vector<double> resample_col_linear(
    const std::vector<double>& t_src,
    const std::vector<double>& y_src,
    const std::vector<double>& t_dst
);

void resample_table_linear(
    const std::vector<double>& t_src,
    const std::vector<std::vector<double>>& cols_src,
    const std::vector<double>& t_dst,
    std::vector<std::vector<double>>& cols_dst
);

// --------------------
// FineMean utilities (legacy)
// --------------------
void accumulate_sum(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& n
);

std::vector<double> finalize_mean(
    const std::vector<double>& sum,
    const std::vector<int>& n
);

// --------------------
// resume/upscale-only readers
// --------------------
bool parse_keyval_file(const std::string& path, std::map<std::string, std::string>& kv);

bool read_mean_params_txt(
    const std::string& path,
    double& lc_mean, double& lx_mean, double& ly_mean, double& dt_mean, double& nu_x_mean, double& nu_y_mean
);

bool read_mean_inverse_cdf_csv(
    const std::string& path,
    TimeSeries<double>& mean_cdf
);

bool read_xy_table(const std::string& path, std::vector<double>& x, std::vector<double>& y);

double interp1_linear(const std::vector<double>& x, const std::vector<double>& y, double xv);

bool read_fine_params_all_csv(
    const std::string& path,
    std::vector<double>& lc,
    std::vector<double>& lx,
    std::vector<double>& ly,
    std::vector<double>& dt
);

std::vector<std::pair<int,std::string>> list_fine_folders(const std::string& run_dir);

std::string fine_qx_cdf_path(const std::string& fine_dir, int r);

bool accumulate_inverse_cdf_on_grid(
    const std::string& fine_dir, int r,
    double du, int nU,
    std::vector<double>& invcdf_sum
);

// --------------------
// BTC calibration helpers
// --------------------

// Reads RED curve (Upscaled column) from x=0.50BTC_Compare.csv etc.
bool read_upscaled_from_btc_compare(
    const std::string& path,
    std::vector<double>& t_out,
    std::vector<double>& c_out
);

double rmse_ignore_nan(const std::vector<double>& a, const std::vector<double>& b);

// Main scoring: mean RMSE over selected locations
double score_upscaled_vs_black_mean_from_compare(
    const std::string& black_btc_mean_csv,
    const std::string& run_dir,
    const std::vector<double>& xLocations
);

// Backward compatible: uses x=0.50,1.50,2.50
double score_upscaled_vs_black_mean_from_compare(
    const std::string& black_btc_mean_csv,
    const std::string& run_dir
);

// --------------------
// gnuplot helpers
// --------------------
bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label,
                                                bool skip_base_t = true);

int run_gnuplot_script(const std::string& gp_path);

// --------------------
// formatting
// --------------------
inline std::string fmt_x(double x) {
    std::ostringstream ss;
    ss << "x=" << std::fixed << std::setprecision(2) << x;
    return ss.str();
}
