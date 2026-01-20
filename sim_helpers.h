#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <limits>
#include <petscsys.h>   // for PetscInt
#include <iomanip>

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

// Non-inline: implemented in helpers.cpp
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
// Reads table with header; first column is time.
// Works for comma/tab/semicolon/space separated.
// Returns times, colnames (excluding time), cols (each column vector).
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
// FineMean utilities (mean of what is on disk; no std)
// --------------------
// Accumulate finite values into sum/count arrays (size matched to x).
void accumulate_sum(
    const std::vector<double>& x,
    std::vector<double>& sum,
    std::vector<int>& n
);

// Finalize mean = sum / n (NaN where n==0).
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
    double& lc_mean, double& lx_mean, double& ly_mean, double& dt_mean
);

// mean_qx_inverse_cdf.txt: header "u,v" then data
bool read_mean_inverse_cdf_csv(
    const std::string& path,
    std::vector<double>& u,
    std::vector<double>& v
);

// 2-column numeric table, any delimiter, no strict header required
bool read_xy_table(const std::string& path, std::vector<double>& x, std::vector<double>& y);

double interp1_linear(const std::vector<double>& x, const std::vector<double>& y, double xv);

// fine_params_all.csv
bool read_fine_params_all_csv(
    const std::string& path,
    std::vector<double>& lc,
    std::vector<double>& lx,
    std::vector<double>& ly,
    std::vector<double>& dt
);

// Scan run_dir for fine_r#### folders => list of (rid, fullpath)
std::vector<std::pair<int,std::string>> list_fine_folders(const std::string& run_dir);

// Build per-realization qx_inverse_cdf filename in fine folder
std::string fine_qx_cdf_path(const std::string& fine_dir, int r);

// Accumulate per-realization qx inverse cdf onto du grid
bool accumulate_inverse_cdf_on_grid(
    const std::string& fine_dir, int r,
    double du, int nU,
    std::vector<double>& invcdf_sum
);

/*
// Writes a python script that plots "Fine realizations (gray) + FineMean (black) + Upscaled_mean (red dashed)"
// for each series column found in the comparison CSV.
// Output: multiple PNGs saved next to the CSV.
bool write_btc_compare_plot_py(const std::string& py_path,
                              const std::string& csv_path,
                              const std::string& fig_prefix,
                              const std::string& y_label);

// Optional: run python script (best-effort)
int run_python_script(const std::string& py_path);

// gnuplot
bool write_btc_compare_plot_gnuplot(const std::string& gp_path,
                                    const std::string& csv_path,
                                    const std::string& fig_prefix,
                                    const std::string& y_label);

int run_gnuplot_script(const std::string& gp_path);

// simple gnuplot
bool write_btc_compare_plot_gnuplot_simple(const std::string& gp_path,
                                           const std::string& csv_path,
                                           const std::string& fig_prefix,
                                           const std::string& y_label,
                                           int max_series = 50);

int run_gnuplot_script(const std::string& gp_path);

// gnuplot by name
bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label);

int run_gnuplot_script(const std::string& gp_path);
*/

// gnuplot by name (one figure per base series: x=0.50, x=1.50, ...)
bool write_btc_compare_plot_gnuplot_by_basename(const std::string& gp_path,
                                                const std::string& csv_path,
                                                const std::string& fig_prefix,
                                                const std::string& y_label,
                                                bool skip_base_t = true);

int run_gnuplot_script(const std::string& gp_path);

inline std::string fmt_x(double x) {
    std::ostringstream ss;
    ss << "x=" << std::fixed << std::setprecision(2) << x;
    return ss.str();
}
