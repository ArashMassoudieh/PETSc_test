#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>

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
