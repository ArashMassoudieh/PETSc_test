//sim_runner.cpp
#include "sim_runner.h"
#include "sim_helpers.h"

#include <mpi.h>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <sstream>
#include <cctype>
#include <cmath>
#include <limits>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <random>

#include "grid.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline double du_from_nu(int nu)
{
    return (nu > 0) ? (1.0 / double(nu)) : 1.0;
}

// ============================================================================
// Forward declarations for file-local helpers (static)
// ============================================================================

static bool build_mean_for_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_dir,
    const std::string& resume_run_dir,
    FineScaleOutputs& fine_outputs,
    TimeSeries<double>& invcdf_mean,
    int rank);

static bool run_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    double lc_mean,
    double lx_mean,
    double ly_mean,
    double nu_x_mean,
    double nu_y_mean,
    double dt_mean,
    const TimeSeries<double>& invcdf_mean,
    std::string& up_dir_out,
    std::string& up_btc_path_out,
    std::string& up_btc_deriv_path_out,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs_transport_pdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs_transport_cdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_pdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_cdf,
    int rank);

// NEW: rebuild/recovery helpers
static bool rebuild_from_existing_btc_files(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    RunOutputs& out,
    int rank);

// ============================================================================
// prepare_run_dir_mpi
// ============================================================================

std::string prepare_run_dir_mpi(
    const std::string& output_dir,
    const std::string& resume_run_dir,
    const RunOptions& opts,
    int rank,
    const std::string& run_tag)
{
    std::string run_dir;

    if (rank == 0) {
        createDirectory(output_dir);

        // IMPORTANT:
        // In upscale-only mode, ALWAYS create a NEW output run_dir.
        // resume_run_dir is treated as INPUT SOURCE only.
        if (opts.upscale_only) {

            if (!opts.hardcoded_mean) {
                if (resume_run_dir.empty()) {
                    std::cerr << "ERROR: --upscale-only (without --hardcoded-mean) requires --run-dir=/path/to/input_run\n";
                    MPI_Abort(PETSC_COMM_WORLD, 1);
                }
                if (!dirExists(resume_run_dir)) {
                    std::cerr << "ERROR: input resume_run_dir does not exist: " << resume_run_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 2);
                }
            }

            std::string name = "run_" + makeTimestamp();
            if (!run_tag.empty()) name += "_" + run_tag;
            run_dir = joinPath(output_dir, name);
            createDirectory(run_dir);

            if (opts.hardcoded_mean) {
                std::cout << "Hardcoded-mean upscaled-only: created NEW run_dir: "
                          << run_dir << "\n";
            } else {
                std::cout << "Upscale-only: INPUT source = " << resume_run_dir << "\n";
                std::cout << "Upscale-only: OUTPUT run_dir = " << run_dir << "\n";
            }

        } else {
            std::string name = "run_" + makeTimestamp();
            if (!run_tag.empty()) name += "_" + run_tag;
            run_dir = joinPath(output_dir, name);
            createDirectory(run_dir);
            std::cout << "Run directory: " << run_dir << "\n";
        }
    }

    // broadcast run_dir
    int len = 0;
    if (rank == 0) len = (int)run_dir.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    run_dir.resize(len);
    MPI_Bcast(run_dir.data(), len, MPI_CHAR, 0, PETSC_COMM_WORLD);

    if (!run_dir.empty() && run_dir.back() != '/' && run_dir.back() != '\\') run_dir += "/";
    MPI_Barrier(PETSC_COMM_WORLD);

    return run_dir;
}

// ============================================================================
// small helpers
// ============================================================================

static inline void build_flat_invcdf(TimeSeries<double>& invcdf_mean, double qx_const)
{
    invcdf_mean.clear();
    invcdf_mean.append(0.0, qx_const);
    invcdf_mean.append(1.0, qx_const);
}

static inline std::string trim_copy(std::string s)
{
    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

static inline std::string to_lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (unsigned char)std::tolower(c); });
    return s;
}

static bool try_load_hardcoded_invcdf(TimeSeries<double>& invcdf_mean, const std::string& path)
{
    if (path.empty()) return false;
    if (!fileExists(path)) return false;

    TimeSeries<double> tmp;
    if (!read_inverse_cdf_any_format(path, tmp)) return false;
    if (tmp.size() < 2) return false;

    invcdf_mean = tmp;
    return true;
}

static inline std::string dirname_of(std::string path)
{
    const auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return std::string();
    return path.substr(0, pos);
}

static bool build_mean_qx_inverse_cdf_from_multi(
    const std::string& in_multi_path,
    const std::string& out_mean_path
){
    std::ifstream f(in_multi_path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;

    const char delim = detect_delim(header);

    std::vector<double> u_out;
    std::vector<double> v_out;
    u_out.reserve(4096);
    v_out.reserve(4096);

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        const auto tok = split_line_delim(line, delim);
        if (tok.size() < 2) continue;

        const size_t npairs = tok.size() / 2;
        if (npairs < 1) continue;

        double t_ref = 0.0;
        double sum_q = 0.0;
        int    cnt_q = 0;

        for (size_t k = 0; k < npairs; ++k) {
            double t = 0.0, q = 0.0;
            if (!try_parse_double(tok[2*k + 0], t)) continue;
            if (!try_parse_double(tok[2*k + 1], q)) continue;

            if (k == 0) t_ref = t;
            if (std::abs(t - t_ref) > 1e-10) continue;

            if (is_finite_number(q)) {
                sum_q += q;
                cnt_q++;
            }
        }

        if (cnt_q <= 0) continue;

        u_out.push_back(t_ref);
        v_out.push_back(sum_q / (double)cnt_q);
    }

    if (u_out.size() < 2) return false;

    std::ofstream o(out_mean_path);
    if (!o) return false;

    o << "u,v\n";
    o << std::setprecision(15);
    for (size_t i = 0; i < u_out.size(); ++i) {
        o << u_out[i] << "," << v_out[i] << "\n";
    }
    return true;
}

// NEW: file-loading helpers for --btc-from-files mode
static bool read_tsset_file_if_exists(const std::string& path, TimeSeriesSet<double>& out)
{
    if (!fileExists(path)) return false;
    return out.read(path);
}

static std::vector<std::string> list_fine_dirs_sorted(const std::string& run_dir)
{
    namespace fs = std::filesystem;
    std::vector<std::string> out;
    std::error_code ec;

    if (!fs::exists(fs::path(run_dir), ec)) return out;

    for (const auto& ent : fs::directory_iterator(fs::path(run_dir), ec)) {
        if (ec) break;
        if (!ent.is_directory()) continue;

        const std::string name = ent.path().filename().string();
        if (name.rfind("fine_r", 0) == 0) {
            out.push_back(ent.path().string());
        }
    }

    std::sort(out.begin(), out.end());
    return out;
}

static std::string basename_only(const std::string& path)
{
    return std::filesystem::path(path).filename().string();
}

static bool fine_dir_to_rlab(const std::string& fine_dir, std::string& rlab)
{
    const std::string base = basename_only(fine_dir); // fine_r12
    if (base.rfind("fine_r", 0) != 0) return false;

    const std::string suffix = base.substr(std::string("fine_r").size());
    if (suffix.empty()) return false;

    rlab = "r" + suffix;
    return true;
}

// ============================================================================
// mean_ts switch (driven by RunOptions::mean_ts_mode)
//   - First   : reference grid = first series
//   - Longest : reference grid = longest series
//   - Union   : reference grid = union of all times (unique), interpolate each
//               series only inside its support, then average NaN-safely
// ============================================================================

static inline bool ts_has_support_at_t(const TimeSeries<double>& ts, double t)
{
    const int n = (int)ts.size();
    if (n < 2) return false;

    const double t0 = ts.getTime(0);
    const double t1 = ts.getTime(n - 1);

    const double tmin = std::min(t0, t1);
    const double tmax = std::max(t0, t1);

    return (t >= tmin && t <= tmax);
}

static TimeSeries<double> mean_ts_by_mode(const TimeSeriesSet<double>& set, RunOptions::MeanTSMode mode)
{
    TimeSeries<double> out;

    const int ns = (int)set.size();
    if (ns <= 0) return out;

    auto pick_ref_index_first = [&]() -> int { return 0; };

    auto pick_ref_index_longest = [&]() -> int {
        int best = 0;
        int best_n = -1;
        for (int i = 0; i < ns; ++i) {
            const int n = (int)set[i].size();
            if (n > best_n) { best_n = n; best = i; }
        }
        return best;
    };

    if (mode == RunOptions::MeanTSMode::Union) {

        std::vector<double> grid;
        grid.reserve(4096);

        for (int i = 0; i < ns; ++i) {
            const int n = (int)set[i].size();
            for (int k = 0; k < n; ++k) {
                const double t = set[i].getTime(k);
                if (std::isfinite(t)) grid.push_back(t);
            }
        }

        if (grid.size() < 2) return out;

        std::sort(grid.begin(), grid.end());

        // unique with tolerance
        const double eps = 1e-12;
        std::vector<double> uniq;
        uniq.reserve(grid.size());
        for (size_t i = 0; i < grid.size(); ++i) {
            if (uniq.empty() || std::abs(grid[i] - uniq.back()) > eps) uniq.push_back(grid[i]);
        }

        for (size_t it = 0; it < uniq.size(); ++it) {
            const double t = uniq[it];

            double sum = 0.0;
            int cnt = 0;

            for (int j = 0; j < ns; ++j) {
                const auto& ts = set[j];
                if (!ts_has_support_at_t(ts, t)) continue;

                const double v = ts.interpol(t);
                if (std::isfinite(v)) { sum += v; cnt++; }
            }

            if (cnt > 0) out.append(t, sum / (double)cnt);
        }

        return out;
    }

    // Reference-grid (First/Longest)
    const int ref =
        (mode == RunOptions::MeanTSMode::First) ? pick_ref_index_first()
                                               : pick_ref_index_longest();

    const auto& ref_ts = set[ref];
    const int nref = (int)ref_ts.size();
    if (nref < 2) return out;

    for (int k = 0; k < nref; ++k) {
        const double t = ref_ts.getTime(k);
        if (!std::isfinite(t)) continue;

        double sum = 0.0;
        int cnt = 0;

        for (int j = 0; j < ns; ++j) {
            const auto& ts = set[j];
            if (!ts_has_support_at_t(ts, t)) continue;

            const double v = ts.interpol(t);
            if (std::isfinite(v)) { sum += v; cnt++; }
        }

        if (cnt > 0) out.append(t, sum / (double)cnt);
    }

    return out;
}

static inline TimeSeries<double> mean_ts_by_opts(const TimeSeriesSet<double>& set, const RunOptions& opts)
{
    return mean_ts_by_mode(set, opts.mean_ts_mode);
}

// ============================================================================
// Robust realization reject helpers (delete folder + log)
// ============================================================================

static bool remove_dir_recursive(const std::string& p)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    if (p.empty()) return true;
    if (!fs::exists(fs::path(p), ec)) return true;
    fs::remove_all(fs::path(p), ec);
    return !ec;
}

static void append_reject_log(const std::string& run_dir,
                              int realization,
                              int attempt,
                              const std::string& stage,
                              const std::string& what)
{
    const std::string path = joinPath(run_dir, "rejected_realizations.csv");
    const bool exists = fileExists(path);

    std::ofstream f(path, std::ios::app);
    if (!f) return;

    if (!exists) f << "realization,attempt,stage,message\n";

    auto esc = [](std::string s){
        std::replace(s.begin(), s.end(), '\n', ' ');
        return s;
    };

    f << realization << "," << attempt << "," << esc(stage) << "," << esc(what) << "\n";
}

// ============================================================================
// PT mean-arrival helpers
// ============================================================================

static inline double trapz_integral(const TimeSeries<double>& ts)
{
    const int n = (int)ts.size();
    if (n < 2) return 0.0;

    double area = 0.0;
    for (int i = 1; i < n; ++i) {
        const double t0 = ts.getTime(i - 1);
        const double t1 = ts.getTime(i);
        const double y0 = ts.getValue(i - 1);
        const double y1 = ts.getValue(i);
        const double dt = (t1 - t0);
        if (!std::isfinite(dt) || dt <= 0.0) continue;
        if (!std::isfinite(y0) || !std::isfinite(y1)) continue;
        area += 0.5 * (y0 + y1) * dt;
    }
    return area;
}

static inline double trapz_integral_t_times_y(const TimeSeries<double>& ts)
{
    const int n = (int)ts.size();
    if (n < 2) return 0.0;

    double area = 0.0;
    for (int i = 1; i < n; ++i) {
        const double t0 = ts.getTime(i - 1);
        const double t1 = ts.getTime(i);
        const double y0 = ts.getValue(i - 1);
        const double y1 = ts.getValue(i);
        const double dt = (t1 - t0);
        if (!std::isfinite(dt) || dt <= 0.0) continue;
        if (!std::isfinite(y0) || !std::isfinite(y1)) continue;
        area += 0.5 * (t0 * y0 + t1 * y1) * dt;
    }
    return area;
}

static inline double mean_time_from_pdf(const TimeSeries<double>& pdf_ts)
{
    const double den = trapz_integral(pdf_ts);
    if (!(den > 0.0) || !std::isfinite(den)) return std::numeric_limits<double>::quiet_NaN();
    const double num = trapz_integral_t_times_y(pdf_ts);
    return num / den;
}

static bool read_pt_mean_csv(const std::string& path, std::vector<double>& x, std::vector<double>& mu)
{
    x.clear(); mu.clear();

    std::ifstream f(path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    const char delim = detect_delim(header);

    std::string line;
    while (std::getline(f, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        const auto tok = split_line_delim(line, delim);
        if (tok.size() < 2) continue;

        double xv = 0.0, muv = 0.0;
        if (!try_parse_double(tok[0], xv)) continue;
        if (!try_parse_double(tok[1], muv)) continue;

        x.push_back(xv);
        mu.push_back(muv);
    }

    return !x.empty();
}

// ============================================================================
// qx-rank helpers (velocity-rank copula diagnostics)
// ============================================================================

struct QxRankStats
{
    double delta_x = 0.0;
    int n_pairs = 0;
    double corr_qx = std::numeric_limits<double>::quiet_NaN();
    double corr_rank = std::numeric_limits<double>::quiet_NaN();
    double gaussian_copula_rho = std::numeric_limits<double>::quiet_NaN();
    double kendall_tau = std::numeric_limits<double>::quiet_NaN();
    double rho_from_tau = std::numeric_limits<double>::quiet_NaN();
    double mardia_skewness = std::numeric_limits<double>::quiet_NaN();
    double mardia_kurtosis = std::numeric_limits<double>::quiet_NaN();
    double gaussian_copula_gof_stat = std::numeric_limits<double>::quiet_NaN();
    double gaussian_copula_gof_pvalue = std::numeric_limits<double>::quiet_NaN();
};

static double pearson_corr(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size() || a.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    const int n = (int)a.size();
    double ma = 0.0, mb = 0.0;
    for (int i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
    ma /= (double)n; mb /= (double)n;

    double num = 0.0, va = 0.0, vb = 0.0;
    for (int i = 0; i < n; ++i) {
        const double da = a[i] - ma;
        const double db = b[i] - mb;
        num += da * db;
        va += da * da;
        vb += db * db;
    }
    if (!(va > 0.0) || !(vb > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    return num / std::sqrt(va * vb);
}

static std::vector<double> fractional_ranks_01(const std::vector<double>& x)
{
    const int n = (int)x.size();
    std::vector<double> r(n, std::numeric_limits<double>::quiet_NaN());
    if (n <= 0) return r;

    std::vector<std::pair<double,int>> xv;
    xv.reserve(n);
    for (int i = 0; i < n; ++i) xv.push_back({x[i], i});
    std::sort(xv.begin(), xv.end(), [](const auto& a, const auto& b){ return a.first < b.first; });

    int i = 0;
    while (i < n) {
        int j = i + 1;
        while (j < n && std::abs(xv[j].first - xv[i].first) <= 1e-14) ++j;
        const double rank_avg = 0.5 * (double(i + 1) + double(j)); // 1-based average rank for ties
        for (int k = i; k < j; ++k) {
            // plotting position in (0,1): avoids exact 0 and 1 for normal-score map
            r[xv[k].second] = rank_avg / (double(n) + 1.0);
        }
        i = j;
    }
    return r;
}

// Acklam-style inverse normal CDF approximation.
static double norminv_approx(double p)
{
    const double plow = 0.02425;
    const double phigh = 1.0 - plow;
    if (!(p > 0.0 && p < 1.0)) return std::numeric_limits<double>::quiet_NaN();

    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                                2.445134137142996e+00, 3.754408661907416e+00};

    if (p < plow) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
    if (p > phigh) {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                 ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }

    const double q = p - 0.5;
    const double r = q * q;
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
           (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
}

static double normal_cdf(double z)
{
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

static double gaussian_copula_cdf_approx(double u, double v, double rho)
{
    constexpr double kPi = 3.14159265358979323846;
    if (!(u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0)) return std::numeric_limits<double>::quiet_NaN();
    const double a = norminv_approx(u);
    const double b = norminv_approx(v);
    if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(rho)) return std::numeric_limits<double>::quiet_NaN();
    rho = std::max(-0.999, std::min(0.999, rho));

    const double s = std::sqrt(1.0 - rho * rho);
    const double lo = -8.0;
    if (a <= lo) return 0.0;
    const double hi = std::min(8.0, a);
    const int N = 80;
    const double h = (hi - lo) / double(N);
    auto phi = [kPi](double z) { return std::exp(-0.5 * z * z) / std::sqrt(2.0 * kPi); };

    double sum = 0.0;
    for (int i = 0; i <= N; ++i) {
        const double x = lo + h * i;
        const double w = (i == 0 || i == N) ? 1.0 : ((i % 2 == 0) ? 2.0 : 4.0);
        const double arg = (b - rho * x) / s;
        sum += w * phi(x) * normal_cdf(arg);
    }
    return std::max(0.0, std::min(1.0, sum * h / 3.0));
}

static double kendall_tau_naive(const std::vector<double>& x, const std::vector<double>& y)
{
    if (x.size() != y.size() || x.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    long long concordant = 0;
    long long discordant = 0;
    const int n = (int)x.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double dx = x[j] - x[i];
            const double dy = y[j] - y[i];
            const double s = dx * dy;
            if (s > 0.0) ++concordant;
            else if (s < 0.0) ++discordant;
        }
    }
    const long long den = (long long)n * (n - 1) / 2;
    if (den <= 0) return std::numeric_limits<double>::quiet_NaN();
    return double(concordant - discordant) / double(den);
}

static void mardia_bivariate_moments(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     double& skew_out,
                                     double& kurt_out)
{
    skew_out = std::numeric_limits<double>::quiet_NaN();
    kurt_out = std::numeric_limits<double>::quiet_NaN();
    if (x.size() != y.size() || x.size() < 3) return;
    const int n = (int)x.size();
    double mx = 0.0, my = 0.0;
    for (int i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;

    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (int i = 0; i < n; ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    sxx /= (n - 1); syy /= (n - 1); sxy /= (n - 1);
    const double det = sxx * syy - sxy * sxy;
    if (!(det > 0.0)) return;

    const double inv00 =  syy / det;
    const double inv01 = -sxy / det;
    const double inv11 =  sxx / det;

    std::vector<double> d2(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        d2[i] = dx * (inv00 * dx + inv01 * dy) + dy * (inv01 * dx + inv11 * dy);
    }

    double b2p = 0.0;
    for (double v : d2) b2p += v * v;
    b2p /= n;
    kurt_out = b2p;

    double b1p = 0.0;
    for (int i = 0; i < n; ++i) {
        const double dxi = x[i] - mx;
        const double dyi = y[i] - my;
        for (int j = 0; j < n; ++j) {
            const double dxj = x[j] - mx;
            const double dyj = y[j] - my;
            const double dij = dxi * (inv00 * dxj + inv01 * dyj) + dyi * (inv01 * dxj + inv11 * dyj);
            b1p += dij * dij * dij;
        }
    }
    b1p /= (double(n) * n);
    skew_out = b1p;
}

static double gaussian_copula_cvm_statistic(const std::vector<double>& u,
                                            const std::vector<double>& v,
                                            double rho)
{
    if (u.size() != v.size() || u.empty()) return std::numeric_limits<double>::quiet_NaN();
    const int n = (int)u.size();
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        int count = 0;
        for (int j = 0; j < n; ++j) {
            if (u[j] <= u[i] && v[j] <= v[i]) ++count;
        }
        const double cn = double(count) / double(n);
        const double cg = gaussian_copula_cdf_approx(u[i], v[i], rho);
        if (std::isfinite(cg)) {
            const double d = cn - cg;
            s += d * d;
        }
    }
    return s;
}

static double estimate_gaussian_copula_gof_pvalue(const std::vector<double>& u,
                                                  const std::vector<double>& v,
                                                  double rho_hat,
                                                  int B,
                                                  unsigned long seed,
                                                  double& stat_obs)
{
    constexpr double kPi = 3.14159265358979323846;
    stat_obs = gaussian_copula_cvm_statistic(u, v, rho_hat);
    if (!std::isfinite(stat_obs) || B <= 0) return std::numeric_limits<double>::quiet_NaN();

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> N01(0.0, 1.0);
    int ge = 0;
    const int n = (int)u.size();
    const double rho = std::max(-0.999, std::min(0.999, rho_hat));
    const double s = std::sqrt(1.0 - rho * rho);
    std::vector<double> ub(n), vb(n);
    for (int b = 0; b < B; ++b) {
        for (int i = 0; i < n; ++i) {
            const double z1 = N01(rng);
            const double z2 = rho * z1 + s * N01(rng);
            ub[i] = std::min(1.0 - 1e-12, std::max(1e-12, normal_cdf(z1)));
            vb[i] = std::min(1.0 - 1e-12, std::max(1e-12, normal_cdf(z2)));
        }
        const double tau_b = kendall_tau_naive(ub, vb);
        const double rho_b = std::sin(0.5 * kPi * tau_b);
        const double s_b = gaussian_copula_cvm_statistic(ub, vb, rho_b);
        if (std::isfinite(s_b) && s_b >= stat_obs) ++ge;
    }
    return double(ge + 1) / double(B + 1);
}

static QxRankStats analyze_and_write_rank_pairs(
    const PathwaySet& pair_set,
    double delta_x,
    const std::string& scatter_csv_path)
{
    QxRankStats out;
    out.delta_x = delta_x;

    if (pair_set.size() < 2) return out;
    const auto& p1 = pair_set[0];
    const auto& p2 = pair_set[1];
    const int n = (int)std::min(p1.size(), p2.size());
    if (n < 2) return out;

    std::vector<double> q1(n), q2(n);
    for (int i = 0; i < n; ++i) {
        q1[i] = p1[i].qx();
        q2[i] = p2[i].qx();
    }
    const std::vector<double> u1 = fractional_ranks_01(q1);
    const std::vector<double> u2 = fractional_ranks_01(q2);

    std::vector<double> z1; z1.reserve(n);
    std::vector<double> z2; z2.reserve(n);
    std::ofstream f(scatter_csv_path);
    if (f) {
        f << "u1,u2,qx1,qx2,qx1_normal_score,qx2_normal_score\n";
        f << std::setprecision(15);
    }

    for (int i = 0; i < n; ++i) {
        const double nz1 = norminv_approx(u1[i]);
        const double nz2 = norminv_approx(u2[i]);
        if (f) f << u1[i] << "," << u2[i] << "," << q1[i] << "," << q2[i] << "," << nz1 << "," << nz2 << "\n";
        if (std::isfinite(nz1) && std::isfinite(nz2)) {
            z1.push_back(nz1);
            z2.push_back(nz2);
        }
    }

    out.n_pairs = n;
    out.corr_qx = pearson_corr(q1, q2);
    out.corr_rank = pearson_corr(u1, u2);          // Spearman-like via normalized ranks
    out.gaussian_copula_rho = pearson_corr(z1, z2); // rho of Gaussianized ranks
    out.kendall_tau = kendall_tau_naive(q1, q2);
    out.rho_from_tau = std::sin(0.5 * 3.14159265358979323846 * out.kendall_tau);
    mardia_bivariate_moments(q1, q2, out.mardia_skewness, out.mardia_kurtosis);
    out.gaussian_copula_gof_pvalue =
        estimate_gaussian_copula_gof_pvalue(u1, u2, out.rho_from_tau, 200, 777UL + (unsigned long)(1000.0 * delta_x),
                                            out.gaussian_copula_gof_stat);
    return out;
}

// ============================================================================
// Fine-scale helpers
// ============================================================================

static void generateKField(Grid2D& g, const SimParams& P, unsigned long seed)
{
    if (P.CorrelationModel == SimParams::correlationmode::gaussian) {
        g.generateGaussianRandomField("K_normal_score",
                                      P.correlation_ls_x,
                                      P.correlation_ls_y,
                                      0.0, 1.0, seed);
    }
    else {
        g.makeGaussianFieldSGS("K_normal_score", P.correlation_ls_x, P.correlation_ls_y, 10, seed);
    }
}

// ============================================================================
// Wiener diffusion runner (kept)
// ============================================================================

int runDiffusionSimulation(const RunOptions &opts, int realization, const std::string &output_dir)
{
    std::string m = to_lower_copy(trim_copy(opts.wiener_mode));
    WienerParams wp;
    wp.dt   = opts.wiener_dt;
    wp.D    = opts.wiener_Dx;
    wp.rx   = opts.wiener_rx;
    wp.ry   = opts.wiener_ry;
    wp.seed = opts.wiener_seed + (unsigned long)realization;
    if (m == "1dx")      wp.mode = WienerMode::W1D_X;
    else if (m == "1dy") wp.mode = WienerMode::W1D_Y;
    else                 wp.mode = WienerMode::W2D;

    std::cout << "Running diffusion simulation\n"
              << "  Mode: " << opts.wiener_mode << "\n"
              << "  D = " << wp.D << ", dt = " << wp.dt << "\n"
              << "  Ellipse: rx = " << wp.rx << ", ry = " << wp.ry << "\n"
              << "  Particles: " << realization << "\n";

    PathwaySet pathwaySet;
    pathwaySet.InitializeAtOrigin(realization);

    TimeSeries<double> times = pathwaySet.trackDiffusion(wp.dt, wp.rx, wp.ry, wp.D);
    times.writefile(joinPath(output_dir, "diffusiontimes.csv"));

    double mean     = times.mean();
    double stdev    = times.stddev();
    double meanlog  = times.mean_log();
    double stdevlog = times.log().stddev();
    auto [tmin, tmax] = pathwaySet.travelTimeRange();

    TimeSeries<double> fptDist = times.distribution(100, 0);
    fptDist.writefile(joinPath(output_dir, "FirstPassageTimeDistribution.csv"));

    {
        std::ofstream f(joinPath(output_dir, "diffusion_stats.txt"));
        if (f.is_open()) {
            f << "# First Passage Time Statistics\n"
              << "mode=" << opts.wiener_mode << "\n"
              << "D=" << wp.D << "\n"
              << "dt=" << wp.dt << "\n"
              << "rx=" << wp.rx << "\n"
              << "ry=" << wp.ry << "\n"
              << "seed=" << wp.seed << "\n"
              << "num_particles=" << realization << "\n"
              << "mean=" << mean << "\n"
              << "stdev=" << stdev << "\n"
              << "mean_log=" << meanlog << "\n"
              << "stdev_log=" << stdevlog << "\n"
              << "min=" << tmin << "\n"
              << "max=" << tmax << "\n";
        }
    }

    return 0;
}

// ============================================================================
// Mean correlations + stats helpers
// ============================================================================

static void writeMeanCorrelations(
    const std::string& run_dir,
    const SimParams& P,
    FineScaleOutputs& outputs,
    const RunOptions& opts)
{
    outputs.advective_correlations.write(joinPath(run_dir, "advective_correlations.txt"));
    outputs.velocity_x_correlations.write(joinPath(run_dir, "diffusion_x_correlations.txt"));
    outputs.velocity_y_correlations.write(joinPath(run_dir, "diffusion_y_correlations.txt"));
    outputs.K_x_correlations.write(joinPath(run_dir, "K_x_correlations.txt"));
    outputs.K_y_correlations.write(joinPath(run_dir, "K_y_correlations.txt"));

    // mean_ts() switch (from main via opts.mean_ts_mode)
    auto mean_corr_a    = mean_ts_by_opts(outputs.advective_correlations, opts);
    auto mean_corr_x    = mean_ts_by_opts(outputs.velocity_x_correlations, opts);
    auto mean_corr_y    = mean_ts_by_opts(outputs.velocity_y_correlations, opts);
    auto mean_K_corr_x  = mean_ts_by_opts(outputs.K_x_correlations, opts);
    auto mean_K_corr_y  = mean_ts_by_opts(outputs.K_y_correlations, opts);

    mean_corr_a.writefile(joinPath(run_dir, "advective_correlations_mean.txt"));
    mean_corr_x.writefile(joinPath(run_dir, "diffusion_x_correlations_mean.txt"));
    mean_corr_y.writefile(joinPath(run_dir, "diffusion_y_correlations_mean.txt"));
    mean_K_corr_x.writefile(joinPath(run_dir, "K_x_correlations_mean.txt"));
    mean_K_corr_y.writefile(joinPath(run_dir, "K_y_correlations_mean.txt"));

    double l_a = mean_corr_a.fitExponentialDecay();
    TimeSeries<double> fitted_a;
    for (size_t i = 0; i < mean_corr_a.size(); ++i) {
        double r = mean_corr_a.getTime((int)i);
        fitted_a.append(r, std::exp(-r / l_a));
    }
    fitted_a.writefile(joinPath(run_dir, "advective_correlations_mean_fitted.txt"));

    auto fitAndWrite = [&](const TimeSeries<double>& mean_corr,
                           const std::string& filename,
                           double& nu_out) {
        double ell = 0, nu = 0.0;
        if (P.CorrelationModel == SimParams::correlationmode::exponentialfit)
            ell = mean_corr.fitExponentialDecay();
        else if (P.CorrelationModel == SimParams::correlationmode::derivative)
            ell = mean_corr.getTime(0) / (1.0 - mean_corr.getValue(0));
        else if (P.CorrelationModel == SimParams::correlationmode::gaussian)
            ell = mean_corr.fitGaussianDecay();
        else {
            auto pr = mean_corr.fitMaternDecay();
            nu = pr.first; ell = pr.second;
        }

        TimeSeries<double> fitted;
        for (size_t i = 0; i < mean_corr.size(); ++i) {
            double r = mean_corr.getTime((int)i);
            double rho =
                (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                    std::exp(-(r / ell) * (r / ell)) :
                (P.CorrelationModel == SimParams::correlationmode::matern) ?
                    matern(r, nu, ell) :
                    std::exp(-r / ell);
            fitted.append(r, rho);
        }
        fitted.writefile(joinPath(run_dir, filename));
        nu_out = nu;
        return ell;
    };

    fitAndWrite(mean_corr_x, "diffusion_x_correlations_mean_fitted.txt", outputs.nu_x_mean);
    fitAndWrite(mean_corr_y, "diffusion_y_correlations_mean_fitted.txt", outputs.nu_y_mean);

    auto fitKAndWrite = [&](const TimeSeries<double>& mean_K_corr, const std::string& filename) {
        double ell_K =
            (P.CorrelationModel == SimParams::correlationmode::gaussian)       ? mean_K_corr.fitGaussianDecay() :
            (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? mean_K_corr.fitExponentialDecay() :
            (P.CorrelationModel == SimParams::correlationmode::derivative)      ? mean_K_corr.getTime(0) / (1.0 - mean_K_corr.getValue(0)) :
                                                                                mean_K_corr.fitMaternDecay().second;

        TimeSeries<double> fitted_K;
        for (size_t i = 0; i < mean_K_corr.size(); ++i) {
            double r = mean_K_corr.getTime((int)i);
            double rho =
                (P.CorrelationModel == SimParams::correlationmode::gaussian) ?
                    std::exp(-(r / ell_K) * (r / ell_K)) :
                    std::exp(-r / ell_K);
            fitted_K.append(r, rho);
        }
        fitted_K.writefile(joinPath(run_dir, filename));
        return ell_K;
    };

    outputs.lambda_K_x_mean = fitKAndWrite(mean_K_corr_x, "K_x_correlations_mean_fitted.txt");
    outputs.lambda_K_y_mean = fitKAndWrite(mean_K_corr_y, "K_y_correlations_mean_fitted.txt");
}

static void computeFinalMeans(FineScaleOutputs& outputs)
{
    outputs.lc_mean          = mean_of(outputs.lc_all);
    outputs.velocity_lx_mean = mean_of(outputs.velocity_lx_all);
    outputs.velocity_ly_mean = mean_of(outputs.velocity_ly_all);
    outputs.dt_mean          = mean_of(outputs.dt_all);

    // NOTE:
    // outputs.lambda_K_x_mean / lambda_K_y_mean are set in writeMeanCorrelations().
    // In hardcoded mode they may remain NaN unless you explicitly set them.
}

// ============================================================================
// Fine loop (RETRY/REJECT)
//
// Collects:
//   - transport CDF stacks (raw BTCs)
//   - transport PDF stacks (derivatives of BTCs)
//   - PT PDF/CDF stacks
//
// IMPORTANT transport convention used by the latest code:
//
//   Legacy transport outputs written later in main.cpp:
//       x=...BTC_Compare.csv
//       BTC_mean.csv
//
//   now store TRANSPORT PDF, not raw/CDF.
//
//   New explicit transport files also written later:
//       x=...BTC_Compare_pdf.csv
//       x=...BTC_Compare_cdf.csv
//       BTC_mean_pdf.csv
//       BTC_mean_cdf.csv
//
// Internally, this function therefore stores BOTH transport forms:
//   - BTCs_transport_cdf : raw BTC / cumulative-style transport curves
//   - BTCs_transport_pdf : derivative of those BTC curves
// ============================================================================

static bool run_fine_loop_collect(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    FineScaleOutputs& outputs,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;

    const std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,lambda_K_x,lambda_K_y,dt_opt\n";
    }

    // Ensure stacks exist
    outputs.BTCs_transport_pdf.resize(P.xLocations.size());
    outputs.BTCs_transport_cdf.resize(P.xLocations.size());
    outputs.PT_pdfs.resize(P.xLocations.size());
    outputs.PT_cdfs.resize(P.xLocations.size());

    const int nReal = P.nReal_default;

    // cap to avoid infinite loop if PETSc is fundamentally failing
    const int MAX_RETRY_PER_REAL = 50;

    int success = 0;
    while (success < nReal) {

        const int r = success + 1; // realization number to produce
        int attempt = 0;

        bool done_this_real = false;

        while (!done_this_real) {
            attempt++;

            const std::string rlab = makeRealLabel(r);
            std::string fine_dir = joinPath(run_dir, makeFineFolder(r));

            if (rank == 0) {
                if (dirExists(fine_dir)) remove_dir_recursive(fine_dir);
                createDirectory(fine_dir);
            }
            MPI_Barrier(PETSC_COMM_WORLD);

            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

            const std::string pfx = rlab + "_";

            // ---- locals (do NOT pollute outputs unless we succeed) ----
            TimeSeries<double> K_corr_x, K_corr_y;
            TimeSeries<double> vel_corr_x, vel_corr_y;
            TimeSeries<double> qx_corr_adv;
            TimeSeries<double> qx_corr_adv_rank_copula;

            double lambda_K_x_emp = std::numeric_limits<double>::quiet_NaN();
            double lambda_K_y_emp = std::numeric_limits<double>::quiet_NaN();
            double lambda_x_emp   = std::numeric_limits<double>::quiet_NaN();
            double lambda_y_emp   = std::numeric_limits<double>::quiet_NaN();
            double lc_emp         = std::numeric_limits<double>::quiet_NaN();
            double lc_emp_raw_qx  = std::numeric_limits<double>::quiet_NaN();
            double lc_emp_rank_copula = std::numeric_limits<double>::quiet_NaN();
            double dt_optimal     = std::numeric_limits<double>::quiet_NaN();

            TimeSeries<double> qx_inverse_cdf;
            TimeSeries<double> qx_pdf;

            TimeSeriesSet<double> BTCs_FineScaled;
            TimeSeriesSet<double> PT_pdf_local;
            TimeSeriesSet<double> PT_cdf_local;

            try {
                Grid2D g(nx, ny, Lx, Ly);

                // IMPORTANT: change seed on retries so you don't repeat the same failing field
                unsigned long seed =
                    P.run_seed
                    + 1000UL * (unsigned long)r
                    + 1000000UL * (unsigned long)attempt
                    + (unsigned long)rank;

                generateKField(g, P, seed);

                // ---- K correlation sampling (local) ----
                {
                    const int num_deltas = 30;
                    const int num_samples_per_delta = 10000;

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_x_range.first *
                                       std::pow(P.correlation_x_range.second / P.correlation_x_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "K_normal_score", Grid2D::ArrayKind::Cell,
                                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                            K_corr_x.append(delta, samples.correlation_tc());
                        } catch (...) {}
                    }
                    K_corr_x.writefile(joinPath(fine_dir, pfx + "K_correlation_x.txt"));

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_y_range.first *
                                       std::pow(P.correlation_y_range.second / P.correlation_y_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "K_normal_score", Grid2D::ArrayKind::Cell,
                                num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                            K_corr_y.append(delta, samples.correlation_tc());
                        } catch (...) {}
                    }
                    K_corr_y.writefile(joinPath(fine_dir, pfx + "K_correlation_y.txt"));

                    lambda_K_x_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? K_corr_x.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? K_corr_x.getTime(0)/(1.0-K_corr_x.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? K_corr_x.fitGaussianDecay() :
                                                                                             K_corr_x.fitMaternDecay().second;

                    lambda_K_y_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? K_corr_y.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? K_corr_y.getTime(0)/(1.0-K_corr_y.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? K_corr_y.fitGaussianDecay() :
                                                                                             K_corr_y.fitMaternDecay().second;
                }

                g.normalizeField("K_normal_score", Grid2D::ArrayKind::Cell);
                g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
                g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.vti"));

                g.createExponentialField("K_normal_score", P.stdev, P.g_mean, "K");
                g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));

                // ---- PETSc solve: can throw (KSP did not converge) ----
                g.DarcySolve(Lx, 0, "K", "K");

                g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.vti"));
                g.writeNamedVTI("qx",   Grid2D::ArrayKind::Fx,   joinPath(fine_dir, pfx + "qx.vti"));
                g.writeNamedVTI("qy",   Grid2D::ArrayKind::Fy,   joinPath(fine_dir, pfx + "qy.vti"));

                // ---- velocity normal score + correlations (local) ----
                {
                    TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
                    TimeSeries<double> QxRanks = AllQxValues.ConvertToRanks();
                    TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
                    g.assignFromTimeSeries(QxRanks, "qx_ranks", Grid2D::ArrayKind::Fx);
                    g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);
                    g.writeNamedVTI("qx_ranks", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx_ranks.vti"));
                    g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx_normal_score.vti"));

                    const int num_deltas = 30;
                    const int num_samples_per_delta = 10000;
                    TimeSeries<double> vel_rank_corr_x;
                    TimeSeries<double> vel_rank_corr_y;

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_x_range.first *
                                       std::pow(P.correlation_x_range.second / P.correlation_x_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "qx_normal_score", Grid2D::ArrayKind::Fx,
                                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                            vel_corr_x.append(delta, samples.correlation_tc());
                            TimeSeries<double> rank_samples = g.sampleGaussianPerturbation(
                                "qx_ranks", Grid2D::ArrayKind::Fx,
                                num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                            vel_rank_corr_x.append(delta, rank_samples.correlation_tc());
                        } catch (...) {}
                    }
                    vel_corr_x.writefile(joinPath(fine_dir, pfx + "velocity_correlation_x.txt"));
                    vel_rank_corr_x.writefile(joinPath(fine_dir, pfx + "velocity_rank_correlation_x.txt"));

                    for (int i = 0; i < num_deltas; ++i) {
                        double exponent = static_cast<double>(i) / (num_deltas - 1);
                        double delta = P.correlation_y_range.first *
                                       std::pow(P.correlation_y_range.second / P.correlation_y_range.first, exponent);
                        try {
                            TimeSeries<double> samples = g.sampleGaussianPerturbation(
                                "qx_normal_score", Grid2D::ArrayKind::Fx,
                                num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                            vel_corr_y.append(delta, samples.correlation_tc());
                            TimeSeries<double> rank_samples = g.sampleGaussianPerturbation(
                                "qx_ranks", Grid2D::ArrayKind::Fx,
                                num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                            vel_rank_corr_y.append(delta, rank_samples.correlation_tc());
                        } catch (...) {}
                    }
                    vel_corr_y.writefile(joinPath(fine_dir, pfx + "velocity_correlation_y.txt"));
                    vel_rank_corr_y.writefile(joinPath(fine_dir, pfx + "velocity_rank_correlation_y.txt"));

                    lambda_x_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? vel_corr_x.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? vel_corr_x.getTime(0)/(1.0-vel_corr_x.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? vel_corr_x.fitGaussianDecay() :
                                                                                             vel_corr_x.fitMaternDecay().second;

                    lambda_y_emp =
                        (P.CorrelationModel == SimParams::correlationmode::exponentialfit) ? vel_corr_y.fitExponentialDecay() :
                        (P.CorrelationModel == SimParams::correlationmode::derivative)      ? vel_corr_y.getTime(0)/(1.0-vel_corr_y.getValue(0)) :
                        (P.CorrelationModel == SimParams::correlationmode::gaussian)        ? vel_corr_y.fitGaussianDecay() :
                                                                                             vel_corr_y.fitMaternDecay().second;
                }

                qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100, 1e-6);
                qx_pdf         = g.extractFieldPDF("qx", Grid2D::ArrayKind::Fx, 50,  1e-6, true);

                qx_inverse_cdf = qx_inverse_cdf.make_uniform(du_from_nu(P.nu));
                qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));
                qx_pdf.writefile(joinPath(fine_dir, pfx + "qx_pdf.txt"));

                dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;

                // ---- Transport (can also throw; reject if it does) ----
                //
                // BTCs_FineScaled is the RAW transport BTC returned by SolveTransport.
                // In the current output convention:
                //   - raw BTC  -> stored into transport CDF stack
                //   - derivative -> stored into transport PDF stack
                //
                // Legacy transport files written later use the PDF stack.
                g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
                g.SetVal("diffusion", P.Diffusion_coefficient);
                g.SetVal("porosity", 1.0);
                g.SetVal("c_left", 1.0);

                g.setBTCLocations(P.xLocations);

                if (opts.solve_fine_scale_transport) {
                    g.SolveTransport(P.t_end_pdf, std::min(dt_optimal, 0.5 / 10.0),
                                     "transport_", 500, fine_dir, "C", &BTCs_FineScaled, r);

                    // Per-realization diagnostic files:
                    //   BTC_FineScaled.csv            -> raw/CDF-style transport BTC
                    //   BTC_FineScaled_derivative.csv -> PDF-style transport BTC
                    BTCs_FineScaled.write(joinPath(fine_dir, pfx + "BTC_FineScaled.csv"));
                    BTCs_FineScaled.derivative().write(joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv"));
                }

                // ---- Advective correlation + PT locals ----
                const bool need_pathway_tracking =
                    opts.perform_particle_tracking || opts.analyze_qx_ranks;

                if (need_pathway_tracking) {
                    PathwaySet pathways;
                    pathways.Initialize(10000, PathwaySet::Weighting::FluxWeighted, &g);
                    pathways.trackAllPathwaysWithDiffusion(&g, 0.01, g.getDiffusion());

                    const std::vector<double>& btc_locations = g.getBTCLocations();
                    if (!btc_locations.empty()) {
                        TimeSeriesSet<double> btc_pdf = pathways.getBreakthroughCurve(
                            btc_locations, true, PathwaySet::BTCType::PDF, 200);
                        btc_pdf.write(joinPath(fine_dir, pfx + "btc_pdf.csv"));

                        TimeSeriesSet<double> btc_cdf = pathways.getBreakthroughCurve(
                            btc_locations, true, PathwaySet::BTCType::CDF);
                        btc_cdf.write(joinPath(fine_dir, pfx + "btc_cdf.csv"));

                        if (opts.perform_particle_tracking) {
                            PT_pdf_local = btc_pdf;
                            PT_cdf_local = btc_cdf;

                            const std::string mean_csv = joinPath(fine_dir, pfx + "PT_mean_arrival_time.csv");
                            std::ofstream mf(mean_csv);
                            mf << "x,mean_arrival_time\n";
                            mf << std::setprecision(15);
                            for (int i = 0; i < (int)btc_locations.size() && i < (int)btc_pdf.size(); ++i) {
                                const double mu = mean_time_from_pdf(btc_pdf[i]);
                                mf << btc_locations[i] << "," << mu << "\n";
                            }
                        }
                    }

                    double Delta_x_min = 0.001, Delta_x_max = 0.5;
                    int num_Delta_x = 30;
                    int num_samples_per_Delta_x = 10000;
                    if (opts.analyze_qx_ranks) {
                        Delta_x_min = opts.qx_rank_delta_x_min;
                        Delta_x_max = opts.qx_rank_delta_x_max;
                        num_Delta_x = std::max(2, opts.qx_rank_num_delta);
                        num_samples_per_Delta_x = std::max(100, opts.qx_rank_num_samples);
                    }

                    std::vector<QxRankStats> rank_stats;
                    rank_stats.reserve((size_t)num_Delta_x);

                    for (int i = 0; i < num_Delta_x; ++i) {
                        double exponent = static_cast<double>(i) / (num_Delta_x - 1);
                        double Delta_x = Delta_x_min * std::pow(Delta_x_max / Delta_x_min, exponent);
                        try {
                            PathwaySet particle_pairs = pathways.sampleParticlePairs(Delta_x, num_samples_per_Delta_x);
                            double correlation = particle_pairs.calculateCorrelation(0, 1, "qx");
                            qx_corr_adv.append(Delta_x, correlation);

                            if (opts.analyze_qx_ranks) {
                                std::ostringstream fn;
                                fn << pfx << "qx_rank_pairs_dx_" << std::fixed << std::setprecision(6) << Delta_x << ".csv";
                                const QxRankStats st =
                                    analyze_and_write_rank_pairs(particle_pairs, Delta_x, joinPath(fine_dir, fn.str()));
                                rank_stats.push_back(st);
                                qx_corr_adv_rank_copula.append(Delta_x, st.gaussian_copula_rho);
                            }
                        } catch (...) {}
                    }

                    qx_corr_adv.writefile(joinPath(fine_dir, pfx + "qx_correlation_vs_distance.txt"));
                    lc_emp = qx_corr_adv.fitExponentialDecay();
                    lc_emp_raw_qx = lc_emp;

                    if (opts.analyze_qx_ranks) {
                        qx_corr_adv_rank_copula.writefile(
                            joinPath(fine_dir, pfx + "qx_rank_gaussian_copula_correlation_vs_distance.txt"));
                        lc_emp_rank_copula = qx_corr_adv_rank_copula.fitExponentialDecay();
                        if (std::isfinite(lc_emp_rank_copula)) {
                            lc_emp = lc_emp_rank_copula;
                        }
                        if (opts.lc_source == RunOptions::LcSource::RawQx) {
                            lc_emp = lc_emp_raw_qx;
                        } else if (opts.lc_source == RunOptions::LcSource::RankCopula) {
                            if (std::isfinite(lc_emp_rank_copula)) lc_emp = lc_emp_rank_copula;
                            else lc_emp = lc_emp_raw_qx;
                        } else {
                            if (std::isfinite(lc_emp_rank_copula)) lc_emp = lc_emp_rank_copula;
                            else lc_emp = lc_emp_raw_qx;
                        }

                        std::ofstream sf(joinPath(fine_dir, pfx + "qx_rank_copula_summary.csv"));
                        sf << "delta_x,n_pairs,corr_qx,corr_rank,gaussian_copula_rho,kendall_tau,rho_from_tau,mardia_skewness,mardia_kurtosis,gaussian_copula_gof_stat,gaussian_copula_gof_pvalue\n";
                        sf << std::setprecision(15);
                        for (const auto& st : rank_stats) {
                            sf << st.delta_x << ","
                               << st.n_pairs << ","
                               << st.corr_qx << ","
                               << st.corr_rank << ","
                               << st.gaussian_copula_rho << ","
                               << st.kendall_tau << ","
                               << st.rho_from_tau << ","
                               << st.mardia_skewness << ","
                               << st.mardia_kurtosis << ","
                               << st.gaussian_copula_gof_stat << ","
                               << st.gaussian_copula_gof_pvalue << "\n";
                        }

                        std::ofstream lcf(joinPath(fine_dir, pfx + "velocity_correlation_lengths.csv"));
                        lcf << "lc_raw_qx,lc_rank_copula,lc_selected\n";
                        lcf << std::setprecision(15)
                            << lc_emp_raw_qx << ","
                            << lc_emp_rank_copula << ","
                            << lc_emp << "\n";

                        std::ofstream lcsf(joinPath(fine_dir, pfx + "lc_selection_mode.txt"));
                        lcsf << "lc_source=";
                        if (opts.lc_source == RunOptions::LcSource::RawQx) lcsf << "raw_qx\n";
                        else if (opts.lc_source == RunOptions::LcSource::RankCopula) lcsf << "rank_copula\n";
                        else lcsf << "auto_prefer_rank_copula\n";
                    }

                    pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
                    pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"), 1000);
                } else {
                    // No pathway tracing requested by options.
                    // Keep downstream statistics valid by reusing velocity
                    // correlation and its fitted decay as fallback proxies.
                    qx_corr_adv = vel_corr_x;
                    lc_emp = lambda_x_emp;
                }

                // =========================
                // SUCCESS -> COMMIT to outputs
                // =========================
                outputs.K_x_correlations.append(K_corr_x);
                outputs.K_y_correlations.append(K_corr_y);

                outputs.velocity_x_correlations.append(vel_corr_x, "Realization" + aquiutils::numbertostring(r));
                outputs.velocity_y_correlations.append(vel_corr_y, "Realization" + aquiutils::numbertostring(r));
                outputs.advective_correlations.append(qx_corr_adv, "Realization" + aquiutils::numbertostring(r));

                outputs.velocity_lx_all.push_back(lambda_x_emp);
                outputs.velocity_ly_all.push_back(lambda_y_emp);
                outputs.K_lx_all.push_back(lambda_K_x_emp);
                outputs.K_ly_all.push_back(lambda_K_y_emp);
                outputs.dt_all.push_back(dt_optimal);
                outputs.lc_all.push_back(lc_emp);

                if (opts.solve_fine_scale_transport) {
                    TimeSeriesSet<double> BTCs_FineScaled_deriv = BTCs_FineScaled.derivative();

                    // Store RAW transport BTC into explicit CDF stack
                    for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_FineScaled.size(); ++i) {
                        BTCs_FineScaled.setSeriesName(i, fmt_x(P.xLocations[i]));
                        outputs.BTCs_transport_cdf[i].append(BTCs_FineScaled[i], pfx + fmt_x(P.xLocations[i]));
                    }

                    // Store derivative transport BTC into explicit PDF stack
                    for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_FineScaled_deriv.size(); ++i) {
                        BTCs_FineScaled_deriv.setSeriesName(i, fmt_x(P.xLocations[i]));
                        outputs.BTCs_transport_pdf[i].append(BTCs_FineScaled_deriv[i], pfx + fmt_x(P.xLocations[i]));
                    }
                }

                if (opts.perform_particle_tracking) {
                    const int nmin = std::min((int)P.xLocations.size(), (int)PT_pdf_local.size());
                    for (int i = 0; i < nmin; ++i)
                        outputs.PT_pdfs[i].append(PT_pdf_local[i], pfx + fmt_x(P.xLocations[i]));

                    const int nmin2 = std::min((int)P.xLocations.size(), (int)PT_cdf_local.size());
                    for (int i = 0; i < nmin2; ++i)
                        outputs.PT_cdfs[i].append(PT_cdf_local[i], pfx + fmt_x(P.xLocations[i]));
                }

                if (rank == 0) {
                    outputs.inverse_qx_cdfs.append(qx_inverse_cdf, "qx_inverse_cdf" + aquiutils::numbertostring(r));
                    outputs.qx_pdfs.append(qx_pdf, "qx_pdf" + aquiutils::numbertostring(r));

                    std::ofstream f(stats_csv, std::ios::app);
                    f << r << "," << lc_emp << "," << lambda_x_emp << "," << lambda_y_emp << ","
                      << lambda_K_x_emp << "," << lambda_K_y_emp << "," << dt_optimal << "\n";
                }

                done_this_real = true;

            } catch (const std::exception& e) {

                if (rank == 0) {
                    std::cerr << "Rejected realization " << r
                              << " (attempt " << attempt << "): " << e.what() << "\n";
                    append_reject_log(run_dir, r, attempt, "exception", e.what());
                    remove_dir_recursive(fine_dir);
                }

                MPI_Barrier(PETSC_COMM_WORLD);

                if (attempt >= MAX_RETRY_PER_REAL) {
                    if (rank == 0) {
                        std::cerr << "ERROR: realization " << r << " failed "
                                  << attempt << " times. Aborting.\n";
                    }
                    MPI_Abort(PETSC_COMM_WORLD, 123);
                    return false;
                }

                continue;

            } catch (...) {

                if (rank == 0) {
                    std::cerr << "Rejected realization " << r
                              << " (attempt " << attempt << "): unknown exception\n";
                    append_reject_log(run_dir, r, attempt, "unknown", "unknown exception");
                    remove_dir_recursive(fine_dir);
                }

                MPI_Barrier(PETSC_COMM_WORLD);

                if (attempt >= MAX_RETRY_PER_REAL) {
                    MPI_Abort(PETSC_COMM_WORLD, 124);
                    return false;
                }

                continue;
            }

            MPI_Barrier(PETSC_COMM_WORLD);
        }

        // only count after success
        success++;
    }

    // Post-processing
    if (rank == 0) {
        outputs.inverse_qx_cdfs.write(joinPath(run_dir, "qx_inverse_cdfs.txt"));
        outputs.qx_pdfs.write(joinPath(run_dir, "qx_pdfs.txt"));

        // mean_ts switch (from main via opts.mean_ts_mode)
        mean_ts_by_opts(outputs.qx_pdfs, opts).writefile(joinPath(run_dir, "qx_mean_pdf.txt"));

        writeMeanCorrelations(run_dir, P, outputs, opts);
        computeFinalMeans(outputs);

        // mean PT arrival time over realizations
        if (opts.perform_particle_tracking) {
            std::vector<double> sum_mu(P.xLocations.size(), 0.0);
            std::vector<int>    cnt_mu(P.xLocations.size(), 0);

            for (int r = 1; r <= nReal; ++r) {
                const std::string rlab = makeRealLabel(r);
                const std::string pfx  = rlab + "_";
                const std::string fine_dir = joinPath(run_dir, makeFineFolder(r));

                const std::string mean_csv = joinPath(fine_dir, pfx + "PT_mean_arrival_time.csv");
                if (!fileExists(mean_csv)) continue;

                std::vector<double> x, mu;
                if (!read_pt_mean_csv(mean_csv, x, mu)) continue;

                const int nmin = std::min((int)P.xLocations.size(), (int)mu.size());
                for (int i = 0; i < nmin; ++i) {
                    if (std::isfinite(mu[i])) {
                        sum_mu[i] += mu[i];
                        cnt_mu[i] += 1;
                    }
                }
            }

            const std::string out_mean_csv = joinPath(run_dir, "PT_mean_arrival_time_mean.csv");
            std::ofstream outm(out_mean_csv);
            outm << "x,mean_arrival_time\n";
            outm << std::setprecision(15);

            for (int i = 0; i < (int)P.xLocations.size(); ++i) {
                const double m = (cnt_mu[i] > 0) ? (sum_mu[i] / (double)cnt_mu[i])
                                                 : std::numeric_limits<double>::quiet_NaN();
                outm << P.xLocations[i] << "," << m << "\n";
            }
        }
    }

    return true;
}

// ============================================================================
// Mean-building helpers
// ============================================================================

static bool loadHardcodedInverseCDF(
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_cdf_copy_path,
    TimeSeries<double>& invcdf_mean)
{
    bool loaded = false;

    loaded = try_load_hardcoded_invcdf(invcdf_mean, opts.hardcoded_qx_cdf_path);

    if (!loaded && !opts.hardcoded_qx_cdf_path.empty()) {
        const std::string src_dir = dirname_of(opts.hardcoded_qx_cdf_path);
        if (!src_dir.empty()) {
            const std::string meanf = joinPath(src_dir, "mean_qx_inverse_cdf.txt");
            const std::string multi = joinPath(src_dir, "qx_inverse_cdfs.txt");

            if (!fileExists(meanf) && fileExists(multi)) {
                if (build_mean_qx_inverse_cdf_from_multi(multi, meanf)) {
                    loaded = try_load_hardcoded_invcdf(invcdf_mean, meanf);
                }
            }
            if (!loaded && fileExists(meanf)) {
                loaded = try_load_hardcoded_invcdf(invcdf_mean, meanf);
            }
        }
    }

    if (!loaded) {
        build_flat_invcdf(invcdf_mean, H.qx_const);
    }

    invcdf_mean.writefile(run_cdf_copy_path);
    return loaded;
}

static void writeHardcodedMeanParams(
    const std::string& run_dir,
    const FineScaleOutputs& outputs,
    const RunOptions& opts,
    const std::string& run_cdf_copy_path,
    bool cdf_loaded,
    double qx_const)
{
    std::ofstream f(joinPath(run_dir, "mean_params_used.txt"));
    f << "lc_mean=" << outputs.lc_mean << "\n";
    f << "lambda_x_mean=" << outputs.velocity_lx_mean << "\n";
    f << "lambda_y_mean=" << outputs.velocity_ly_mean << "\n";
    f << "lambda_K_x_mean=" << outputs.lambda_K_x_mean << "\n";
    f << "lambda_K_y_mean=" << outputs.lambda_K_y_mean << "\n";
    f << "dt_mean=" << outputs.dt_mean << "\n";
    f << "nu_x_mean=" << outputs.nu_x_mean << "\n";
    f << "nu_y_mean=" << outputs.nu_y_mean << "\n";
    f << "mean_qx_inverse_cdf_written_to=" << run_cdf_copy_path << "\n";
    f << "qx_cdf_preferred_path=" << opts.hardcoded_qx_cdf_path << "\n";
    if (!cdf_loaded) f << "qx_const=" << qx_const << "\n";
}

static void writeComputedMeanParams(
    const std::string& run_dir,
    const SimParams& P,
    const FineScaleOutputs& outputs)
{
    std::ofstream f(joinPath(run_dir, "mean_params.txt"));
    f << "nReal=" << (int)outputs.lc_all.size() << "\n";
    f << "lc_mean=" << outputs.lc_mean << "\n";
    f << "lambda_x_mean=" << outputs.velocity_lx_mean << "\n";
    f << "lambda_y_mean=" << outputs.velocity_ly_mean << "\n";
    f << "lambda_K_x_mean=" << outputs.lambda_K_x_mean << "\n";
    f << "lambda_K_y_mean=" << outputs.lambda_K_y_mean << "\n";
    f << "dt_mean=" << outputs.dt_mean << "\n";
    f << "matern_nu_x_mean=" << outputs.nu_x_mean << "\n";
    f << "matern_nu_y_mean=" << outputs.nu_y_mean << "\n";
    f << "du=" << du_from_nu(P.nu) << "\n";
}

static bool loadMeanParamsForUpscaleOnly_FromInputSource(
    const std::string& input_dir,
    FineScaleOutputs& outputs,
    TimeSeries<double>& invcdf_mean)
{
    const std::string mean_params_path = joinPath(input_dir, "mean_params.txt");
    const std::string mean_cdf_path    = joinPath(input_dir, "mean_qx_inverse_cdf.txt");

    bool ok_params = fileExists(mean_params_path) &&
                     read_mean_params_txt(mean_params_path,
                                          outputs.lc_mean,
                                          outputs.velocity_lx_mean,
                                          outputs.velocity_ly_mean,
                                          outputs.dt_mean,
                                          outputs.nu_x_mean,
                                          outputs.nu_y_mean);

    TimeSeries<double> mean_qx_cdf;
    bool ok_cdf = fileExists(mean_cdf_path) &&
                  read_inverse_cdf_any_format(mean_cdf_path, mean_qx_cdf);

    if (!ok_params || !ok_cdf) {
        std::cerr << "ERROR: --upscale-only requires existing mean_params.txt and mean_qx_inverse_cdf.txt in INPUT folder:\n"
                  << "  " << input_dir << "\n"
                  << "If you don't have them, run with --hardcoded-mean.\n";
        return false;
    }

    invcdf_mean = mean_qx_cdf;
    return true;
}

static bool build_mean_for_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_dir,
    const std::string& resume_run_dir,
    FineScaleOutputs& fine_outputs,
    TimeSeries<double>& invcdf_mean,
    int rank)
{
    if (rank != 0) return true; // rank0 builds; broadcast later

    const std::string run_cdf_copy_path = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

    // Ensure these are always defined (even in hardcoded mode)
    fine_outputs.lambda_K_x_mean = std::numeric_limits<double>::quiet_NaN();
    fine_outputs.lambda_K_y_mean = std::numeric_limits<double>::quiet_NaN();

    // Case 1: hardcoded mean
    if (opts.hardcoded_mean) {
        fine_outputs.lc_mean          = H.lc_mean;
        fine_outputs.velocity_lx_mean = H.lx_mean;
        fine_outputs.velocity_ly_mean = H.ly_mean;
        fine_outputs.dt_mean          = H.dt_mean;
        fine_outputs.nu_x_mean        = H.nu_x;
        fine_outputs.nu_y_mean        = H.nu_y;

        bool loaded = loadHardcodedInverseCDF(opts, H, run_cdf_copy_path, invcdf_mean);
        writeHardcodedMeanParams(run_dir, fine_outputs, opts, run_cdf_copy_path, loaded, H.qx_const);
        return true;
    }

    // Case 2: computed from fine loop
    if (!opts.upscale_only) {
        // mean_ts switch (from main via opts.mean_ts_mode)
        invcdf_mean = mean_ts_by_opts(fine_outputs.inverse_qx_cdfs, opts);
        invcdf_mean.writefile(run_cdf_copy_path);
        writeComputedMeanParams(run_dir, P, fine_outputs);
        return true;
    }

    // Case 3: strict upscale-only load from INPUT SOURCE (resume_run_dir)
    if (resume_run_dir.empty()) {
        std::cerr << "ERROR: strict --upscale-only requires resume_run_dir as INPUT source.\n";
        return false;
    }

    if (!loadMeanParamsForUpscaleOnly_FromInputSource(resume_run_dir, fine_outputs, invcdf_mean)) {
        MPI_Abort(PETSC_COMM_WORLD, 99);
        return false;
    }

    // copy mean CDF into OUTPUT run folder
    invcdf_mean.writefile(run_cdf_copy_path);

    // provenance
    {
        std::ofstream f(joinPath(run_dir, "input_source.txt"));
        f << "input_source_run_dir=" << resume_run_dir << "\n";
    }

    return true;
}

static void computeMeanBTCs(
    const std::vector<double>& xLocations,
    const std::vector<TimeSeriesSet<double>>& stacks,
    TimeSeriesSet<double>& mean_out,
    const RunOptions& opts)
{
    mean_out.clear();
    for (int i = 0; i < (int)xLocations.size(); ++i) {
        const auto m = mean_ts_by_opts(stacks[i], opts);
        mean_out.append(m, fmt_x(xLocations[i]));
    }
}

// -----------------------------------------------------------------------------
// IMPORTANT FIX:
// TimeSeries<double> is NOT guaranteed to be a contiguous double buffer.
// So we broadcast as (times[], values[]) and reconstruct.
// -----------------------------------------------------------------------------
static void broadcastTimeSeries(TimeSeries<double>& ts, int rank)
{
    int n = 0;
    if (rank == 0) n = (int)ts.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    std::vector<double> tbuf, vbuf;
    if (rank == 0) {
        tbuf.resize(n);
        vbuf.resize(n);
        for (int i = 0; i < n; ++i) {
            tbuf[i] = ts.getTime(i);
            vbuf[i] = ts.getValue(i);
        }
    } else {
        tbuf.resize(n);
        vbuf.resize(n);
    }

    if (n > 0) {
        MPI_Bcast(tbuf.data(), n, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
        MPI_Bcast(vbuf.data(), n, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    }

    if (rank != 0) {
        ts.clear();
        for (int i = 0; i < n; ++i) ts.append(tbuf[i], vbuf[i]);
    }
}

static void broadcastMeanParameters(
    FineScaleOutputs& outputs,
    TimeSeries<double>& invcdf_mean,
    int rank)
{
    MPI_Bcast(&outputs.lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.velocity_lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.velocity_ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.nu_x_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.nu_y_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.lambda_K_x_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&outputs.lambda_K_y_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    broadcastTimeSeries(invcdf_mean, rank);
}

// ============================================================================
// Upscaled run (transport + PT separate)
//
// Transport convention used here:
//
//   Raw upscaled BTC   -> explicit transport CDF stack
//   Derivative of BTC  -> explicit transport PDF stack
//
// Later in main.cpp:
//   - legacy transport files use the PDF stack
//   - explicit _pdf/_cdf files use the corresponding explicit stacks
// ============================================================================

static bool run_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    double lc_mean,
    double lx_mean,
    double ly_mean,
    double nu_x_mean,
    double nu_y_mean,
    double dt_mean,
    const TimeSeries<double>& invcdf_mean,
    std::string& up_dir_out,
    std::string& up_btc_path_out,
    std::string& up_btc_deriv_path_out,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs_transport_pdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs_transport_cdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_pdf,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_PT_cdf,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;

    std::string up_dir = joinPath(run_dir, "upscaled_mean");
    if (rank == 0) createDirectory(up_dir);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (!up_dir.empty() && up_dir.back() != '/' && up_dir.back() != '\\') up_dir += "/";

    const std::string up_pfx = "upscaled_";

    Grid2D g_u(nx, ny, Lx, Ly);

    switch (opts.upscaled_mixing_model) {
    case RunOptions::UpscaledMixingModel::Exponential:
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::exponential;
        break;
    case RunOptions::UpscaledMixingModel::ExponentialVDep:
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::exponential_vdep;
        break;
    case RunOptions::UpscaledMixingModel::Gaussian:
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::gaussian;
        break;
    case RunOptions::UpscaledMixingModel::Matern:
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::matern;
        break;
    case RunOptions::UpscaledMixingModel::OneOverSum:
    default:
        g_u.VelocityCorrelationModel = Grid2D::velocity_correlation_model::oneoversum;
        break;
    }

    g_u.setDiffusionFactor(P.diffusion_factor);

    int qx_size = nu * (nx + 1);
    int qy_size = (nu + 1) * nx;

    auto& qx_u = g_u.flux("qx");
    auto& qy_u = g_u.flux("qy");

    qx_u.resize(qx_size, 0.0);
    qy_u.resize(qy_size, 0.0);

    for (int j = 0; j < nu; ++j) {
        double u = (nu <= 1) ? 0.0 : static_cast<double>(j) / (nu - 1);
        double v_at_u = invcdf_mean.interpol(u);

        for (int i = 0; i < nx + 1; ++i) {
            int id = j * (nx + 1) + i;
            if (id >= qx_size) return false;
            qx_u[id] = v_at_u;
        }
    }

    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));

    g_u.setMixingParams(lc_mean, lx_mean, ly_mean, nu_x_mean, nu_y_mean);

    g_u.SetVal("diffusion", P.Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));

    // FIX: SolveTransport below uses "Cu" -> make sure it exists/initialized
    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0.0);

    g_u.setBTCLocations(P.xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    if (opts.solve_upscale_transport) {
        g_u.SolveTransport(P.t_end_pdf, dt_mean, "transport_", 500, up_dir, "C", &BTCs_Upscaled);

        TimeSeriesSet<double> BTCs_Upscaled_deriv = BTCs_Upscaled.derivative();

        // Append RAW upscaled transport BTC to explicit CDF stack
        for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_Upscaled.size(); ++i) {
            BTCs_Upscaled.setSeriesName(i, fmt_x(P.xLocations[i]));
            Fine_Scale_BTCs_transport_cdf[i].append(BTCs_Upscaled[i], "Upscaled" + fmt_x(P.xLocations[i]));
        }

        // Append derivative upscaled transport BTC to explicit PDF stack
        for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_Upscaled_deriv.size(); ++i) {
            BTCs_Upscaled_deriv.setSeriesName(i, fmt_x(P.xLocations[i]));
            Fine_Scale_BTCs_transport_pdf[i].append(BTCs_Upscaled_deriv[i], "Upscaled" + fmt_x(P.xLocations[i]));
        }

        // Per-upscaled-run transport files:
        //   BTC_Upscaled.csv            -> raw/CDF-style transport BTC
        //   BTC_Upscaled_derivative.csv -> PDF-style transport BTC
        BTCs_Upscaled.write(up_btc_path);
        BTCs_Upscaled.derivative().write(up_btc_deriv_path);
    }

    // PT (SEPARATE)
    if (opts.perform_upscaled_PT) {
        const int n_particles = 100000;
        const double dx_step = 0.01;

        PathwaySet pathways;
        pathways.Initialize(n_particles, PathwaySet::Weighting::FluxWeighted, &g_u);
        pathways.trackAllPathwaysWithDiffusion(&g_u, dx_step, g_u.getDiffusion());

        const std::vector<double>& btc_locations = g_u.getBTCLocations();
        if (!btc_locations.empty()) {
            TimeSeriesSet<double> btc_pdf = pathways.getBreakthroughCurve(
                btc_locations, true, PathwaySet::BTCType::PDF, 200);
            btc_pdf.write(joinPath(up_dir, up_pfx + "PT_btc_pdf.csv"));

            TimeSeriesSet<double> btc_cdf = pathways.getBreakthroughCurve(
                btc_locations, true, PathwaySet::BTCType::CDF);
            btc_cdf.write(joinPath(up_dir, up_pfx + "PT_btc_cdf.csv"));

            const int nmin = std::min((int)P.xLocations.size(), (int)btc_pdf.size());
            for (int i = 0; i < nmin; ++i)
                Fine_Scale_PT_pdf[i].append(btc_pdf[i], "Upscaled_PT" + fmt_x(P.xLocations[i]));

            const int nmin2 = std::min((int)P.xLocations.size(), (int)btc_cdf.size());
            for (int i = 0; i < nmin2; ++i)
                Fine_Scale_PT_cdf[i].append(btc_cdf[i], "Upscaled_PT" + fmt_x(P.xLocations[i]));
        }
    }

    up_dir_out = up_dir;
    up_btc_path_out = up_btc_path;
    up_btc_deriv_path_out = up_btc_deriv_path;
    return true;
}

// ============================================================================
// Rebuild from existing BTC/PT files (no re-solve)
// ============================================================================

static bool rebuild_from_existing_btc_files(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    RunOutputs& out,
    int rank)
{
    if (rank == 0) {
        std::cout << "Rebuild mode: reading BTC/PT from existing files in:\n"
                  << "  " << run_dir << "\n";
    }

    // -----------------------------
    // Fine-scale folders
    // -----------------------------
    const std::vector<std::string> fine_dirs = list_fine_dirs_sorted(run_dir);

    for (const auto& fine_dir_raw : fine_dirs) {
        std::string rlab;
        if (!fine_dir_to_rlab(fine_dir_raw, rlab)) continue;

        std::string fine_dir = fine_dir_raw;
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        const std::string pfx = rlab + "_";

        TimeSeriesSet<double> btc_raw;
        TimeSeriesSet<double> btc_deriv;

        const std::string btc_raw_path   = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
        const std::string btc_deriv_path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

        const bool has_btc_raw   = read_tsset_file_if_exists(btc_raw_path, btc_raw);
        const bool has_btc_deriv = read_tsset_file_if_exists(btc_deriv_path, btc_deriv);

        if (has_btc_raw || has_btc_deriv) {
            if (!has_btc_raw && has_btc_deriv) {
                if (rank == 0) {
                    std::cerr << "WARNING: raw fine BTC missing; skipping CDF append for " << fine_dir << "\n";
                }
            }
            if (has_btc_raw && !has_btc_deriv) {
                btc_deriv = btc_raw.derivative();
            }

            if (has_btc_raw) {
                for (int i = 0; i < (int)P.xLocations.size() && i < (int)btc_raw.size(); ++i) {
                    btc_raw.setSeriesName(i, fmt_x(P.xLocations[i]));
                    out.Fine_Scale_BTCs_cdf[i].append(btc_raw[i], pfx + fmt_x(P.xLocations[i]));
                }
            }

            for (int i = 0; i < (int)P.xLocations.size() && i < (int)btc_deriv.size(); ++i) {
                btc_deriv.setSeriesName(i, fmt_x(P.xLocations[i]));
                out.Fine_Scale_BTCs_pdf[i].append(btc_deriv[i], pfx + fmt_x(P.xLocations[i]));
            }
        }

        // optional PT
        TimeSeriesSet<double> pt_pdf, pt_cdf;
        const std::string pt_pdf_path = joinPath(fine_dir, pfx + "btc_pdf.csv");
        const std::string pt_cdf_path = joinPath(fine_dir, pfx + "btc_cdf.csv");

        if (read_tsset_file_if_exists(pt_pdf_path, pt_pdf)) {
            for (int i = 0; i < (int)P.xLocations.size() && i < (int)pt_pdf.size(); ++i) {
                pt_pdf.setSeriesName(i, fmt_x(P.xLocations[i]));
                out.Fine_Scale_PT_pdf[i].append(pt_pdf[i], pfx + fmt_x(P.xLocations[i]));
            }
        }

        if (read_tsset_file_if_exists(pt_cdf_path, pt_cdf)) {
            for (int i = 0; i < (int)P.xLocations.size() && i < (int)pt_cdf.size(); ++i) {
                pt_cdf.setSeriesName(i, fmt_x(P.xLocations[i]));
                out.Fine_Scale_PT_cdf[i].append(pt_cdf[i], pfx + fmt_x(P.xLocations[i]));
            }
        }
    }

    // -----------------------------
    // Upscaled files
    // -----------------------------
    const std::string up_dir = joinPath(run_dir, "upscaled_mean");
    const std::string up_btc_path       = joinPath(up_dir, "upscaled_BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, "upscaled_BTC_Upscaled_derivative.csv");

    TimeSeriesSet<double> up_btc_raw;
    TimeSeriesSet<double> up_btc_deriv;

    const bool has_up_btc_raw   = read_tsset_file_if_exists(up_btc_path, up_btc_raw);
    const bool has_up_btc_deriv = read_tsset_file_if_exists(up_btc_deriv_path, up_btc_deriv);

    if (has_up_btc_raw || has_up_btc_deriv) {
        if (!has_up_btc_raw && has_up_btc_deriv) {
            if (rank == 0) {
                std::cerr << "WARNING: raw upscaled BTC missing; skipping upscaled CDF append.\n";
            }
        }
        if (has_up_btc_raw && !has_up_btc_deriv) {
            up_btc_deriv = up_btc_raw.derivative();
        }

        if (has_up_btc_raw) {
            for (int i = 0; i < (int)P.xLocations.size() && i < (int)up_btc_raw.size(); ++i) {
                up_btc_raw.setSeriesName(i, fmt_x(P.xLocations[i]));
                out.Fine_Scale_BTCs_cdf[i].append(up_btc_raw[i], "Upscaled" + fmt_x(P.xLocations[i]));
            }
        }

        for (int i = 0; i < (int)P.xLocations.size() && i < (int)up_btc_deriv.size(); ++i) {
            up_btc_deriv.setSeriesName(i, fmt_x(P.xLocations[i]));
            out.Fine_Scale_BTCs_pdf[i].append(up_btc_deriv[i], "Upscaled" + fmt_x(P.xLocations[i]));
        }
    }

    // optional upscaled PT
    {
        TimeSeriesSet<double> up_pt_pdf, up_pt_cdf;
        const std::string up_pt_pdf_path = joinPath(up_dir, "upscaled_PT_btc_pdf.csv");
        const std::string up_pt_cdf_path = joinPath(up_dir, "upscaled_PT_btc_cdf.csv");

        if (read_tsset_file_if_exists(up_pt_pdf_path, up_pt_pdf)) {
            for (int i = 0; i < (int)P.xLocations.size() && i < (int)up_pt_pdf.size(); ++i) {
                up_pt_pdf.setSeriesName(i, fmt_x(P.xLocations[i]));
                out.Fine_Scale_PT_pdf[i].append(up_pt_pdf[i], "Upscaled_PT" + fmt_x(P.xLocations[i]));
            }
        }

        if (read_tsset_file_if_exists(up_pt_cdf_path, up_pt_cdf)) {
            for (int i = 0; i < (int)P.xLocations.size() && i < (int)up_pt_cdf.size(); ++i) {
                up_pt_cdf.setSeriesName(i, fmt_x(P.xLocations[i]));
                out.Fine_Scale_PT_cdf[i].append(up_pt_cdf[i], "Upscaled_PT" + fmt_x(P.xLocations[i]));
            }
        }
    }

    out.Fine_Scale_BTCs = out.Fine_Scale_BTCs_pdf;

    computeMeanBTCs(P.xLocations, out.Fine_Scale_BTCs_pdf, out.mean_transport_pdf, opts);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_BTCs_cdf, out.mean_transport_cdf, opts);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_PT_pdf,   out.mean_pt_pdf,        opts);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_PT_cdf,   out.mean_pt_cdf,        opts);

    out.mean_transport_full = out.mean_transport_pdf;
    out.mean_BTCs = out.mean_transport_pdf;

    out.up_dir = up_dir;
    out.up_btc_path = up_btc_path;
    out.up_btc_deriv_path = up_btc_deriv_path;

    return true;
}

// ============================================================================
// Main orchestrator
//
// Transport convention propagated to RunOutputs:
//
//   Explicit stacks:
//     Fine_Scale_BTCs_pdf  -> PDF transport compare stack
//     Fine_Scale_BTCs_cdf  -> CDF/raw transport compare stack
//
//   Legacy stack:
//     Fine_Scale_BTCs      -> alias to PDF transport stack
//
//   Explicit means:
//     mean_transport_pdf
//     mean_transport_cdf
//
//   Legacy means:
//     mean_BTCs            -> alias to PDF transport mean
//     mean_transport_full  -> currently also alias to PDF transport mean
//
// This preserves old filenames while changing their content to PDF.
// ============================================================================

bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& resume_run_dir,
    RunOutputs& out,
    int rank)
{
    // Transport compare (legacy names should map to PDF later in main)
    out.Fine_Scale_BTCs.clear();
    out.Fine_Scale_BTCs.resize(P.xLocations.size());

    // Explicit transport split
    out.Fine_Scale_BTCs_pdf.clear();
    out.Fine_Scale_BTCs_cdf.clear();
    out.Fine_Scale_BTCs_pdf.resize(P.xLocations.size());
    out.Fine_Scale_BTCs_cdf.resize(P.xLocations.size());

    // PT compare (separate)
    out.Fine_Scale_PT_pdf.clear();
    out.Fine_Scale_PT_cdf.clear();
    out.Fine_Scale_PT_pdf.resize(P.xLocations.size());
    out.Fine_Scale_PT_cdf.resize(P.xLocations.size());

    // NEW: rebuild/recovery mode from already-written files
    if (opts.read_btc_from_files) {
        const bool ok_rebuild = rebuild_from_existing_btc_files(P, opts, resume_run_dir, out, rank);
        if (!ok_rebuild) return false;

        if (rank == 0) {
            for (int i = 0; i < (int)P.xLocations.size(); ++i) {
                out.Fine_Scale_PT_pdf[i].write(joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "PT_Compare_pdf.csv"));
                out.Fine_Scale_PT_cdf[i].write(joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "PT_Compare_cdf.csv"));
            }
        }
        return true;
    }

    FineScaleOutputs fine_outputs;
    fine_outputs.BTCs_transport_pdf.resize(P.xLocations.size());
    fine_outputs.BTCs_transport_cdf.resize(P.xLocations.size());
    fine_outputs.PT_pdfs.resize(P.xLocations.size());
    fine_outputs.PT_cdfs.resize(P.xLocations.size());

    if (!opts.upscale_only) {
        const bool ok = run_fine_loop_collect(P, opts, out.run_dir, fine_outputs, rank);
        if (!ok) return false;
    }

    out.Fine_Scale_BTCs_pdf = fine_outputs.BTCs_transport_pdf;
    out.Fine_Scale_BTCs_cdf = fine_outputs.BTCs_transport_cdf;

    // Legacy transport compare now points to PDF stack
    out.Fine_Scale_BTCs = out.Fine_Scale_BTCs_pdf;

    out.Fine_Scale_PT_pdf = fine_outputs.PT_pdfs;
    out.Fine_Scale_PT_cdf = fine_outputs.PT_cdfs;

    TimeSeries<double> invcdf_mean;
    const bool ok_mean = build_mean_for_upscaled(P, opts, H, out.run_dir, resume_run_dir, fine_outputs, invcdf_mean, rank);
    if (!ok_mean) return false;

    broadcastMeanParameters(fine_outputs, invcdf_mean, rank);

    std::string up_dir, up_btc_path, up_btc_deriv_path;
    const bool ok_up = run_upscaled(
        P, opts, out.run_dir,
        fine_outputs.lc_mean,
        fine_outputs.velocity_lx_mean,
        fine_outputs.velocity_ly_mean,
        fine_outputs.nu_x_mean,
        fine_outputs.nu_y_mean,
        fine_outputs.dt_mean,
        invcdf_mean,
        up_dir, up_btc_path, up_btc_deriv_path,
        out.Fine_Scale_BTCs_pdf,
        out.Fine_Scale_BTCs_cdf,
        out.Fine_Scale_PT_pdf,
        out.Fine_Scale_PT_cdf,
        rank
    );
    if (!ok_up) return false;

    // Legacy transport compare now points to PDF stack after upscaled append too
    out.Fine_Scale_BTCs = out.Fine_Scale_BTCs_pdf;

    out.up_dir = up_dir;
    out.up_btc_path = up_btc_path;
    out.up_btc_deriv_path = up_btc_deriv_path;

    // Means (SEPARATED) — switchable via opts.mean_ts_mode
    computeMeanBTCs(P.xLocations, out.Fine_Scale_BTCs_pdf, out.mean_transport_pdf, opts);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_BTCs_cdf, out.mean_transport_cdf, opts);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_PT_pdf,   out.mean_pt_pdf,        opts);
    computeMeanBTCs(P.xLocations, out.Fine_Scale_PT_cdf,   out.mean_pt_cdf,        opts);

    // Legacy transport mean now points to PDF mean
    out.mean_transport_full = out.mean_transport_pdf;
    out.mean_BTCs = out.mean_transport_pdf;

    // Write per-x PT compare files
    //
    // Transport compare files are intentionally NOT written here.
    // They are written in main.cpp so that:
    //   - legacy transport names can point to PDF
    //   - explicit transport _pdf/_cdf files can both be written
    if (rank == 0) {
        for (int i = 0; i < (int)P.xLocations.size(); ++i) {
            out.Fine_Scale_PT_pdf[i].write(joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "PT_Compare_pdf.csv"));
            out.Fine_Scale_PT_cdf[i].write(joinPath(out.run_dir, fmt_x(P.xLocations[i]) + "PT_Compare_cdf.csv"));
        }
    }

    return true;
}
