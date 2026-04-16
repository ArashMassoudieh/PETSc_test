#include "copula_analysis.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <utility>

#include <gsl/gsl_cdf.h>

// ============================================================================
// copula_analysis.cpp
//
// Workflow-facing copula analysis layer used by sim_runner.
// This file provides:
//   - empirical copula diagnostics and bootstrap tests
//   - Gaussian copula diagnostics and GOF testing
//   - rank/pseudo-observation transforms
//   - helpers to build CopulaPairData from PathwaySet samples
//   - CSV / VTK export for copula point clouds and binned matrices
// ============================================================================

// ============================================================================
// File-local empirical copula helpers
// These stay internal to this translation unit.
// ============================================================================
namespace
{
    // Empirical copula CDF evaluated by direct counting:
    //   C(u,v) = P(U <= u, V <= v)
    // Complexity: O(n)
    double empirical_copula_cdf_naive(
        const std::vector<double>& u1,
        const std::vector<double>& u2,
        double u,
        double v)
    {
        const int n = (int)std::min(u1.size(), u2.size());
        if (n <= 0) return std::numeric_limits<double>::quiet_NaN();

        int cnt = 0;
        for (int i = 0; i < n; ++i) {
            if (u1[i] <= u && u2[i] <= v) ++cnt;
        }
        return double(cnt) / double(n);
    }

    // Cramér–von Mises-type statistic versus the independence copula C0(u,v)=u*v.
    // Larger values indicate stronger departure from independence.
    double empirical_copula_cvm_vs_independence(
        const std::vector<double>& u1,
        const std::vector<double>& u2)
    {
        const int n = (int)std::min(u1.size(), u2.size());
        if (n <= 1) return std::numeric_limits<double>::quiet_NaN();

        double s = 0.0;
        for (int i = 0; i < n; ++i) {
            const double cn = empirical_copula_cdf_naive(u1, u2, u1[i], u2[i]);
            const double c0 = u1[i] * u2[i];
            if (std::isfinite(cn) && std::isfinite(c0)) {
                const double d = cn - c0;
                s += d * d;
            }
        }
        return s / double(n);
    }

    // Fraction of samples simultaneously in the upper tail:
    //   U > q and V > q
    double empirical_copula_upper_tail_fraction(
        const std::vector<double>& u1,
        const std::vector<double>& u2,
        double q)
    {
        const int n = (int)std::min(u1.size(), u2.size());
        if (n <= 0) return std::numeric_limits<double>::quiet_NaN();

        int cnt = 0;
        for (int i = 0; i < n; ++i) {
            if (u1[i] > q && u2[i] > q) ++cnt;
        }
        return double(cnt) / double(n);
    }

    // L1 deviation of the empirical copula along the diagonal u=v.
    // Useful as a simple summary of departure from perfect dependence.
    double empirical_copula_diagonal_l1(
        const std::vector<double>& u1,
        const std::vector<double>& u2)
    {
        const int n = (int)std::min(u1.size(), u2.size());
        if (n <= 1) return std::numeric_limits<double>::quiet_NaN();

        double s = 0.0;
        for (int k = 1; k <= n; ++k) {
            const double t = (double(k) - 0.5) / double(n);
            const double cn = empirical_copula_cdf_naive(u1, u2, t, t);
            s += std::abs(cn - t);
        }
        return s / double(n);
    }

    // Bootstrap p-value for the null hypothesis of independence.
    // Simulates independent uniforms and compares the empirical CvM statistic.
    double estimate_empirical_copula_independence_pvalue(
        const std::vector<double>& u1,
        const std::vector<double>& u2,
        int B,
        unsigned long seed,
        double& stat_obs_out)
    {
        stat_obs_out = empirical_copula_cvm_vs_independence(u1, u2);
        if (!std::isfinite(stat_obs_out) || B <= 0) return std::numeric_limits<double>::quiet_NaN();

        const int n = (int)std::min(u1.size(), u2.size());
        if (n <= 1) return std::numeric_limits<double>::quiet_NaN();

        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> U01(1e-12, 1.0 - 1e-12);

        int ge = 0;
        std::vector<double> a(n), b(n);
        for (int bb = 0; bb < B; ++bb) {
            for (int i = 0; i < n; ++i) {
                a[i] = U01(rng);
                b[i] = U01(rng);
            }
            const double stat_b = empirical_copula_cvm_vs_independence(a, b);
            if (std::isfinite(stat_b) && stat_b >= stat_obs_out) ++ge;
        }
        return double(ge + 1) / double(B + 1);
    }
}

// ============================================================================
// Basic statistical helpers
// ============================================================================

// Pearson correlation on two equal-length vectors.
double pearson_corr(const std::vector<double>& a, const std::vector<double>& b)
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

// Convert raw values to fractional ranks in (0,1).
// Ties receive the average rank.
std::vector<double> fractional_ranks_01(const std::vector<double>& x)
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
        const double rank_avg = 0.5 * (double(i + 1) + double(j));
        for (int k = i; k < j; ++k) {
            r[xv[k].second] = rank_avg / (double(n) + 1.0);
        }
        i = j;
    }
    return r;
}

// Approximate inverse standard-normal CDF.
// Used to map pseudo-observations to Gaussian normal scores.
double norminv_approx(double p)
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

// Approximate Gaussian copula CDF C(u,v; rho).
// Uses numerical integration over the bivariate normal representation.
double gaussian_copula_cdf_approx(double u, double v, double rho)
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
        sum += w * phi(x) * gsl_cdf_ugaussian_P(arg);
    }
    return std::max(0.0, std::min(1.0, sum * h / 3.0));
}

// Naive O(n^2) Kendall tau implementation.
double kendall_tau_naive(const std::vector<double>& x, const std::vector<double>& y)
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

// Mardia skewness and kurtosis for a bivariate sample.
// These are diagnostic summaries for departure from multivariate normality.
void mardia_bivariate_moments(const std::vector<double>& x,
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

// CvM-type goodness-of-fit statistic for a Gaussian copula with parameter rho.
double gaussian_copula_cvm_statistic(const std::vector<double>& u,
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

// Bootstrap p-value for Gaussian copula goodness of fit.
double estimate_gaussian_copula_gof_pvalue(const std::vector<double>& u,
                                           const std::vector<double>& v,
                                           double rho,
                                           int B,
                                           unsigned long seed,
                                           double& stat_obs_out)
{
    constexpr double kPi = 3.14159265358979323846;
    stat_obs_out = gaussian_copula_cvm_statistic(u, v, rho);
    if (!std::isfinite(stat_obs_out) || B <= 0) return std::numeric_limits<double>::quiet_NaN();

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> N01(0.0, 1.0);
    int ge = 0;
    const int n = (int)u.size();
    const double rho_clamped = std::max(-0.999, std::min(0.999, rho));
    const double s = std::sqrt(1.0 - rho_clamped * rho_clamped);
    std::vector<double> ub(n), vb(n);
    for (int b = 0; b < B; ++b) {
        for (int i = 0; i < n; ++i) {
            const double z1 = N01(rng);
            const double z2 = rho_clamped * z1 + s * N01(rng);
            ub[i] = std::min(1.0 - 1e-12, std::max(1e-12, gsl_cdf_ugaussian_P(z1)));
            vb[i] = std::min(1.0 - 1e-12, std::max(1e-12, gsl_cdf_ugaussian_P(z2)));
        }
        const double tau_b = kendall_tau_naive(ub, vb);
        const double rho_b = std::sin(0.5 * kPi * tau_b);
        const double s_b = gaussian_copula_cvm_statistic(ub, vb, rho_b);
        if (std::isfinite(s_b) && s_b >= stat_obs_out) ++ge;
    }
    return double(ge + 1) / double(B + 1);
}

// Evenly downsample a vector to at most max_points entries.
// Used to keep diagnostics/bootstraps from getting too expensive.
std::vector<double> downsample_evenly(const std::vector<double>& a, int max_points)
{
    if (max_points <= 0 || (int)a.size() <= max_points) return a;
    std::vector<double> out;
    out.reserve(max_points);
    const int n = (int)a.size();
    for (int k = 0; k < max_points; ++k) {
        int idx = (int)std::llround((double)k * (n - 1) / (max_points - 1));
        idx = std::max(0, std::min(n - 1, idx));
        out.push_back(a[idx]);
    }
    return out;
}

// Build raw values, pseudo-observations, and Gaussian normal scores
// from a sampled pair of pathways.
CopulaPairData build_copula_pair_data(const PathwaySet& pair_set)
{
    CopulaPairData data;
    if (pair_set.size() < 2) return data;

    const auto& p1 = pair_set[0];
    const auto& p2 = pair_set[1];
    const int n = (int)std::min(p1.size(), p2.size());
    if (n < 2) return data;

    data.q1.resize(n);
    data.q2.resize(n);
    for (int i = 0; i < n; ++i) {
        data.q1[i] = p1[i].qx();
        data.q2[i] = p2[i].qx();
    }

    data.u1 = fractional_ranks_01(data.q1);
    data.u2 = fractional_ranks_01(data.q2);

    data.z1.resize(n, std::numeric_limits<double>::quiet_NaN());
    data.z2.resize(n, std::numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < n; ++i) {
        data.z1[i] = norminv_approx(data.u1[i]);
        data.z2[i] = norminv_approx(data.u2[i]);
    }

    data.n_pairs = n;
    return data;
}

// Build a normalized binned empirical copula matrix on the unit square.
CopulaBinnedMatrix build_empirical_binned_copula(
    const std::vector<double>& u1,
    const std::vector<double>& u2,
    int nBins)
{
    CopulaBinnedMatrix M(nBins);
    if (nBins <= 0) return M;
    M.accumulateUnitPairs(u1, u2);
    M.normalizeToUnitMass();
    return M;
}

// Run the selected dependence diagnostics on prepared copula pair data.
CopulaDiagnostics analyze_copula_pair_data(
    const CopulaPairData& data,
    double delta_x,
    const CopulaAnalysisOptions& opts)
{
    CopulaDiagnostics out;
    out.delta_x = delta_x;
    const int n = (int)std::min(data.u1.size(), data.u2.size());
    if (n < 2) return out;

    std::vector<double> z1_finite;
    std::vector<double> z2_finite;
    z1_finite.reserve(n);
    z2_finite.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (std::isfinite(data.z1[i]) && std::isfinite(data.z2[i])) {
            z1_finite.push_back(data.z1[i]);
            z2_finite.push_back(data.z2[i]);
        }
    }

    out.n_pairs = n;
    out.corr_qx = pearson_corr(data.q1, data.q2);
    out.corr_rank = pearson_corr(data.u1, data.u2);
    out.gaussian_copula_rho = pearson_corr(z1_finite, z2_finite);
    out.selected_rank_dependence =
        (opts.dependence_model == CopulaDependenceModel::Empirical)
            ? out.corr_rank
            : out.gaussian_copula_rho;

    if (opts.compute_diagnostics) {
        const std::vector<double> q1s = downsample_evenly(data.q1, opts.max_points);
        const std::vector<double> q2s = downsample_evenly(data.q2, opts.max_points);
        const std::vector<double> u1s = downsample_evenly(data.u1, opts.max_points);
        const std::vector<double> u2s = downsample_evenly(data.u2, opts.max_points);

        out.kendall_tau = kendall_tau_naive(q1s, q2s);
        out.rho_from_tau = std::sin(0.5 * 3.14159265358979323846 * out.kendall_tau);
        mardia_bivariate_moments(q1s, q2s, out.mardia_skewness, out.mardia_kurtosis);

        out.gaussian_copula_gof_pvalue =
            estimate_gaussian_copula_gof_pvalue(
                u1s, u2s, out.rho_from_tau,
                std::max(10, opts.bootstrap_B),
                777UL + (unsigned long)(1000.0 * delta_x),
                out.gaussian_copula_gof_stat);

        out.empirical_copula_pvalue =
            estimate_empirical_copula_independence_pvalue(
                u1s, u2s,
                std::max(10, opts.bootstrap_B),
                1777UL + (unsigned long)(1000.0 * delta_x),
                out.empirical_copula_stat);

        out.empirical_upper_tail_frac_90 = empirical_copula_upper_tail_fraction(u1s, u2s, 0.9);
        out.empirical_diagonal_l1 = empirical_copula_diagonal_l1(u1s, u2s);
        out.selected_rank_dependence =
            (opts.dependence_model == CopulaDependenceModel::Empirical)
                ? out.corr_rank
                : out.gaussian_copula_rho;
    }

    return out;
}

// Write raw/rank/Gaussian-transformed pair data to CSV.
bool write_rank_pairs_csv(const CopulaPairData& data, const std::string& filename)
{
    const int n = (int)std::min({data.q1.size(), data.q2.size(), data.u1.size(), data.u2.size(), data.z1.size(), data.z2.size()});
    if (n <= 0) return false;

    std::ofstream f(filename.c_str());
    if (!f.good()) return false;

    f << "u1,u2,qx1,qx2,qx1_normal_score,qx2_normal_score\n";
    f << std::setprecision(15);
    for (int i = 0; i < n; ++i) {
        f << data.u1[i] << "," << data.u2[i] << ","
          << data.q1[i] << "," << data.q2[i] << ","
          << data.z1[i] << "," << data.z2[i] << "\n";
    }
    return true;
}

// Write a 2D matrix to VTK ImageData (.vti) on the unit square.
bool write_matrix_as_vti_2d(const CMatrix& M,
                            const std::string& filename,
                            const std::string& array_name,
                            bool point_data)
{
    const int nx = M.getnumrows();
    const int ny = M.getnumcols();
    if (nx <= 0 || ny <= 0) return false;

    std::ofstream f(filename.c_str());
    if (!f.good()) return false;

    const int ex_i1 = point_data ? (nx - 1) : nx;
    const int ex_j1 = point_data ? (ny - 1) : ny;
    const double dx = 1.0 / static_cast<double>(nx);
    const double dy = 1.0 / static_cast<double>(ny);

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <ImageData Origin=\"0 0 0\" "
      << "Spacing=\"" << dx << " " << dy << " 1\" "
      << "WholeExtent=\"0 " << ex_i1 << " 0 " << ex_j1 << " 0 0\">\n";
    f << "    <Piece Extent=\"0 " << ex_i1 << " 0 " << ex_j1 << " 0 0\">\n";

    if (point_data) f << "      <PointData Scalars=\"" << array_name << "\">\n";
    else            f << "      <CellData Scalars=\"" << array_name << "\">\n";

    f << "        <DataArray type=\"Float64\" Name=\"" << array_name << "\" format=\"ascii\">\n";
    f << std::setprecision(17);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double v = M[i][j];
            if (!std::isfinite(v)) v = 0.0;
            f << v << "\n";
        }
    }
    f << "        </DataArray>\n";

    if (point_data) f << "      </PointData>\n";
    else            f << "      </CellData>\n";

    f << "      <FieldData>\n";
    f << "        <DataArray type=\"Float64\" Name=\"xmin\" NumberOfTuples=\"1\" format=\"ascii\">0</DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"xmax\" NumberOfTuples=\"1\" format=\"ascii\">1</DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"ymin\" NumberOfTuples=\"1\" format=\"ascii\">0</DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"ymax\" NumberOfTuples=\"1\" format=\"ascii\">1</DataArray>\n";
    f << "      </FieldData>\n";

    f << "    </Piece>\n";
    f << "  </ImageData>\n";
    f << "</VTKFile>\n";
    return true;
}

// Write copula pseudo-observation points to VTK PolyData (.vtp).
// Optional scalar arrays can be attached for raw qx values and normal scores.
bool write_rank_points_as_vtp(const std::vector<double>& u1,
                              const std::vector<double>& u2,
                              const std::vector<double>* qx1,
                              const std::vector<double>* qx2,
                              const std::vector<double>* z1,
                              const std::vector<double>* z2,
                              const std::string& filename)
{
    const int n = (int)std::min(u1.size(), u2.size());
    if (n <= 0) return false;

    auto same_or_null = [n](const std::vector<double>* v) {
        return (!v) || ((int)v->size() >= n);
    };
    if (!same_or_null(qx1) || !same_or_null(qx2) || !same_or_null(z1) || !same_or_null(z2))
        return false;

    std::ofstream f(filename.c_str());
    if (!f.good()) return false;

    auto write_scalar = [&](const std::string& name, const std::vector<double>& a)
    {
        f << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        f << std::setprecision(17);
        for (int i = 0; i < n; ++i) {
            double v = a[i];
            if (!std::isfinite(v)) v = 0.0;
            f << v << "\n";
        }
        f << "        </DataArray>\n";
    };

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <PolyData>\n";
    f << "    <Piece NumberOfPoints=\"" << n
      << "\" NumberOfVerts=\"" << n
      << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    f << "      <PointData Scalars=\"u1\">\n";
    write_scalar("u1", u1);
    write_scalar("u2", u2);
    if (qx1) write_scalar("qx1", *qx1);
    if (qx2) write_scalar("qx2", *qx2);
    if (z1)  write_scalar("qx1_normal_score", *z1);
    if (z2)  write_scalar("qx2_normal_score", *z2);

    std::vector<double> diagdist(n, 0.0), avg_u(n, 0.0);
    for (int i = 0; i < n; ++i) {
        diagdist[i] = std::abs(u1[i] - u2[i]);
        avg_u[i] = 0.5 * (u1[i] + u2[i]);
    }
    write_scalar("diag_distance", diagdist);
    write_scalar("avg_u", avg_u);
    f << "      </PointData>\n";
    f << "      <CellData/>\n";

    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    f << std::setprecision(17);
    for (int i = 0; i < n; ++i) {
        const double x = std::isfinite(u1[i]) ? u1[i] : 0.0;
        const double y = std::isfinite(u2[i]) ? u2[i] : 0.0;
        f << x << " " << y << " 0\n";
    }
    f << "        </DataArray>\n";
    f << "      </Points>\n";

    f << "      <Verts>\n";
    f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < n; ++i) f << i << "\n";
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < n; ++i) f << (i + 1) << "\n";
    f << "        </DataArray>\n";
    f << "      </Verts>\n";

    f << "    </Piece>\n";
    f << "  </PolyData>\n";
    f << "</VTKFile>\n";
    return true;
}

// Convenience overload using a CopulaPairData bundle.
bool write_rank_points_as_vtp(const CopulaPairData& data,
                              const std::string& filename)
{
    return write_rank_points_as_vtp(data.u1, data.u2, &data.q1, &data.q2, &data.z1, &data.z2, filename);
}

// High-level helper used by sim_runner:
//   1) build copula data from a sampled pathway pair set
//   2) optionally build a binned empirical copula matrix
//   3) write CSV + VTP outputs
//   4) return diagnostics
CopulaDiagnostics analyze_and_write_rank_pairs(
    const PathwaySet& pair_set,
    double delta_x,
    const std::string& scatter_csv_path,
    const CopulaAnalysisOptions& opts,
    CopulaBinnedMatrix* empirical_binned_out)
{
    CopulaDiagnostics out;
    CopulaPairData data = build_copula_pair_data(pair_set);
    out.delta_x = delta_x;
    if (data.n_pairs < 2) return out;

    if (empirical_binned_out && opts.empirical_copula_bins > 0) {
        *empirical_binned_out = build_empirical_binned_copula(data.u1, data.u2, opts.empirical_copula_bins);
    }

    write_rank_pairs_csv(data, scatter_csv_path);

    std::string vtp_path = scatter_csv_path;
    const std::string ext = ".csv";
    if (vtp_path.size() >= ext.size() &&
        vtp_path.substr(vtp_path.size() - ext.size()) == ext)
        vtp_path.replace(vtp_path.size() - ext.size(), ext.size(), ".vtp");
    else
        vtp_path += ".vtp";

    write_rank_points_as_vtp(data, vtp_path);
    return analyze_copula_pair_data(data, delta_x, opts);
}
