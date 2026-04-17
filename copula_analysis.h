#pragma once

#include <string>
#include <vector>
#include <limits>

#include "Matrix.h"
#include "PathwaySet.h"
#include "TimeSeries.h"

// ============================================================================
// Holds paired data for copula analysis between two variables (e.g., qx1, qx2)
// Includes:
//   - raw values (q1, q2)
//   - rank-transformed values in [0,1] (u1, u2)
//   - Gaussian-transformed values (z1, z2) via inverse normal CDF
// ============================================================================
struct CopulaPairData
{
    std::vector<double> q1;   // Raw variable 1 values (e.g., qx at location A)
    std::vector<double> q2;   // Raw variable 2 values (e.g., qx at location B)

    std::vector<double> u1;   // Rank-transformed values of q1 in (0,1)
    std::vector<double> u2;   // Rank-transformed values of q2 in (0,1)

    std::vector<double> z1;   // Gaussianized values: z1 = Phi^{-1}(u1)
    std::vector<double> z2;   // Gaussianized values: z2 = Phi^{-1}(u2)

    int n_pairs = 0;          // Number of valid paired samples
};

// ============================================================================
// Stores diagnostics and goodness-of-fit metrics for a copula analysis
// ============================================================================
struct CopulaDiagnostics
{
    double delta_x = 0.0;     // Spatial separation (or lag distance)

    int n_pairs = 0;          // Number of data pairs used

    // --- Basic correlations ---
    double corr_qx = std::numeric_limits<double>::quiet_NaN();     // Pearson correlation on raw values
    double corr_rank = std::numeric_limits<double>::quiet_NaN();   // Pearson correlation on rank-transformed values

    // --- Gaussian copula fit ---
    double gaussian_copula_rho = std::numeric_limits<double>::quiet_NaN();         // Estimated correlation parameter
    double gaussian_copula_gof_stat = std::numeric_limits<double>::quiet_NaN();    // CvM or similar GOF statistic
    double gaussian_copula_gof_pvalue = std::numeric_limits<double>::quiet_NaN();  // Bootstrap p-value

    // --- Rank-based dependence ---
    double kendall_tau = std::numeric_limits<double>::quiet_NaN(); // Kendall’s tau
    double rho_from_tau = std::numeric_limits<double>::quiet_NaN();// Implied Gaussian rho from tau

    // --- Multivariate normality diagnostics (on z1, z2) ---
    double mardia_skewness = std::numeric_limits<double>::quiet_NaN();  // Mardia skewness
    double mardia_kurtosis = std::numeric_limits<double>::quiet_NaN();  // Mardia kurtosis

    // --- Empirical copula diagnostics ---
    double empirical_copula_stat = std::numeric_limits<double>::quiet_NaN();       // GOF stat vs empirical copula
    double empirical_copula_pvalue = std::numeric_limits<double>::quiet_NaN();     // Bootstrap p-value

    double empirical_upper_tail_frac_90 = std::numeric_limits<double>::quiet_NaN(); // Fraction of points in upper tail (u>0.9, v>0.9)
    double empirical_diagonal_l1 = std::numeric_limits<double>::quiet_NaN();        // L1 deviation from diagonal (dependence strength)

    // --- Selected dependence measure (based on chosen model) ---
    double selected_rank_dependence = std::numeric_limits<double>::quiet_NaN();
};

// ============================================================================
// Choice of dependence model used in analysis
// ============================================================================
enum class CopulaDependenceModel
{
    GaussianRank,  // Use Gaussian copula fitted on rank data
    Empirical      // Use empirical copula directly
};

// ============================================================================
// Options controlling copula analysis behavior
// ============================================================================
struct CopulaAnalysisOptions
{
    bool compute_diagnostics = false;     // Whether to compute full diagnostics (costly)

    int bootstrap_B = 50;                 // Number of bootstrap samples for GOF tests
    int max_points = 1200;               // Maximum number of points (downsampling limit)
    int empirical_copula_bins = 20;      // Resolution for binned empirical copula

    CopulaDependenceModel dependence_model = CopulaDependenceModel::GaussianRank;
};

// ============================================================================
// Utility functions
// ============================================================================

// Compute Pearson correlation between two vectors
double pearson_corr(const std::vector<double>& a, const std::vector<double>& b);

// Convert values to fractional ranks in (0,1)
std::vector<double> fractional_ranks_01(const std::vector<double>& x);

// Approximate inverse normal CDF (Phi^{-1})
double norminv_approx(double p);

// Approximate Gaussian copula CDF C(u,v; rho)
double gaussian_copula_cdf_approx(double u, double v, double rho);

// Compute Kendall’s tau (naive O(N^2) implementation)
double kendall_tau_naive(const std::vector<double>& x, const std::vector<double>& y);

// Compute Mardia’s skewness and kurtosis for bivariate data
void mardia_bivariate_moments(const std::vector<double>& x,
                              const std::vector<double>& y,
                              double& skew_out,
                              double& kurt_out);

// Cramér–von Mises statistic for Gaussian copula fit
double gaussian_copula_cvm_statistic(const std::vector<double>& u,
                                     const std::vector<double>& v,
                                     double rho);

// Estimate p-value via bootstrap for Gaussian copula GOF
double estimate_gaussian_copula_gof_pvalue(const std::vector<double>& u,
                                           const std::vector<double>& v,
                                           double rho,
                                           int B,
                                           unsigned long seed,
                                           double& stat_obs_out);

// Downsample vector evenly to limit size
std::vector<double> downsample_evenly(const std::vector<double>& a, int max_points);

// ============================================================================
// Core copula processing functions
// ============================================================================

// Build CopulaPairData from a PathwaySet (extract qx pairs, ranks, Gaussian transform)
CopulaPairData build_copula_pair_data(const PathwaySet& pair_set);

// Build CopulaPairData from generic paired samples stored as (t,c) in a TimeSeries
CopulaPairData build_copula_pair_data(const TimeSeries<double>& pairs);

// Build a binned empirical copula matrix (2D histogram in [0,1]^2)
CopulaBinnedMatrix build_empirical_binned_copula(
    const std::vector<double>& u1,
    const std::vector<double>& u2,
    int nBins);

// Run full copula diagnostics on prepared data
CopulaDiagnostics analyze_copula_pair_data(
    const CopulaPairData& data,
    double delta_x,
    const CopulaAnalysisOptions& opts);

// Convenience function: build data, analyze, and write CSV of rank pairs
CopulaDiagnostics analyze_and_write_rank_pairs(
    const PathwaySet& pair_set,
    double delta_x,
    const std::string& scatter_csv_path,
    const CopulaAnalysisOptions& opts,
    CopulaBinnedMatrix* empirical_binned_out = nullptr);

// Convenience function: analyze generic paired samples and write CSV/VTP
CopulaDiagnostics analyze_and_write_sample_pairs(
    const TimeSeries<double>& pairs,
    double delta_x,
    const std::string& scatter_csv_path,
    const CopulaAnalysisOptions& opts,
    CopulaBinnedMatrix* empirical_binned_out = nullptr);

// ============================================================================
// Output / visualization helpers
// ============================================================================

// Write raw and transformed pair data to CSV
bool write_rank_pairs_csv(const CopulaPairData& data, const std::string& filename);

// Write 2D matrix (e.g., empirical copula density) as VTK ImageData (.vti)
bool write_matrix_as_vti_2d(const CMatrix& M,
                            const std::string& filename,
                            const std::string& array_name = "copula",
                            bool point_data = false);

// Write copula points as VTK PolyData (.vtp) with optional attributes
bool write_rank_points_as_vtp(const std::vector<double>& u1,
                              const std::vector<double>& u2,
                              const std::vector<double>* qx1,
                              const std::vector<double>* qx2,
                              const std::vector<double>* z1,
                              const std::vector<double>* z2,
                              const std::string& filename);

// Overload using CopulaPairData directly
bool write_rank_points_as_vtp(const CopulaPairData& data,
                              const std::string& filename);
