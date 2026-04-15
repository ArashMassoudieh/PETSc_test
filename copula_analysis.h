#pragma once

#include <string>
#include <vector>
#include <limits>

#include "Matrix.h"
#include "PathwaySet.h"

struct CopulaPairData
{
    std::vector<double> q1;
    std::vector<double> q2;
    std::vector<double> u1;
    std::vector<double> u2;
    std::vector<double> z1;
    std::vector<double> z2;

    int n_pairs = 0;
};

struct CopulaDiagnostics
{
    double delta_x = 0.0;
    int n_pairs = 0;
    double corr_qx = std::numeric_limits<double>::quiet_NaN();
    double corr_rank = std::numeric_limits<double>::quiet_NaN();

    double gaussian_copula_rho = std::numeric_limits<double>::quiet_NaN();
    double gaussian_copula_gof_stat = std::numeric_limits<double>::quiet_NaN();
    double gaussian_copula_gof_pvalue = std::numeric_limits<double>::quiet_NaN();

    double kendall_tau = std::numeric_limits<double>::quiet_NaN();
    double rho_from_tau = std::numeric_limits<double>::quiet_NaN();
    double mardia_skewness = std::numeric_limits<double>::quiet_NaN();
    double mardia_kurtosis = std::numeric_limits<double>::quiet_NaN();

    double empirical_copula_stat = std::numeric_limits<double>::quiet_NaN();
    double empirical_copula_pvalue = std::numeric_limits<double>::quiet_NaN();
    double empirical_upper_tail_frac_90 = std::numeric_limits<double>::quiet_NaN();
    double empirical_diagonal_l1 = std::numeric_limits<double>::quiet_NaN();

    double selected_rank_dependence = std::numeric_limits<double>::quiet_NaN();
};

enum class CopulaDependenceModel
{
    GaussianRank,
    Empirical
};

struct CopulaAnalysisOptions
{
    bool compute_diagnostics = false;
    int bootstrap_B = 50;
    int max_points = 1200;
    int empirical_copula_bins = 20;
    CopulaDependenceModel dependence_model = CopulaDependenceModel::GaussianRank;
};

double pearson_corr(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> fractional_ranks_01(const std::vector<double>& x);
double norminv_approx(double p);
double gaussian_copula_cdf_approx(double u, double v, double rho);
double kendall_tau_naive(const std::vector<double>& x, const std::vector<double>& y);
void mardia_bivariate_moments(const std::vector<double>& x,
                              const std::vector<double>& y,
                              double& skew_out,
                              double& kurt_out);
double gaussian_copula_cvm_statistic(const std::vector<double>& u,
                                     const std::vector<double>& v,
                                     double rho);
double estimate_gaussian_copula_gof_pvalue(const std::vector<double>& u,
                                           const std::vector<double>& v,
                                           double rho,
                                           int B,
                                           unsigned long seed,
                                           double& stat_obs_out);
std::vector<double> downsample_evenly(const std::vector<double>& a, int max_points);

CopulaPairData build_copula_pair_data(const PathwaySet& pair_set);

CopulaBinnedMatrix build_empirical_binned_copula(
    const std::vector<double>& u1,
    const std::vector<double>& u2,
    int nBins);

CopulaDiagnostics analyze_copula_pair_data(
    const CopulaPairData& data,
    double delta_x,
    const CopulaAnalysisOptions& opts);

CopulaDiagnostics analyze_and_write_rank_pairs(
    const PathwaySet& pair_set,
    double delta_x,
    const std::string& scatter_csv_path,
    const CopulaAnalysisOptions& opts,
    CopulaBinnedMatrix* empirical_binned_out = nullptr);

bool write_rank_pairs_csv(const CopulaPairData& data, const std::string& filename);

bool write_matrix_as_vti_2d(const CMatrix& M,
                            const std::string& filename,
                            const std::string& array_name = "copula",
                            bool point_data = false);

bool write_rank_points_as_vtp(const std::vector<double>& u1,
                              const std::vector<double>& u2,
                              const std::vector<double>* qx1,
                              const std::vector<double>* qx2,
                              const std::vector<double>* z1,
                              const std::vector<double>* z2,
                              const std::string& filename);

bool write_rank_points_as_vtp(const CopulaPairData& data,
                              const std::string& filename);
