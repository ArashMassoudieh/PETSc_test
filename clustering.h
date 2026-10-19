#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>
#include "sim_runner.h"

// ============================================================================
// Representative Clustering for Ensemble Reduction (k-medoids)
//
// Clusters realizations by velocity-CDF / correlation-length similarity and
// selects one representative (medoid) per cluster, so Phase 2 can upscale
// only the representatives and weight-average their results instead of
// upscaling every realization (Phase 1).
//
// Named k-medoids, not k-means: cluster centers are always actual
// realizations (medoids), never a synthetic average of several
// realizations' curves/parameters.
// ============================================================================

struct ClusterInfo {
    int realization_id;                // 1-based representative (medoid)
    std::vector<int> cluster_members;  // All realizations in this cluster (1-based)
    double weight;                     // Cluster weight = size / total
};

struct ClusterMetrics {
    int num_clusters = 0;              // ACTUAL number of nonempty clusters produced
    int requested_clusters = 0;        // K as requested by the caller
    int total_realizations = 0;
    double avg_cluster_size = 0.0;
    std::vector<ClusterInfo> clusters;

    // Diagnostics (Task 4)
    std::string metric;
    int iterations_run = 0;
    bool converged = false;
};

// ============================================================================
// Distance Metrics
// ============================================================================

// True (approximate) 1-Wasserstein distance between two velocity inverse-CDFs,
// i.e. quantile functions F1^{-1}(u), F2^{-1}(u) with u in [0,1]:
//
//   W1(F1, F2) = integral_0^1 |F1^{-1}(u) - F2^{-1}(u)| du
//
// Approximated by evaluating both quantile functions (via TimeSeries::interpol,
// which is robust to each series having its own native u-grid) on a shared
// uniform grid and integrating the absolute difference with the trapezoidal
// rule. This replaces the previous implementation, which was mislabeled: it
// only compared points whose native times matched almost exactly (silently
// dropping the rest when grids differed) and took their max absolute
// difference (an L-infinity distance), not an integral.
inline double wasserstein_distance(
    const TimeSeries<double>& invcdf1,
    const TimeSeries<double>& invcdf2,
    int n_quantile_points = 101)
{
    if (invcdf1.size() == 0 || invcdf2.size() == 0)
        return std::numeric_limits<double>::max();
    if (n_quantile_points < 2) n_quantile_points = 2;

    double integral = 0.0;
    double v_prev = std::numeric_limits<double>::quiet_NaN();
    double u_prev = 0.0;

    for (int k = 0; k < n_quantile_points; ++k) {
        const double u = static_cast<double>(k) / (n_quantile_points - 1);
        const double v1 = invcdf1.interpol(u);
        const double v2 = invcdf2.interpol(u);

        if (!std::isfinite(v1) || !std::isfinite(v2)) {
            // A non-finite quantile value anywhere makes the comparison
            // meaningless; report maximal (but finite) dissimilarity instead
            // of propagating NaN into every downstream distance comparison.
            return std::numeric_limits<double>::max() / 4.0;
        }

        const double v = std::abs(v1 - v2);
        if (k > 0) {
            integral += 0.5 * (v + v_prev) * (u - u_prev);
        }
        v_prev = v;
        u_prev = u;
    }
    return integral;
}

// Combined distance: velocity CDF (Wasserstein) + lc + lambda_x/y (relative
// differences). Falls back to a large-but-finite distance if any scalar
// feature is non-finite, instead of letting NaN propagate into comparisons.
inline double combined_distance(
    const TimeSeries<double>& cdf1,
    const TimeSeries<double>& cdf2,
    double lc1, double lc2,
    double lx1, double lx2,
    double ly1, double ly2)
{
    const double sentinel = std::numeric_limits<double>::max() / 4.0;

    auto relative_diff = [sentinel](double a, double b) -> double {
        if (!std::isfinite(a) || !std::isfinite(b)) return sentinel;
        return std::abs(a - b) / (std::abs(a) + std::abs(b) + 1e-10);
    };

    const double w_cdf = wasserstein_distance(cdf1, cdf2);
    const double w_lc  = relative_diff(lc1, lc2);
    const double w_lx  = relative_diff(lx1, lx2);
    const double w_ly  = relative_diff(ly1, ly2);

    // Weighted combination: 50% CDF, 50% spatial params (lc, lx, ly equally).
    return 0.5 * w_cdf + 0.1667 * (w_lc + w_lx + w_ly);
}

// ============================================================================
// K-medoids Clustering
// ============================================================================

class KMedoidsClustering {
public:
    static ClusterMetrics cluster(
        const FineScaleOutputs& fine_outputs,
        int num_clusters,
        const std::string& metric = "velocity_lc")
    {
        const int n_real = static_cast<int>(fine_outputs.qx_inverse_cdf_all.size());

        ClusterMetrics result;
        result.total_realizations = n_real;
        result.requested_clusters = num_clusters;
        result.metric = metric;

        if (n_real <= 0 || num_clusters <= 0) {
            result.num_clusters = 0;
            return result;
        }

        if (n_real <= num_clusters) {
            // Each realization is its own cluster -- nothing to reduce.
            for (int i = 0; i < n_real; ++i) {
                ClusterInfo ci;
                ci.realization_id = i + 1;
                ci.cluster_members.push_back(i + 1);
                ci.weight = 1.0 / n_real;
                result.clusters.push_back(ci);
            }
            result.num_clusters = n_real;
            result.avg_cluster_size = 1.0;
            result.converged = true;
            result.iterations_run = 0;
            return result;
        }

        // ---- Deterministic, spread-out initialization (farthest-first) ----
        // Start from realization 0 (deterministic, no RNG), then repeatedly
        // add the realization whose distance to its nearest already-chosen
        // center is largest. This spreads the initial medoids across the
        // population instead of clumping them in the first K realizations
        // (which, for realizations generated in sequence, are not a random
        // sample of the whole set).
        std::vector<int> center_ids;
        center_ids.push_back(0);
        while (static_cast<int>(center_ids.size()) < num_clusters) {
            int farthest_point = -1;
            double farthest_dist = -1.0;
            for (int i = 0; i < n_real; ++i) {
                if (std::find(center_ids.begin(), center_ids.end(), i) != center_ids.end())
                    continue;
                double min_dist_to_center = std::numeric_limits<double>::max();
                for (int c : center_ids) {
                    const double d = compute_distance(fine_outputs, i, c, metric);
                    min_dist_to_center = std::min(min_dist_to_center, d);
                }
                if (min_dist_to_center > farthest_dist) {
                    farthest_dist = min_dist_to_center;
                    farthest_point = i;
                }
            }
            if (farthest_point < 0) break;  // fewer distinct points than requested K
            center_ids.push_back(farthest_point);
        }

        std::vector<int> assignments(n_real, 0);
        auto assign_all = [&](const std::vector<int>& centers) {
            for (int i = 0; i < n_real; ++i) {
                double min_dist = std::numeric_limits<double>::max();
                int best_cluster = 0;
                for (int c = 0; c < (int)centers.size(); ++c) {
                    const double dist = compute_distance(fine_outputs, i, centers[c], metric);
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_cluster = c;
                    }
                }
                assignments[i] = best_cluster;
            }
        };

        bool converged = false;
        const int max_iterations = 20;
        int iteration = 0;

        for (; iteration < max_iterations; ++iteration) {
            assign_all(center_ids);

            // Update cluster centers (pick medoid of each cluster).
            std::vector<int> new_center_ids = center_ids;
            converged = true;

            for (int c = 0; c < (int)center_ids.size(); ++c) {
                std::vector<int> members;
                for (int i = 0; i < n_real; ++i)
                    if (assignments[i] == c) members.push_back(i);

                if (members.empty()) continue;  // handled after reassignment below

                int best_medoid = members[0];
                double best_sum_dist = std::numeric_limits<double>::max();
                for (int medoid : members) {
                    double sum_dist = 0.0;
                    for (int other : members)
                        if (medoid != other) sum_dist += compute_distance(fine_outputs, medoid, other, metric);
                    if (sum_dist < best_sum_dist) {
                        best_sum_dist = sum_dist;
                        best_medoid = medoid;
                    }
                }

                if (new_center_ids[c] != best_medoid) {
                    new_center_ids[c] = best_medoid;
                    converged = false;
                }
            }

            center_ids = new_center_ids;
            if (converged) { ++iteration; break; }
        }

        // Recompute final assignments against the FINAL centers. Without this,
        // `assignments` reflects the centers from the start of the last loop
        // iteration, which can be stale relative to `center_ids` whenever the
        // medoids changed on that last pass (in particular whenever the loop
        // exits by hitting max_iterations rather than by converging).
        assign_all(center_ids);

        result.converged = converged;
        result.iterations_run = iteration;

        // Build final clusters, dropping any center that ended up with no
        // members (can happen with farthest-point initialization when two
        // chosen centers turn out to be mutually closer to a third point than
        // to each other after reassignment).
        std::vector<std::vector<int>> final_clusters(center_ids.size());
        for (int i = 0; i < n_real; ++i)
            final_clusters[assignments[i]].push_back(i + 1);

        for (int c = 0; c < (int)center_ids.size(); ++c) {
            if (final_clusters[c].empty()) continue;
            ClusterInfo ci;
            ci.realization_id = center_ids[c] + 1;
            ci.cluster_members = final_clusters[c];
            ci.weight = static_cast<double>(final_clusters[c].size()) / n_real;
            result.clusters.push_back(ci);
        }

        // The reported cluster count is the number of NONEMPTY clusters
        // actually produced, which can be less than the requested K.
        result.num_clusters = static_cast<int>(result.clusters.size());
        result.avg_cluster_size = result.num_clusters > 0
            ? static_cast<double>(n_real) / result.num_clusters
            : 0.0;

        return result;
    }

private:
    static double compute_distance(
        const FineScaleOutputs& fine_outputs,
        int real_i,
        int real_j,
        const std::string& metric)
    {
        if (real_i == real_j) return 0.0;
        if (real_i < 0 || real_j < 0 ||
            real_i >= (int)fine_outputs.qx_inverse_cdf_all.size() ||
            real_j >= (int)fine_outputs.qx_inverse_cdf_all.size()) {
            return std::numeric_limits<double>::max() / 4.0;
        }

        const TimeSeries<double>& cdf_i = fine_outputs.qx_inverse_cdf_all[real_i];
        const TimeSeries<double>& cdf_j = fine_outputs.qx_inverse_cdf_all[real_j];

        if (metric == "velocity_lc") {
            return combined_distance(
                cdf_i, cdf_j,
                fine_outputs.lc_all[real_i], fine_outputs.lc_all[real_j],
                fine_outputs.velocity_lx_all[real_i], fine_outputs.velocity_lx_all[real_j],
                fine_outputs.velocity_ly_all[real_i], fine_outputs.velocity_ly_all[real_j]);
        } else if (metric == "lc_only") {
            const double a = fine_outputs.lc_all[real_i];
            const double b = fine_outputs.lc_all[real_j];
            if (!std::isfinite(a) || !std::isfinite(b))
                return std::numeric_limits<double>::max() / 4.0;
            return std::abs(a - b) / (std::abs(a) + std::abs(b) + 1e-10);
        } else {
            // "velocity" (default): velocity CDF only.
            return wasserstein_distance(cdf_i, cdf_j);
        }
    }
};

// Backward-compatible alias for the one call site that predates the rename.
using KMeansClustering = KMedoidsClustering;

// ============================================================================
// Clustering Output (Task 4: extended diagnostics)
// ============================================================================

inline void write_clustering_summary(
    const std::string& run_dir,
    const ClusterMetrics& metrics,
    double runtime_seconds = -1.0)
{
    const std::string path = run_dir + "clustering_summary.txt";
    std::ofstream f(path);
    if (!f) return;

    double weight_sum = 0.0;
    for (const auto& c : metrics.clusters) weight_sum += c.weight;

    f << "# Realization Clustering Summary\n";
    f << "metric=" << metrics.metric << "\n";
    f << "requested_clusters=" << metrics.requested_clusters << "\n";
    f << "num_clusters=" << metrics.num_clusters << " (actual, nonempty)\n";
    f << "total_realizations=" << metrics.total_realizations << "\n";
    f << "avg_cluster_size=" << metrics.avg_cluster_size << "\n";
    f << "iterations_run=" << metrics.iterations_run << "\n";
    f << "converged=" << (metrics.converged ? "true" : "false") << "\n";
    f << "weight_sum=" << weight_sum
      << " (should be 1.0 within floating-point tolerance)\n";
    if (runtime_seconds >= 0.0)
        f << "phase2_runtime_seconds=" << runtime_seconds << "\n";

    f << "\n# Cluster Assignments (representative_id, weight, members)\n";
    for (const auto& cluster : metrics.clusters) {
        f << "cluster_representative=" << cluster.realization_id << "\n";
        f << "cluster_weight=" << cluster.weight << "\n";
        f << "cluster_size=" << cluster.cluster_members.size() << "\n";
        f << "members=";
        for (size_t i = 0; i < cluster.cluster_members.size(); ++i) {
            if (i > 0) f << ",";
            f << cluster.cluster_members[i];
        }
        f << "\n\n";
    }
}
