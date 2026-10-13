#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "sim_runner.h"

// ============================================================================
// Representative Clustering for Ensemble Reduction
// ============================================================================

struct ClusterInfo {
    int realization_id;           // 1-based realization number
    std::vector<int> cluster_members;  // All realizations in this cluster
    double weight;                // Cluster weight = size / total
};

struct ClusterMetrics {
    int num_clusters;
    int total_realizations;
    double avg_cluster_size;
    std::vector<ClusterInfo> clusters;
};

// ============================================================================
// Distance Metrics
// ============================================================================

// Wasserstein distance between two inverse CDFs
inline double wasserstein_distance(
    const TimeSeries<double>& cdf1,
    const TimeSeries<double>& cdf2)
{
    if (cdf1.size() == 0 || cdf2.size() == 0)
        return std::numeric_limits<double>::max();

    double max_diff = 0.0;
    size_t i = 0, j = 0;

    while (i < cdf1.size() && j < cdf2.size()) {
        double t1 = cdf1.getTime(i);
        double t2 = cdf2.getTime(j);

        if (std::abs(t1 - t2) < 1e-10) {
            max_diff = std::max(max_diff, std::abs(cdf1.getValue(i) - cdf2.getValue(j)));
            ++i; ++j;
        } else if (t1 < t2) {
            ++i;
        } else {
            ++j;
        }
    }
    return max_diff;
}

// Combined distance: velocity CDF + lc + lambda_x/y
inline double combined_distance(
    const TimeSeries<double>& cdf1,
    const TimeSeries<double>& cdf2,
    double lc1, double lc2,
    double lx1, double lx2,
    double ly1, double ly2)
{
    double w_cdf = wasserstein_distance(cdf1, cdf2);
    double w_lc = std::abs(lc1 - lc2) / (std::abs(lc1) + std::abs(lc2) + 1e-10);
    double w_lx = std::abs(lx1 - lx2) / (std::abs(lx1) + std::abs(lx2) + 1e-10);
    double w_ly = std::abs(ly1 - ly2) / (std::abs(ly1) + std::abs(ly2) + 1e-10);

    // Weighted combination: 50% CDF, 50% spatial params
    return 0.5 * w_cdf + 0.1667 * (w_lc + w_lx + w_ly);
}

// ============================================================================
// K-means Clustering
// ============================================================================

class KMeansClustering {
public:
    static ClusterMetrics cluster(
        const FineScaleOutputs& fine_outputs,
        int num_clusters,
        const std::string& metric = "velocity_lc")
    {
        const int n_real = fine_outputs.qx_inverse_cdf_all.size();

        ClusterMetrics result;
        result.total_realizations = n_real;
        result.num_clusters = num_clusters;

        if (n_real <= num_clusters) {
            // Each realization is its own cluster
            for (int i = 0; i < n_real; ++i) {
                ClusterInfo ci;
                ci.realization_id = i + 1;
                ci.cluster_members.push_back(i + 1);
                ci.weight = 1.0 / n_real;
                result.clusters.push_back(ci);
            }
            result.num_clusters = n_real;
            result.avg_cluster_size = 1.0;
            return result;
        }

        // Initialize cluster centers (pick first num_clusters realizations)
        std::vector<int> center_ids;
        for (int i = 0; i < num_clusters; ++i) {
            center_ids.push_back(i);
        }

        std::vector<int> assignments(n_real);
        bool converged = false;
        int max_iterations = 20;
        int iteration = 0;

        while (!converged && iteration < max_iterations) {
            // Assign each realization to nearest cluster
            for (int i = 0; i < n_real; ++i) {
                double min_dist = std::numeric_limits<double>::max();
                int best_cluster = 0;

                for (int c = 0; c < num_clusters; ++c) {
                    int center_id = center_ids[c];
                    double dist = compute_distance(fine_outputs, i, center_id, metric);
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_cluster = c;
                    }
                }
                assignments[i] = best_cluster;
            }

            // Update cluster centers (pick medoid of each cluster)
            std::vector<int> new_center_ids = center_ids;
            converged = true;

            for (int c = 0; c < num_clusters; ++c) {
                std::vector<int> cluster_members;
                for (int i = 0; i < n_real; ++i) {
                    if (assignments[i] == c) {
                        cluster_members.push_back(i);
                    }
                }

                if (!cluster_members.empty()) {
                    // Find medoid (member with min sum of distances to others)
                    int best_medoid = cluster_members[0];
                    double best_sum_dist = std::numeric_limits<double>::max();

                    for (int medoid : cluster_members) {
                        double sum_dist = 0.0;
                        for (int other : cluster_members) {
                            if (medoid != other) {
                                sum_dist += compute_distance(fine_outputs, medoid, other, metric);
                            }
                        }
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
            }

            center_ids = new_center_ids;
            iteration++;
        }

        // Build final clusters
        std::vector<std::vector<int>> final_clusters(num_clusters);
        for (int i = 0; i < n_real; ++i) {
            final_clusters[assignments[i]].push_back(i + 1);
        }

        result.avg_cluster_size = (double)n_real / num_clusters;
        for (int c = 0; c < num_clusters; ++c) {
            if (!final_clusters[c].empty()) {
                ClusterInfo ci;
                ci.realization_id = center_ids[c] + 1;  // Representative
                ci.cluster_members = final_clusters[c];
                ci.weight = (double)final_clusters[c].size() / n_real;
                result.clusters.push_back(ci);
            }
        }

        return result;
    }

private:
    static double compute_distance(
        const FineScaleOutputs& fine_outputs,
        int real_i,
        int real_j,
        const std::string& metric)
    {
        if (real_i >= (int)fine_outputs.qx_inverse_cdf_all.size() ||
            real_j >= (int)fine_outputs.qx_inverse_cdf_all.size()) {
            return std::numeric_limits<double>::max();
        }

        const TimeSeries<double>& cdf_i = fine_outputs.qx_inverse_cdf_all[real_i];
        const TimeSeries<double>& cdf_j = fine_outputs.qx_inverse_cdf_all[real_j];

        if (metric == "velocity_lc") {
            // Velocity CDF + correlation length
            return combined_distance(
                cdf_i, cdf_j,
                fine_outputs.lc_all[real_i],
                fine_outputs.lc_all[real_j],
                fine_outputs.velocity_lx_all[real_i],
                fine_outputs.velocity_lx_all[real_j],
                fine_outputs.velocity_ly_all[real_i],
                fine_outputs.velocity_ly_all[real_j]);
        } else if (metric == "lc_only") {
            // Only correlation length
            return std::abs(fine_outputs.lc_all[real_i] - fine_outputs.lc_all[real_j]) /
                   (std::abs(fine_outputs.lc_all[real_i]) + std::abs(fine_outputs.lc_all[real_j]) + 1e-10);
        } else {
            // Default: velocity CDF only
            return wasserstein_distance(cdf_i, cdf_j);
        }
    }
};

// ============================================================================
// Clustering Output
// ============================================================================

inline void write_clustering_summary(
    const std::string& run_dir,
    const ClusterMetrics& metrics)
{
    const std::string path = run_dir + "clustering_summary.txt";
    std::ofstream f(path);
    if (!f) return;

    f << "# Realization Clustering Summary\n";
    f << "num_clusters=" << metrics.num_clusters << "\n";
    f << "total_realizations=" << metrics.total_realizations << "\n";
    f << "avg_cluster_size=" << metrics.avg_cluster_size << "\n";
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
