// Focused synthetic tests for clustering.h's k-medoids implementation.
//
// Standalone: only needs sim_runner.h + clustering.h (both PETSc/VTK-free at
// the header level) plus GSL + Armadillo (pulled in transitively by
// Utilities/TimeSeries.h). No MPI, no PETSc, no full simulation pipeline.
//
// Build (from the repository root):
//   g++ -std=c++17 -I. -IUtilities \
//       -I/usr/include/x86_64-linux-gnu/qt6 \
//       tests/test_clustering.cpp -lgsl -lgslcblas -larmadillo \
//       -o /tmp/test_clustering
//   /tmp/test_clustering
//
// Exits 0 and prints "ALL TESTS PASSED" on success; aborts (via assert) on
// the first failing check otherwise.

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include "sim_runner.h"
#include "clustering.h"

namespace {

TimeSeries<double> make_invcdf(double v_lo, double v_hi, int n = 11)
{
    TimeSeries<double> ts;
    for (int k = 0; k < n; ++k) {
        const double u = static_cast<double>(k) / (n - 1);
        ts.addPoint(u, v_lo + (v_hi - v_lo) * u);
    }
    return ts;
}

FineScaleOutputs make_outputs(const std::vector<std::pair<double,double>>& v_ranges,
                              const std::vector<double>& lc,
                              const std::vector<double>& lx,
                              const std::vector<double>& ly)
{
    FineScaleOutputs fo;
    for (size_t i = 0; i < v_ranges.size(); ++i)
        fo.qx_inverse_cdf_all.push_back(make_invcdf(v_ranges[i].first, v_ranges[i].second));
    fo.lc_all = lc;
    fo.velocity_lx_all = lx;
    fo.velocity_ly_all = ly;
    return fo;
}

double weight_sum(const ClusterMetrics& m)
{
    double s = 0.0;
    for (const auto& c : m.clusters) s += c.weight;
    return s;
}

void assert_weights_sum_to_one(const ClusterMetrics& m, double tol = 1e-9)
{
    if (m.clusters.empty()) return;
    const double s = weight_sum(m);
    assert(std::abs(s - 1.0) < tol && "cluster weights must sum to 1 within float tolerance");
}

void assert_no_nan(const ClusterMetrics& m)
{
    assert(std::isfinite(m.avg_cluster_size));
    for (const auto& c : m.clusters) assert(std::isfinite(c.weight));
}

void test_deterministic_assignments()
{
    auto fo = make_outputs(
        {{0.1,1.0},{0.2,1.1},{5.0,6.0},{5.1,6.1},{9.0,10.0},{9.2,10.2}},
        {0.5,0.52,0.6,0.61,0.7,0.72},
        {1.0,1.02,1.1,1.12,1.2,1.22},
        {0.1,0.11,0.12,0.13,0.14,0.15});

    const ClusterMetrics a = KMedoidsClustering::cluster(fo, 3, "velocity_lc");
    const ClusterMetrics b = KMedoidsClustering::cluster(fo, 3, "velocity_lc");

    assert(a.num_clusters == b.num_clusters);
    assert(a.clusters.size() == b.clusters.size());
    for (size_t i = 0; i < a.clusters.size(); ++i) {
        assert(a.clusters[i].realization_id == b.clusters[i].realization_id);
        assert(a.clusters[i].cluster_members == b.clusters[i].cluster_members);
        assert(std::abs(a.clusters[i].weight - b.clusters[i].weight) < 1e-12);
    }
    assert_weights_sum_to_one(a);
    std::cout << "test_deterministic_assignments: OK (" << a.num_clusters << " clusters)\n";
}

void test_k_equals_1()
{
    auto fo = make_outputs(
        {{0.1,1.0},{5.0,6.0},{9.0,10.0},{3.0,4.0}},
        {0.5,0.6,0.7,0.55}, {1.0,1.1,1.2,1.05}, {0.1,0.12,0.14,0.11});

    const ClusterMetrics m = KMedoidsClustering::cluster(fo, 1, "velocity_lc");
    assert(m.num_clusters == 1);
    assert(m.clusters.size() == 1);
    assert(m.clusters[0].cluster_members.size() == 4);
    assert(std::abs(m.clusters[0].weight - 1.0) < 1e-12);
    assert_weights_sum_to_one(m);
    std::cout << "test_k_equals_1: OK\n";
}

void test_k_equals_n()
{
    auto fo = make_outputs(
        {{0.1,1.0},{5.0,6.0},{9.0,10.0},{3.0,4.0},{2.0,3.0}},
        {0.5,0.6,0.7,0.55,0.52}, {1.0,1.1,1.2,1.05,1.02}, {0.1,0.12,0.14,0.11,0.105});
    const int n = 5;

    const ClusterMetrics m = KMedoidsClustering::cluster(fo, n, "velocity_lc");
    assert(m.num_clusters == n);
    assert(m.clusters.size() == static_cast<size_t>(n));
    for (const auto& c : m.clusters) {
        assert(c.cluster_members.size() == 1);
        assert(std::abs(c.weight - 1.0 / n) < 1e-12);
    }
    assert_weights_sum_to_one(m);
    std::cout << "test_k_equals_n: OK\n";
}

void test_separated_groups()
{
    // Two well-separated groups of 4 realizations each: "slow" (v in
    // [0.1,1.1]-ish) and "fast" (v in [19,21]-ish). K=2 must recover both
    // groups intact, regardless of which member becomes each medoid.
    auto fo = make_outputs(
        {{0.10,1.00},{0.11,1.05},{0.12,1.02},{0.09,0.98},
         {19.0,21.0},{19.2,21.3},{18.9,20.8},{19.1,21.1}},
        {0.50,0.51,0.52,0.49, 0.80,0.81,0.79,0.82},
        {1.00,1.01,1.02,0.99, 2.00,2.02,1.98,2.01},
        {0.10,0.11,0.09,0.10, 0.30,0.31,0.29,0.30});

    const ClusterMetrics m = KMedoidsClustering::cluster(fo, 2, "velocity_lc");
    assert(m.num_clusters == 2);
    assert(m.clusters.size() == 2);

    std::vector<int> slow_group = {1, 2, 3, 4};
    std::vector<int> fast_group = {5, 6, 7, 8};
    for (const auto& c : m.clusters) {
        std::vector<int> members = c.cluster_members;
        std::sort(members.begin(), members.end());
        const bool is_slow = (members == slow_group);
        const bool is_fast = (members == fast_group);
        assert((is_slow || is_fast) &&
               "k-medoids must recover the two well-separated groups exactly");
    }
    assert_weights_sum_to_one(m);
    std::cout << "test_separated_groups: OK\n";
}

void test_malformed_inputs()
{
    // One realization has a NaN correlation length, another has an empty
    // inverse-CDF. Neither should crash the clustering or poison the
    // weights/assignments of the well-formed realizations with NaN.
    FineScaleOutputs fo = make_outputs(
        {{0.1,1.0},{0.2,1.1},{5.0,6.0},{5.1,6.1}},
        {0.5, std::numeric_limits<double>::quiet_NaN(), 0.6, 0.61},
        {1.0, 1.02, 1.1, 1.12},
        {0.1, 0.11, 0.12, 0.13});
    fo.qx_inverse_cdf_all[3] = TimeSeries<double>();  // empty

    const ClusterMetrics m = KMedoidsClustering::cluster(fo, 2, "velocity_lc");
    assert_no_nan(m);
    assert_weights_sum_to_one(m);
    // Every one of the 4 realizations must still end up assigned somewhere.
    int total_members = 0;
    for (const auto& c : m.clusters) total_members += (int)c.cluster_members.size();
    assert(total_members == 4);
    std::cout << "test_malformed_inputs: OK (" << m.num_clusters << " clusters)\n";
}

void test_single_and_empty_inputs()
{
    // Zero realizations.
    FineScaleOutputs empty_fo;
    const ClusterMetrics m0 = KMedoidsClustering::cluster(empty_fo, 3, "velocity_lc");
    assert(m0.num_clusters == 0);
    assert(m0.clusters.empty());

    // Requested K <= 0 must not crash.
    auto fo = make_outputs({{0.1,1.0},{5.0,6.0}}, {0.5,0.6}, {1.0,1.1}, {0.1,0.12});
    const ClusterMetrics m1 = KMedoidsClustering::cluster(fo, 0, "velocity_lc");
    assert(m1.num_clusters == 0);
    std::cout << "test_single_and_empty_inputs: OK\n";
}

} // namespace

int main()
{
    test_deterministic_assignments();
    test_k_equals_1();
    test_k_equals_n();
    test_separated_groups();
    test_malformed_inputs();
    test_single_and_empty_inputs();
    std::cout << "ALL TESTS PASSED\n";
    return 0;
}
