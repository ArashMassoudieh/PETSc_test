// main.cpp
#include "petsc_init.h"
#include "petscmatrix.h"
#include "petscvector.h"
#include <cmath>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <cstdlib>
#include "grid.h"
#include "TimeSeries.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline PetscInt idx(PetscInt i, PetscInt j, PetscInt nx) { return j*nx + i; }

static inline bool on_bc(PetscInt i, PetscInt j, PetscInt nx, PetscInt ny,
                         double x, double y, double Lx, double Ly) {
    const double eps = 1e-14;
    return (i==0) || (j==0) || (std::abs(x - Lx) < eps) || (std::abs(y - Ly) < eps);
}

// Helper: directory creation
bool createDirectory(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        if (info.st_mode & S_IFDIR) return true;
        std::cerr << "Error: Path exists but is not a directory: " << path << "\n";
        return false;
    }
#ifdef _WIN32
    if (mkdir(path.c_str()) != 0) {
#else
    if (mkdir(path.c_str(), 0755) != 0) {
#endif
        if (errno != EEXIST) {
            std::cerr << "Error creating directory " << path << ": " << strerror(errno) << "\n";
            return false;
        }
    }
    std::cout << "Created output directory: " << path << "\n";
    return true;
}

std::string joinPath(const std::string& dir, const std::string& filename) {
    if (dir.empty()) return filename;
    if (dir.back() == '/' || dir.back() == '\\') return dir + filename;
    return dir + "/" + filename;
}

// Timestamp: YYYYMMDD_HHMMSS
static std::string makeTimestamp() {
    std::time_t now = std::time(nullptr);
    std::tm tm_now;
    localtime_r(&now, &tm_now);
    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

// Realization label starting at 1: r001, r002, ...
static std::string makeRealLabel(int r1) {
    std::ostringstream oss;
    oss << "r" << std::setw(4) << std::setfill('0') << r1;
    return oss.str();
}

// For folder name: fine_r001, ...
static std::string makeFineFolder(int r1) {
    return "fine_" + makeRealLabel(r1);
}

static double mean_of(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / (double)v.size();
}

// ===============================
// CSV utilities for aggregation
// ===============================

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') { in_quotes = !in_quotes; continue; }
        if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

static bool try_parse_double(const std::string& s, double& v) {
    char* endp = nullptr;
    v = std::strtod(s.c_str(), &endp);
    return endp != s.c_str() && *endp == '\0';
}

// Reads CSV with header. Assumes first column is time (numeric).
// Returns: times, column_names (excluding time), columns data (each column is vector<double>)
static bool read_time_series_table_csv(
    const std::string& path,
    std::vector<double>& times,
    std::vector<std::string>& colnames,
    std::vector<std::vector<double>>& cols
) {
    std::ifstream f(path);
    if (!f) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    auto h = split_csv_line(header);
    if (h.size() < 2) return false;

    // header[0] is time label (e.g., t or time)
    colnames.assign(h.begin() + 1, h.end());
    cols.assign(colnames.size(), std::vector<double>{});
    times.clear();

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto parts = split_csv_line(line);
        if (parts.size() < 1) continue;

        double t;
        if (!try_parse_double(parts[0], t)) {
            // skip non-numeric rows
            continue;
        }
        times.push_back(t);

        // fill data columns if present, else NaN
        for (size_t j = 0; j < colnames.size(); ++j) {
            double val = std::numeric_limits<double>::quiet_NaN();
            if (j + 1 < parts.size()) {
                double tmp;
                if (try_parse_double(parts[j + 1], tmp)) val = tmp;
            }
            cols[j].push_back(val);
        }
    }

    // sanity: all columns same length
    for (auto& c : cols) {
        if (c.size() != times.size()) return false;
    }
    return true;
}

// Writes aggregated comparison CSV: time + many columns.
static bool write_comparison_csv(
    const std::string& out_path,
    const std::vector<double>& base_time,
    const std::vector<std::string>& out_colnames,
    const std::vector<std::vector<double>>& out_cols
) {
    if (out_colnames.size() != out_cols.size()) return false;
    for (auto& c : out_cols) if (c.size() != base_time.size()) return false;

    std::ofstream f(out_path);
    if (!f) return false;

    f << "t";
    for (auto& n : out_colnames) f << "," << n;
    f << "\n";

    for (size_t i = 0; i < base_time.size(); ++i) {
        f << std::setprecision(12) << base_time[i];
        for (size_t j = 0; j < out_cols.size(); ++j) {
            f << "," << std::setprecision(12) << out_cols[j][i];
        }
        f << "\n";
    }
    return true;
}

// ================================================================
// NEW: Resampling utilities (fixed uniform time grid for aggregation)
// ================================================================
static double lerp(double x0, double y0, double x1, double y1, double x) {
    if (x1 == x0) return y0;
    const double a = (x - x0) / (x1 - x0);
    return y0 + a * (y1 - y0);
}

// Linear resample y(t_src) onto t_dst. Out-of-range clamped to endpoints.
static std::vector<double> resample_col_linear(
    const std::vector<double>& t_src,
    const std::vector<double>& y_src,
    const std::vector<double>& t_dst
) {
    std::vector<double> y_dst;
    y_dst.reserve(t_dst.size());

    if (t_src.empty() || y_src.size() != t_src.size()) {
        y_dst.assign(t_dst.size(), std::numeric_limits<double>::quiet_NaN());
        return y_dst;
    }

    size_t k = 0;
    for (double td : t_dst) {
        if (td <= t_src.front()) { y_dst.push_back(y_src.front()); continue; }
        if (td >= t_src.back())  { y_dst.push_back(y_src.back());  continue; }

        while (k + 1 < t_src.size() && t_src[k + 1] < td) ++k;

        y_dst.push_back(lerp(t_src[k], y_src[k], t_src[k+1], y_src[k+1], td));
    }
    return y_dst;
}

static void resample_table_linear(
    const std::vector<double>& t_src,
    const std::vector<std::vector<double>>& cols_src,
    const std::vector<double>& t_dst,
    std::vector<std::vector<double>>& cols_dst
) {
    cols_dst.clear();
    cols_dst.reserve(cols_src.size());
    for (const auto& c : cols_src) {
        cols_dst.push_back(resample_col_linear(t_src, c, t_dst));
    }
}
// ================================================================

int main(int argc, char** argv) {
    PETScInit petsc(argc, argv);

#ifdef Arash
    std::string output_dir = "/home/arash/Projects/UpscalingResults";
#elif PowerEdge
    std::string output_dir = "/mnt/3rd900/Projects/PETSc_test/Results";
#elif Behzad
    std::string output_dir = "/home/behzad/Projects/PETSc_test/Results";
#elif SligoCreek
    std::string output_dir = "/media/arash/E/Projects/PETSc_test/Results";
#else
    std::string output_dir = "./Results";
#endif

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::string run_dir;
    if (rank == 0) {
        createDirectory(output_dir);
        run_dir = joinPath(output_dir, "run_" + makeTimestamp());
        createDirectory(run_dir);
        std::cout << "Run directory: " << run_dir << "\n";
    }

    // broadcast run_dir
    int len = 0;
    if (rank == 0) len = (int)run_dir.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    run_dir.resize(len);
    MPI_Bcast(run_dir.data(), len, MPI_CHAR, 0, PETSC_COMM_WORLD);

    if (!run_dir.empty() && run_dir.back() != '/' && run_dir.back() != '\\') run_dir += "/";
    MPI_Barrier(PETSC_COMM_WORLD);

    // -----------------------------
    // Domain / grid resolution
    // -----------------------------
    int nx = 300;
    int nu = 100;
    int ny = 100;
    double Lx = 3.0;
    double Ly = 1.0;
    double correlation_ls_x = 1;
    double correlation_ls_y = 0.1;
    double stdev = 1.0;
    double g_mean = 0;
    double Diffusion_coefficient = 0;

    // -----------------------------
    // Monte Carlo / realizations
    // -----------------------------
    const int nReal = 20;     // number of fine realizations
    const double du = 0.01;   // u-grid for averaging inverse CDF
    const int nU = (int)std::round(1.0 / du) + 1;

    std::vector<double> lc_all, lx_all, ly_all, dt_all;
    std::vector<double> invcdf_sum(nU, 0.0);

    // file lists for end aggregation (BTC + derivative)
    std::vector<std::string> fine_btc_files;
    std::vector<std::string> fine_btc_deriv_files;

    std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,dt_opt\n";
    }

    // =====================================================================
    // FINE-SCALE LOOP
    // =====================================================================
    for (int r = 1; r <= nReal; ++r) {
        const std::string rlab = makeRealLabel(r);

        // folder fine_r###/
        std::string fine_dir = joinPath(run_dir, makeFineFolder(r));
        if (rank == 0) createDirectory(fine_dir);
        MPI_Barrier(PETSC_COMM_WORLD);
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        // NEW: prefix all files with r###_
        const std::string pfx = rlab + "_";

        Grid2D g(nx, ny, Lx, Ly);

        PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

        unsigned long run_seed = 20260115UL;  // pick any constant, or read from argv
        unsigned long seed = run_seed + 1000UL*(unsigned long)r + (unsigned long)rank;
        g.makeGaussianFieldSGS("K_normal_score", correlation_ls_x, correlation_ls_y, 10, seed);

        PetscTime(&t_asm0);
        PetscTime(&t_total0);

        g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
        g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "NormalScore.vti"));

        g.createExponentialField("K_normal_score", stdev, g_mean, "K");
        PetscTime(&t_asm1);

        PetscTime(&t_solve0);
        g.DarcySolve(Lx, 0, "K", "K");
        std::cout << "Darcy solved ... " << std::endl;
        g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));
        PetscTime(&t_solve1);

        g.computeMassBalanceError("MassBalanceError");
        g.writeNamedMatrix("MassBalanceError", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "error.txt"));

        g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.vti"));
        g.writeNamedMatrix("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.txt"));
        g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.vti"));
        g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.vti"));

        TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
        TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
        g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);

        g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx_normal_score.vti"));

        std::cout << "Sampling points for derivative ..." << std::endl;
        TimeSeries<double> curvture = g.sampleSecondDerivative("qx_normal_score", Grid2D::ArrayKind::Fx,
                                                               Grid2D::DerivDir::X, 10000, 0.05);
        curvture.writefile(joinPath(fine_dir, pfx + "2nd_deriv.txt"));

        // velocity autocorrelation (X,Y) via perturbation
        double delta_min = 0.001, delta_max = 0.2;
        int num_deltas = 30;
        int num_samples_per_delta = 10000;

        TimeSeries<double> corr_x, corr_y;

        for (int i = 0; i < num_deltas; ++i) {
            double exponent = static_cast<double>(i) / (num_deltas - 1);
            double delta = delta_min * std::pow(delta_max / delta_min, exponent);
            try {
                TimeSeries<double> samples = g.sampleGaussianPerturbation(
                    "qx_normal_score", Grid2D::ArrayKind::Fx,
                    num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                corr_x.append(delta, samples.correlation_tc());
            } catch (...) {}
        }
        corr_x.writefile(joinPath(fine_dir, pfx + "velocity_correlation_x.txt"));
        double lambda_x_emp = corr_x.fitExponentialDecay();

        for (int i = 0; i < num_deltas; ++i) {
            double exponent = static_cast<double>(i) / (num_deltas - 1);
            double delta = delta_min * std::pow(delta_max / delta_min, exponent);
            try {
                TimeSeries<double> samples = g.sampleGaussianPerturbation(
                    "qx_normal_score", Grid2D::ArrayKind::Fx,
                    num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                corr_y.append(delta, samples.correlation_tc());
            } catch (...) {}
        }
        corr_y.writefile(joinPath(fine_dir, pfx + "velocity_correlation_y.txt"));
        double lambda_y_emp = corr_y.fitExponentialDecay();

        // dt
        double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;

        // transport
        g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

        g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.txt"));
        g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.txt"));
        g.writeNamedMatrix("K",  Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.txt"));

        g.SetVal("diffusion", Diffusion_coefficient);
        g.SetVal("porosity", 1);
        g.SetVal("c_left", 1.0);

        std::vector<double> xLocations{0.5, 1.5, 2.5};
        g.setBTCLocations(xLocations);

        TimeSeriesSet<double> BTCs_FineScaled;
        g.SolveTransport(10,
                         std::min(dt_optimal, 0.5 / 10.0),
                         "transport_", 50,
                         fine_dir,
                         "C",
                         &BTCs_FineScaled,
                         r); // <-- realization

        g.writeNamedVTI_Auto("C", joinPath(fine_dir, pfx + "C.vti"));

        // BTC file names with realization prefix
        const std::string btc_path       = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
        const std::string btc_deriv_path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

        BTCs_FineScaled.write(btc_path);
        BTCs_FineScaled.derivative().write(btc_deriv_path);

        // keep list for final aggregation (rank0 will use it)
        if (rank == 0) {
            fine_btc_files.push_back(btc_path);
            fine_btc_deriv_files.push_back(btc_deriv_path);
        }

        PetscTime(&t_total1);

        // pathway correlation -> lc
        PathwaySet pathways;
        pathways.Initialize(1000, PathwaySet::Weighting::FluxWeighted, &g);
        pathways.trackAllPathways(&g, 0.01);

        double Delta_x_min = 0.001, Delta_x_max = 0.5;
        int num_Delta_x = 30;
        int num_samples_per_Delta_x = 10000;

        TimeSeries<double> qx_correlation;
        for (int i = 0; i < num_Delta_x; ++i) {
            double exponent = static_cast<double>(i) / (num_Delta_x - 1);
            double Delta_x = Delta_x_min * std::pow(Delta_x_max / Delta_x_min, exponent);
            try {
                PathwaySet particle_pairs = pathways.sampleParticlePairs(Delta_x, num_samples_per_Delta_x);
                double correlation = particle_pairs.calculateCorrelation(0, 1, "qx");
                qx_correlation.append(Delta_x, correlation);
            } catch (...) {}
        }

        qx_correlation.writefile(joinPath(fine_dir, pfx + "qx_correlation_vs_distance.txt"));
        double advection_correlation_length_scale = qx_correlation.fitExponentialDecay();

        pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
        pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"));

        // inverse CDF
        TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100);
        qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
        qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));

        // meta (configuration)
        if (rank == 0) {
            std::ofstream m(joinPath(fine_dir, pfx + "meta.txt"));
            m << "realization=" << r << "\n";
            m << "nx=" << nx << "\nny=" << ny << "\nnu=" << nu << "\n";
            m << "Lx=" << Lx << "\nLy=" << Ly << "\n";
            m << "correlation_ls_x=" << correlation_ls_x << "\n";
            m << "correlation_ls_y=" << correlation_ls_y << "\n";
            m << "stdev=" << stdev << "\n";
            m << "g_mean=" << g_mean << "\n";
            m << "lc=" << advection_correlation_length_scale << "\n";
            m << "lambda_x=" << lambda_x_emp << "\n";
            m << "lambda_y=" << lambda_y_emp << "\n";
            m << "dt_opt=" << dt_optimal << "\n";
            m << "seed=" << seed << "\n";
        }

        // accumulate means + save summary line
        if (rank == 0) {
            lc_all.push_back(advection_correlation_length_scale);
            lx_all.push_back(lambda_x_emp);
            ly_all.push_back(lambda_y_emp);
            dt_all.push_back(dt_optimal);

            for (int k = 0; k < nU; ++k) {
                double u = k * du;
                invcdf_sum[k] += qx_inverse_cdf.interpol(u);
            }

            std::ofstream f(stats_csv, std::ios::app);
            f << r << "," << advection_correlation_length_scale << ","
              << lambda_x_emp << "," << lambda_y_emp << "," << dt_optimal << "\n";
        }

        if (rank == 0) {
            std::cout << "[Fine " << rlab << "] Assembly time: " << (t_asm1 - t_asm0) << " s\n";
            std::cout << "[Fine " << rlab << "] Solve time:    " << (t_solve1 - t_solve0) << " s\n";
            std::cout << "[Fine " << rlab << "] Total time:   " << (t_total1 - t_total0) << " s\n";
            std::cout << "[Fine " << rlab << "] Outputs saved to: " << fine_dir << "\n";
        }

        MPI_Barrier(PETSC_COMM_WORLD);
    } // end fine loop

    // =====================================================================
    // MEAN PARAMS + UPSCALED RUN
    // =====================================================================
    double lc_mean = 0, lx_mean = 0, ly_mean = 0, dt_mean = 0;
    std::vector<double> invcdf_mean;

    if (rank == 0) {
        lc_mean = mean_of(lc_all);
        lx_mean = mean_of(lx_all);
        ly_mean = mean_of(ly_all);
        dt_mean = mean_of(dt_all);

        invcdf_mean.resize(nU, 0.0);
        for (int k = 0; k < nU; ++k) invcdf_mean[k] = invcdf_sum[k] / (double)nReal;

        {
            std::ofstream f(joinPath(run_dir, "mean_params.txt"));
            f << "nReal=" << nReal << "\n";
            f << "lc_mean=" << lc_mean << "\n";
            f << "lambda_x_mean=" << lx_mean << "\n";
            f << "lambda_y_mean=" << ly_mean << "\n";
            f << "dt_mean=" << dt_mean << "\n";
            f << "du=" << du << "\n";
        }

        {
            std::ofstream f(joinPath(run_dir, "mean_qx_inverse_cdf.txt"));
            f << "u,v\n";
            for (int k = 0; k < nU; ++k) {
                double u = k * du;
                f << u << "," << invcdf_mean[k] << "\n";
            }
        }
    }

    // broadcast scalar means
    MPI_Bcast(&lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // broadcast mean inverse cdf vector
    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // build mean inverse CDF TimeSeries
    TimeSeries<double> qx_inverse_cdf_mean;
    for (int k = 0; k < nU; ++k) {
        double u = k * du;
        qx_inverse_cdf_mean.append(u, invcdf_mean[k]);
    }

    // Upscaled folder + prefix
    std::string up_dir = joinPath(run_dir, "upscaled_mean");
    if (rank == 0) createDirectory(up_dir);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (!up_dir.empty() && up_dir.back() != '/' && up_dir.back() != '\\') up_dir += "/";

    const std::string up_pfx = std::string("upscaled_"); // file prefix for upscaled outputs

    if (rank == 0) {
        std::cout << "\n=== UPSCALED RUN USING MEAN PARAMETERS ===\n";
        std::cout << "lc_mean=" << lc_mean << "\n";
        std::cout << "lambda_x_mean=" << lx_mean << "\n";
        std::cout << "lambda_y_mean=" << ly_mean << "\n";
        std::cout << "dt_mean=" << dt_mean << "\n";
        std::cout << "upscaled output: " << up_dir << "\n";
    }

    // Mixing PDF grid
    Grid2D g_u(nx, ny, Lx, Ly);

    int qx_size = nu * (nx + 1);
    int qy_size = (nu + 1) * nx;

    auto& qx_u = g_u.flux("qx");
    auto& qy_u = g_u.flux("qy");

    qx_u.resize(qx_size, 0.0);
    qy_u.resize(qy_size, 0.0);

    // Assign v(u) from mean inverse CDF
    for (int j = 0; j < nu; ++j) {
        double u = static_cast<double>(j) / (nu - 1);
        double v_at_u = qx_inverse_cdf_mean.interpol(u);
        for (int i = 0; i < nx + 1; ++i) {
            int id = j * (nx + 1) + i;
            if (id >= qx_size) {
                std::cerr << "ERROR: qx index out of bounds: " << id << " >= " << qx_size << "\n";
                return 1;
            }
            qx_u[id] = v_at_u;
        }
    }

    // Save initial fields with prefix
    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.txt"));

    // Mixing parameters (mean)
    g_u.setMixingParams(lc_mean, lx_mean, ly_mean);

    g_u.SetVal("diffusion", Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    double t_end_pdf = 10;
    double dt_pdf = dt_mean;
    int output_interval_pdf = 50;

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

    std::vector<double> xLocations{0.5, 1.5, 2.5};
    g_u.setBTCLocations(xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    g_u.SolveTransport(t_end_pdf, dt_pdf, "transport_", output_interval_pdf, up_dir, "Cu", &BTCs_Upscaled);

    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    BTCs_Upscaled.write(up_btc_path);
    BTCs_Upscaled.derivative().write(up_btc_deriv_path);

    g_u.writeNamedVTI_Auto("C", joinPath(up_dir, "Cu.vti"));
    g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell, joinPath(up_dir, up_pfx + "Cu.txt"));

    // =====================================================================
    // FINAL AGGREGATION CSVs for comparison / plotting
    // =====================================================================
    if (rank == 0) {

        // ------------------------------------------------------------
        // FIXED UNIFORM TIME GRID (best for consistent plotting)
        // ------------------------------------------------------------
        const double t_end_cmp = 10.0;    // must match SolveTransport(t_end)
        const double dt_cmp    = 0.001;   // plotting resolution you want
        std::vector<double> t_base;
        t_base.reserve((size_t)std::ceil(t_end_cmp / dt_cmp) + 1);
        for (double tt = 0.0; tt <= t_end_cmp + 1e-12; tt += dt_cmp) t_base.push_back(tt);

        // 1) BTC comparison
        std::vector<std::string> out_names;
        std::vector<std::vector<double>> out_cols;

        // Read any BTC CSV and RESAMPLE onto fixed time grid t_base
        auto ingest_one = [&](const std::string& path, const std::string& series_prefix) -> bool {
            std::vector<double> t;
            std::vector<std::string> names;
            std::vector<std::vector<double>> cols;
            if (!read_time_series_table_csv(path, t, names, cols)) {
                std::cerr << "WARNING: could not read CSV for aggregation: " << path << "\n";
                return false;
            }

            std::vector<std::vector<double>> cols_rs;
            resample_table_linear(t, cols, t_base, cols_rs);

            for (size_t j = 0; j < names.size(); ++j) {
                out_names.push_back(series_prefix + "_" + names[j]);
                out_cols.push_back(std::move(cols_rs[j]));
            }
            return true;
        };

        // ingest all fine BTCs
        for (int r = 1; r <= nReal; ++r) {
            std::string fine_dir = joinPath(run_dir, makeFineFolder(r)) + "/";
            std::string pfx = makeRealLabel(r) + "_";
            std::string btc = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
            ingest_one(btc, "Fine_" + makeRealLabel(r));
        }

        // ingest upscaled
        ingest_one(up_btc_path, "Upscaled_mean");

        const std::string out_cmp = joinPath(run_dir, "BTC_Compare_Fine_vs_Upscaled.csv");
        if (!write_comparison_csv(out_cmp, t_base, out_names, out_cols)) {
            std::cerr << "WARNING: failed to write " << out_cmp << "\n";
        } else {
            std::cout << "Wrote comparison BTC CSV: " << out_cmp << "\n";
        }

        // 2) Derivative BTC comparison
        out_names.clear();
        out_cols.clear();

        for (int r = 1; r <= nReal; ++r) {
            std::string fine_dir = joinPath(run_dir, makeFineFolder(r)) + "/";
            std::string pfx = makeRealLabel(r) + "_";
            std::string btc = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");
            ingest_one(btc, "FineDeriv_" + makeRealLabel(r));
        }

        ingest_one(up_btc_deriv_path, "UpscaledDeriv_mean");

        const std::string out_cmp_d = joinPath(run_dir, "BTC_Compare_FineDerivative_vs_UpscaledDerivative.csv");
        if (!write_comparison_csv(out_cmp_d, t_base, out_names, out_cols)) {
            std::cerr << "WARNING: failed to write " << out_cmp_d << "\n";
        } else {
            std::cout << "Wrote comparison BTC-derivative CSV: " << out_cmp_d << "\n";
        }

        std::cout << "\nMixing PDF simulation complete!\n";
        std::cout << "All outputs saved to: " << run_dir << "\n";
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
