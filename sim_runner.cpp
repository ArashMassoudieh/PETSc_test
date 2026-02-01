// sim_runner.cpp

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

#include "grid.h"
#include "Pathway.h"
#include "PathwaySet.h"

static inline double du_from_nu(int nu)
{
    return (nu > 0) ? (1.0 / double(nu)) : 1.0;
}

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

        if (opts.upscale_only) {
            if (opts.hardcoded_mean) {
                // ALWAYS create a fresh run directory
                std::string name = "run_" + makeTimestamp();
                if (!run_tag.empty()) name += "_" + run_tag;
                run_dir = joinPath(output_dir, name);
                createDirectory(run_dir);

                std::cout << "Hardcoded-mean upscaled-only: created NEW run_dir: "
                          << run_dir << "\n";
            }
            else {
                // strict upscale-only: must point to existing run dir
                if (resume_run_dir.empty()) {
                    std::cerr << "ERROR: --upscale-only requires --run-dir=/path/to/run_YYYYMMDD_HHMMSS\n";
                    MPI_Abort(PETSC_COMM_WORLD, 1);
                }
                run_dir = resume_run_dir;
                if (!dirExists(run_dir)) {
                    std::cerr << "ERROR: run_dir does not exist: " << run_dir << "\n";
                    MPI_Abort(PETSC_COMM_WORLD, 2);
                }
                std::cout << "Upscale-only: using run_dir: " << run_dir << "\n";
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

// Load an inverse CDF from a file, tolerant to header/no-header.
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

// Build mean inverse-CDF from the multi-series output (qx_inverse_cdfs.txt)
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

        // Expect pairs: (t,q), (t,q), ...
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

static bool run_fine_loop_collect(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    std::vector<double>& lc_all,
    std::vector<double>& lx_all,
    std::vector<double>& ly_all,
    std::vector<double>& dt_all,
    TimeSeriesSet<double>& inverse_qx_cdfs,
    TimeSeriesSet<double>& qx_pdfs,
    TimeSeriesSet<double>& lambda_x_correlations,
    TimeSeriesSet<double>& lambda_y_correlations,
    TimeSeriesSet<double>& lambda_a_correlations,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;
    const double du = du_from_nu(nu);

    const std::string stats_csv = joinPath(run_dir, "fine_params_all.csv");
    if (rank == 0) {
        std::ofstream f(stats_csv);
        f << "realization,lc,lambda_x,lambda_y,dt_opt\n";
    }

    const int nReal = P.nReal_default;

    for (int r = 1; r <= nReal; ++r) {
        const std::string rlab = makeRealLabel(r);
        std::string fine_dir = joinPath(run_dir, makeFineFolder(r));
        if (rank == 0) createDirectory(fine_dir);
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        const std::string pfx = rlab + "_";

        Grid2D g(nx, ny, Lx, Ly);

        PetscLogDouble t_total0, t_total1, t_asm0, t_asm1, t_solve0, t_solve1;

        unsigned long seed = P.run_seed + 1000UL*(unsigned long)r + (unsigned long)rank;

        g.makeGaussianFieldSGS("K_normal_score", P.correlation_ls_x, P.correlation_ls_y, 10, seed);
        g.normalizeField("K_normal_score", Grid2D::ArrayKind::Cell);

        PetscTime(&t_asm0);
        PetscTime(&t_total0);

        g.writeNamedMatrix("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.txt"));
        g.writeNamedVTI("K_normal_score", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K_normal_score.vti"));

        g.createExponentialField("K_normal_score", P.stdev, P.g_mean, "K");
        g.writeNamedVTI("K", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.vti"));

        PetscTime(&t_asm1);

        PetscTime(&t_solve0);
        g.DarcySolve(Lx, 0, "K", "K");
        std::cout << "Darcy solved ... " << std::endl;
        PetscTime(&t_solve1);

        g.computeMassBalanceError("MassBalanceError");
        g.writeNamedMatrix("MassBalanceError", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "error.txt"));

        g.writeNamedVTI("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.vti"));
        g.writeNamedMatrix("head", Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "Head.txt"));
        g.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.vti"));
        g.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.vti"));


        g.computeMassBalanceError("div_q");

        // Find cells with non-zero divergence
        const auto& div_q = g.field("div_q");
        double max_div = 0.0;
        int worst_i = -1, worst_j = -1;

        for (int j = 0; j < g.ny(); ++j) {
            for (int i = 0; i < g.nx(); ++i) {
                double div = div_q[g.cellIndex(i,j)];
                if (std::abs(div) > std::abs(max_div)) {
                    max_div = div;
                    worst_i = i;
                    worst_j = j;
                }
            }
        }

        std::cout << "Max divergence: " << max_div
                  << " at cell (" << worst_i << "," << worst_j << ")" << std::endl;


        TimeSeries<double> AllQxValues = g.exportFieldToTimeSeries("qx", Grid2D::ArrayKind::Fx);
        TimeSeries<double> QxNormalScores = AllQxValues.ConvertToNormalScore();
        g.assignFromTimeSeries(QxNormalScores, "qx_normal_score", Grid2D::ArrayKind::Fx);
        g.writeNamedVTI("qx_normal_score", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx_normal_score.vti"));

        // velocity autocorrelation (X,Y) via perturbation
        int num_deltas = 30;
        int num_samples_per_delta = 10000;

        TimeSeries<double> corr_x, corr_y;

        for (int i = 0; i < num_deltas; ++i) {
            double exponent = static_cast<double>(i) / (num_deltas - 1);
            double delta = P.correlation_x_range.first * std::pow(P.correlation_x_range.second / P.correlation_x_range.first, exponent);
            try {
                TimeSeries<double> samples = g.sampleGaussianPerturbation(
                    "qx_normal_score", Grid2D::ArrayKind::Fx,
                    num_samples_per_delta, delta, 0, PerturbDir::XOnly);
                corr_x.append(delta, samples.correlation_tc());
            } catch (...) {}
        }
        corr_x.writefile(joinPath(fine_dir, pfx + "velocity_correlation_x.txt"));
        lambda_x_correlations.append(corr_x, "Realization" + aquiutils::numbertostring(r + 1));
        double lambda_x_emp;
        if (P.CorrelationModel == SimParams::correlationmode::exponentialfit)
            lambda_x_emp= P.lambda_y_multiplier*corr_x.fitExponentialDecay();
        else if (P.CorrelationModel == SimParams::correlationmode::derivative)
            lambda_x_emp= corr_x.getTime(0)/(1.0-corr_x.getValue(0));
        else
            lambda_x_emp= corr_x.fitGaussianDecay();


        for (int i = 0; i < num_deltas; ++i) {
            double exponent = static_cast<double>(i) / (num_deltas - 1);
            double delta = P.correlation_y_range.first * std::pow(P.correlation_y_range.second / P.correlation_y_range.first, exponent);
            try {
                TimeSeries<double> samples = g.sampleGaussianPerturbation(
                    "qx_normal_score", Grid2D::ArrayKind::Fx,
                    num_samples_per_delta, delta, 0, PerturbDir::YOnly);
                corr_y.append(delta, samples.correlation_tc());
            } catch (...) {}
        }
        corr_y.writefile(joinPath(fine_dir, pfx + "velocity_correlation_y.txt"));
        lambda_y_correlations.append(corr_y, "Realization" + aquiutils::numbertostring(r + 1));

        double lambda_y_emp;
        if (P.CorrelationModel == SimParams::correlationmode::exponentialfit)
            lambda_y_emp= P.lambda_y_multiplier*corr_y.fitExponentialDecay();
        else if (P.CorrelationModel == SimParams::correlationmode::derivative)
            lambda_y_emp= corr_y.getTime(0)/(1.0-corr_y.getValue(0));
        else
            lambda_y_emp= corr_y.fitGaussianDecay();

        // inverse CDF + pdf
        TimeSeries<double> qx_inverse_cdf = g.extractFieldCDF("qx", Grid2D::ArrayKind::Fx, 100, 1e-6);
        TimeSeries<double> qx_pdf = g.extractFieldPDF("qx", Grid2D::ArrayKind::Fx, 50, 1e-6, true);
        qx_inverse_cdf = qx_inverse_cdf.make_uniform(du);
        qx_inverse_cdf.writefile(joinPath(fine_dir, pfx + "qx_inverse_cdf.txt"));
        qx_pdf.writefile(joinPath(fine_dir, pfx + "qx_pdf.txt"));

        // dt
        double dt_optimal = 0.5 * g.dx() / g.fieldMinMax("qx", Grid2D::ArrayKind::Fx).second;

        // transport
        g.assignConstant("C", Grid2D::ArrayKind::Cell, 0);

        g.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(fine_dir, pfx + "qx.txt"));
        g.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "qy.txt"));
        g.writeNamedMatrix("K",  Grid2D::ArrayKind::Cell, joinPath(fine_dir, pfx + "K.txt"));

        g.SetVal("diffusion", P.Diffusion_coefficient);
        g.assignConstant("D_y", Grid2D::ArrayKind::Fy, P.Diffusion_coefficient);
        g.writeNamedVTI("D_y",Grid2D::ArrayKind::Fy, joinPath(fine_dir, pfx + "D_y.vti"));
        g.SetVal("porosity", 1.0);
        g.SetVal("c_left", 1.0);

        TimeSeriesSet<double> BTCs_FineScaled;
        g.setBTCLocations(P.xLocations);

        if (opts.solve_fine_scale_transport) {
            g.SolveTransport(
                P.t_end_pdf,
                std::min(dt_optimal, 0.5 / 10.0),
                "transport_", 500,
                fine_dir,
                "C",
                &BTCs_FineScaled,
                r
            );

            for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_FineScaled.size(); ++i) {
                BTCs_FineScaled.setSeriesName(i, fmt_x(P.xLocations[i]));
                Fine_Scale_BTCs[i].append(BTCs_FineScaled[i].derivative(), pfx + fmt_x(P.xLocations[i]));
            }

            const std::string btc_path       = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
            const std::string btc_deriv_path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");
            BTCs_FineScaled.write(btc_path);
            BTCs_FineScaled.derivative().write(btc_deriv_path);
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
        lambda_a_correlations.append(qx_correlation, "Realization" + aquiutils::numbertostring(r + 1));
        double lc_emp = qx_correlation.fitExponentialDecay();

        pathways.writeToFile(joinPath(fine_dir, pfx + "pathway_summary.txt"));
        pathways.writeCombinedVTK(joinPath(fine_dir, pfx + "all_pathways.vtk"));

        // meta
        if (rank == 0) {
            std::ofstream m(joinPath(fine_dir, pfx + "meta.txt"));
            m << "realization=" << r << "\n";
            m << "nx=" << nx << "\nny=" << ny << "\nnu=" << nu << "\n";
            m << "Lx=" << Lx << "\nLy=" << Ly << "\n";
            m << "D=" << P.Diffusion_coefficient << "\n";
            m << "correlation_ls_x=" << P.correlation_ls_x << "\n";
            m << "correlation_ls_y=" << P.correlation_ls_y << "\n";
            m << "stdev=" << P.stdev << "\n";
            m << "g_mean=" << P.g_mean << "\n";
            m << "lc=" << lc_emp << "\n";
            m << "lambda_x=" << lambda_x_emp << "\n";
            m << "lambda_y=" << lambda_y_emp << "\n";
            m << "dt_opt=" << dt_optimal << "\n";
            m << "seed=" << seed << "\n";
        }

        // accumulate means
        if (rank == 0) {
            lc_all.push_back(lc_emp);
            lx_all.push_back(lambda_x_emp);
            ly_all.push_back(lambda_y_emp);
            dt_all.push_back(dt_optimal);

            inverse_qx_cdfs.append(qx_inverse_cdf, "qx_inverse_cdf" + aquiutils::numbertostring(r + 1));
            qx_pdfs.append(qx_pdf, "qx_pdf" + aquiutils::numbertostring(r + 1));

            std::ofstream f(stats_csv, std::ios::app);
            f << r << "," << lc_emp << "," << lambda_x_emp << "," << lambda_y_emp << "," << dt_optimal << "\n";
        }

        if (rank == 0) {
            std::cout << "[Fine " << rlab << "] Assembly time: " << (t_asm1 - t_asm0) << " s\n";
            std::cout << "[Fine " << rlab << "] Solve time:    " << (t_solve1 - t_solve0) << " s\n";
            std::cout << "[Fine " << rlab << "] Total time:   " << (t_total1 - t_total0) << " s\n";
            std::cout << "[Fine " << rlab << "] Outputs saved to: " << fine_dir << "\n";
        }

        MPI_Barrier(PETSC_COMM_WORLD);
    }

    // write run-level outputs (only when we ran fine)
    if (rank == 0) {
        inverse_qx_cdfs.write(joinPath(run_dir, "qx_inverse_cdfs.txt"));
        qx_pdfs.write(joinPath(run_dir, "qx_pdfs.txt"));
        qx_pdfs.mean_ts().writefile(joinPath(run_dir, "qx_mean_pdf.txt"));

        lambda_a_correlations.write(joinPath(run_dir, "advective_correlations.txt"));
        lambda_x_correlations.write(joinPath(run_dir, "diffusion_x_correlations.txt"));
        lambda_y_correlations.write(joinPath(run_dir, "diffusion_y_correlations.txt"));
    }

    return true;
}

static bool build_mean_for_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    const std::string& run_dir,
    const std::vector<double>& lc_all,
    const std::vector<double>& lx_all,
    const std::vector<double>& ly_all,
    const std::vector<double>& dt_all,
    const TimeSeriesSet<double>& inverse_qx_cdfs,
    double& lc_mean,
    double& lx_mean,
    double& ly_mean,
    double& dt_mean,
    TimeSeries<double>& invcdf_mean,
    int rank)
{
    const double du = du_from_nu(P.nu);
    (void)du;

    if (rank != 0) return true; // only rank0 builds; broadcast later

    const std::string run_cdf_copy_path = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

    if (opts.hardcoded_mean) {
        lc_mean = H.lc_mean;
        lx_mean = H.lx_mean;
        ly_mean = H.ly_mean;
        dt_mean = H.dt_mean;

        // --------------------------------------------
        // qx inverse-CDF ladder:
        // 1) load mean_qx_inverse_cdf.txt (opts.hardcoded_qx_cdf_path)
        // 2) if missing, compute from qx_inverse_cdfs.txt in the same folder
        // 3) else flat constant
        // --------------------------------------------
        bool loaded = false;

        // (1) Try load preferred mean file
        loaded = try_load_hardcoded_invcdf(invcdf_mean, opts.hardcoded_qx_cdf_path);

        // (2) If mean missing, compute from qx_inverse_cdfs.txt (if exists)
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
                // If mean file exists now (maybe created earlier), try load it anyway
                if (!loaded && fileExists(meanf)) {
                    loaded = try_load_hardcoded_invcdf(invcdf_mean, meanf);
                }
            }
        }

        // (3) Flat fallback
        if (!loaded) {
            build_flat_invcdf(invcdf_mean, H.qx_const);
        }

        // Always write a copy into this NEW run_dir
        invcdf_mean.writefile(run_cdf_copy_path);

        std::ofstream f(joinPath(run_dir, "mean_params_used.txt"));
        f << "lc_mean=" << lc_mean << "\n";
        f << "lambda_x_mean=" << lx_mean << "\n";
        f << "lambda_y_mean=" << ly_mean << "\n";
        f << "dt_mean=" << dt_mean << "\n";
        f << "mean_qx_inverse_cdf_written_to=" << run_cdf_copy_path << "\n";
        f << "qx_cdf_preferred_path=" << opts.hardcoded_qx_cdf_path << "\n";
        if (!loaded) f << "qx_const=" << H.qx_const << "\n";

        if (loaded) {
            std::cout << "Using HARD-CODED mean params + qx inverse CDF (loaded/derived)\n";
            std::cout << "Copied to: " << run_cdf_copy_path << "\n";
        } else {
            std::cout << "Using HARD-CODED mean params + CONSTANT qx (no mean/multi files).\n";
            std::cout << "Wrote flat CDF to: " << run_cdf_copy_path << "\n";
        }
        return true;
    }

    if (!opts.upscale_only) {
        // mean from current fine run
        lc_mean = mean_of(lc_all);
        lx_mean = mean_of(lx_all);
        ly_mean = mean_of(ly_all);
        dt_mean = mean_of(dt_all);

        invcdf_mean = inverse_qx_cdfs.mean_ts();
        invcdf_mean.writefile(run_cdf_copy_path);

        std::ofstream f(joinPath(run_dir, "mean_params.txt"));
        f << "nReal=" << (int)lc_all.size() << "\n";
        f << "lc_mean=" << lc_mean << "\n";
        f << "lambda_x_mean=" << lx_mean << "\n";
        f << "lambda_y_mean=" << ly_mean << "\n";
        f << "dt_mean=" << dt_mean << "\n";
        f << "du=" << du_from_nu(P.nu) << "\n";

        return true;
    }

    // strict upscale-only: MUST load mean files; NO folder scanning here.
    const std::string mean_params_path = joinPath(run_dir, "mean_params.txt");
    const std::string mean_cdf_path    = joinPath(run_dir, "mean_qx_inverse_cdf.txt");

    bool ok_params = fileExists(mean_params_path) &&
                     read_mean_params_txt(mean_params_path, lc_mean, lx_mean, ly_mean, dt_mean);

    TimeSeries<double> mean_qx_cdf;
    bool ok_cdf = fileExists(mean_cdf_path) &&
                  read_inverse_cdf_any_format(mean_cdf_path, mean_qx_cdf);

    if (!ok_params || !ok_cdf) {
        std::cerr << "ERROR: --upscale-only requires existing mean_params.txt and mean_qx_inverse_cdf.txt in:\n"
                  << "  " << run_dir << "\n"
                  << "If you don't have them, run with --hardcoded-mean.\n";
        MPI_Abort(PETSC_COMM_WORLD, 99);
    }

    invcdf_mean = mean_qx_cdf;
    std::cout << "Loaded mean_params.txt and mean_qx_inverse_cdf.txt\n";
    return true;
}

static bool run_upscaled(
    const SimParams& P,
    const RunOptions& opts,
    const std::string& run_dir,
    double lc_mean,
    double lx_mean,
    double ly_mean,
    double dt_mean,
    const TimeSeries<double>& invcdf_mean,
    std::string& up_dir_out,
    std::string& up_btc_path_out,
    std::string& up_btc_deriv_path_out,
    std::vector<TimeSeriesSet<double>>& Fine_Scale_BTCs,
    int rank)
{
    const int nx = P.nx, ny = P.ny, nu = P.nu;
    const double Lx = P.Lx, Ly = P.Ly;

    std::string up_dir = joinPath(run_dir, "upscaled_mean");
    if (rank == 0) createDirectory(up_dir);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (!up_dir.empty() && up_dir.back() != '/' && up_dir.back() != '\\') up_dir += "/";

    const std::string up_pfx = std::string("upscaled_");

    if (rank == 0) {
        std::cout << "\n=== UPSCALED RUN USING MEAN PARAMETERS ===\n";
        std::cout << "lc_mean=" << lc_mean << "\n";
        std::cout << "lambda_x_mean=" << lx_mean << "\n";
        std::cout << "lambda_y_mean=" << ly_mean << "\n";
        std::cout << "dt_mean=" << dt_mean << "\n";
        std::cout << "upscaled output: " << up_dir << "\n";
    }

    Grid2D g_u(nx, ny, Lx, Ly);

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
            if (id >= qx_size) {
                std::cerr << "ERROR: qx index out of bounds: " << id << " >= " << qx_size << "\n";
                return false;
            }
            qx_u[id] = v_at_u;
        }
    }

    g_u.writeNamedVTI("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.vti"));
    g_u.writeNamedVTI("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.vti"));
    g_u.writeNamedMatrix("qx", Grid2D::ArrayKind::Fx, joinPath(up_dir, up_pfx + "qx_u_initial.txt"));
    g_u.writeNamedMatrix("qy", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "qy_u_initial.txt"));

    g_u.setMixingParams(lc_mean, lx_mean, ly_mean);

    const double dt_pdf = dt_mean;
    const int output_interval_pdf = 500;

    g_u.SetVal("diffusion", P.Diffusion_coefficient);
    g_u.SetVal("porosity", 1.0);
    g_u.SetVal("c_left", 1.0);

    g_u.computeMixingDiffusionCoefficient();
    g_u.writeNamedVTI("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.vti"));
    g_u.writeNamedMatrix("D_y", Grid2D::ArrayKind::Fy, joinPath(up_dir, up_pfx + "D_y.txt"));

    g_u.assignConstant("C", Grid2D::ArrayKind::Cell, 0);
    g_u.setBTCLocations(P.xLocations);

    TimeSeriesSet<double> BTCs_Upscaled;
    const std::string up_btc_path       = joinPath(up_dir, up_pfx + "BTC_Upscaled.csv");
    const std::string up_btc_deriv_path = joinPath(up_dir, up_pfx + "BTC_Upscaled_derivative.csv");

    if (opts.solve_upscale_transport) {
        g_u.SolveTransport(P.t_end_pdf, dt_pdf, "transport_", output_interval_pdf, up_dir, "Cu", &BTCs_Upscaled);

        for (int i = 0; i < (int)P.xLocations.size() && i < (int)BTCs_Upscaled.size(); ++i) {
            BTCs_Upscaled.setSeriesName(i, fmt_x(P.xLocations[i]));
            Fine_Scale_BTCs[i].append(BTCs_Upscaled[i].derivative(), "Upscaled" + fmt_x(P.xLocations[i]));
        }

        BTCs_Upscaled.write(up_btc_path);
        BTCs_Upscaled.derivative().write(up_btc_deriv_path);

        g_u.writeNamedVTI_Auto("C", joinPath(up_dir, "Cu.vti"));
        g_u.writeNamedMatrix("C", Grid2D::ArrayKind::Cell, joinPath(up_dir, up_pfx + "Cu.txt"));
    }

    up_dir_out = up_dir;
    up_btc_path_out = up_btc_path;
    up_btc_deriv_path_out = up_btc_deriv_path;
    return true;
}

bool run_simulation_blocks(
    const SimParams& P,
    const RunOptions& opts,
    const HardcodedMean& H,
    RunOutputs& out,
    int rank)
{
    out.Fine_Scale_BTCs.clear();
    out.Fine_Scale_BTCs.resize(P.xLocations.size());

    std::vector<double> lc_all, lx_all, ly_all, dt_all;
    TimeSeriesSet<double> inverse_qx_cdfs;
    TimeSeriesSet<double> qx_pdfs;
    TimeSeriesSet<double> lambda_x_correlations;
    TimeSeriesSet<double> lambda_y_correlations;
    TimeSeriesSet<double> lambda_a_correlations;

    if (!opts.upscale_only) {
        const bool ok = run_fine_loop_collect(
            P, opts, out.run_dir,
            lc_all, lx_all, ly_all, dt_all,
            inverse_qx_cdfs, qx_pdfs,
            lambda_x_correlations, lambda_y_correlations, lambda_a_correlations,
            out.Fine_Scale_BTCs,
            rank
        );
        if (!ok) return false;
    }

    out.mean_BTCs.clear();
    for (int i = 0; i < (int)P.xLocations.size(); ++i) {
        out.mean_BTCs.append(out.Fine_Scale_BTCs[i].mean_ts(), fmt_x(P.xLocations[i]));
    }

    double lc_mean=0, lx_mean=0, ly_mean=0, dt_mean=0;
    TimeSeries<double> invcdf_mean;

    build_mean_for_upscaled(
        P, opts, H, out.run_dir,
        lc_all, lx_all, ly_all, dt_all, inverse_qx_cdfs,
        lc_mean, lx_mean, ly_mean, dt_mean, invcdf_mean,
        rank
    );

    MPI_Bcast(&lc_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&lx_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&ly_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dt_mean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    int nU_b = 0;
    if (rank == 0) nU_b = (int)invcdf_mean.size();
    MPI_Bcast(&nU_b, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (rank != 0) invcdf_mean.resize(nU_b);
    MPI_Bcast(invcdf_mean.data(), nU_b, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    std::string up_dir, up_btc_path, up_btc_deriv_path;
    const bool ok_up = run_upscaled(
        P, opts, out.run_dir,
        lc_mean, lx_mean, ly_mean, dt_mean,
        invcdf_mean,
        up_dir, up_btc_path, up_btc_deriv_path,
        out.Fine_Scale_BTCs,
        rank
    );
    if (!ok_up) return false;

    out.up_dir = up_dir;
    out.up_btc_path = up_btc_path;
    out.up_btc_deriv_path = up_btc_deriv_path;

    return true;
}

