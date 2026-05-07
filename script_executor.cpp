// script_executor.cpp
//
// Executes a parsed ScriptProgram by maintaining a ScriptContext (variable map)
// and dispatching each statement to the appropriate sim_runner call.
// -----------------------------------------------------------------------

#include "script_executor.h"
#include "script_parser.h"
#include "script_types.h"

#include "sim_runner.h"
#include "sim_helpers.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <functional>

#include <mpi.h>

// -----------------------------------------------------------------------
// Internal helpers
// -----------------------------------------------------------------------
namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return s;
}

// -----------------------------------------------------------------------
// Build SimParams from context
// -----------------------------------------------------------------------
SimParams buildSimParams(const ScriptContext& ctx)
{
    SimParams P;

    auto getd = [&](const std::string& k, double def) -> double {
        return ctx.has(k) ? std::stod(ctx.expand(ctx.get(k))) : def;
    };
    auto geti = [&](const std::string& k, int def) -> int {
        return ctx.has(k) ? std::stoi(ctx.expand(ctx.get(k))) : def;
    };

    P.nx              = geti("nx",   300);
    P.ny              = geti("ny",   100);
    P.nu              = geti("nu",    20);
    P.Lx              = getd("Lx",   3.0);
    P.Ly              = getd("Ly",   1.0);
    P.nReal_default   = geti("n_realizations", 20);
    P.stdev           = getd("stdev", 1.0);
    P.g_mean          = getd("g_mean", 0.0);
    P.diffusion_factor= getd("diffusion_factor", 1.0);
    P.correlation_ls_x= getd("corr_x", 1.0);
    P.correlation_ls_y= getd("corr_y", 0.1);
    P.Diffusion_coefficient = getd("diffusion", 0.0);
    P.t_end_pdf       = getd("t_end", 20.0);
    P.run_seed        = static_cast<unsigned long>(geti("seed", 42));

    // x_locations: comma-separated list
    if (ctx.has("x_locations")) {
        const std::string raw = ctx.expand(ctx.get("x_locations"));
        std::stringstream ss(raw);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            const auto t = [](std::string s){
                s.erase(0, s.find_first_not_of(" \t"));
                s.erase(s.find_last_not_of(" \t") + 1);
                return s;
            }(tok);
            if (!t.empty()) P.xLocations.push_back(std::stod(t));
        }
    }
    if (P.xLocations.empty()) P.xLocations = {0.5, 1.5, 2.5};

    // Correlation model
    if (ctx.has("correlation_model")) {
        const std::string m = lower(ctx.expand(ctx.get("correlation_model")));
        if      (m == "exponential")  P.CorrelationModel = SimParams::correlationmode::exponentialfit;
        else if (m == "derivative")   P.CorrelationModel = SimParams::correlationmode::derivative;
        else if (m == "gaussian")     P.CorrelationModel = SimParams::correlationmode::gaussian;
        else if (m == "matern")       P.CorrelationModel = SimParams::correlationmode::matern;
        else if (m == "oneoversum")   P.CorrelationModel = SimParams::correlationmode::oneoversums;
    }

    // Correlation sampling ranges  "min,max"
    auto parseRange = [&](const std::string& key,
                          std::pair<double,double>& out,
                          double def_min, double def_max) {
        if (!ctx.has(key)) { out = {def_min, def_max}; return; }
        const std::string raw = ctx.expand(ctx.get(key));
        const auto comma = raw.find(',');
        if (comma == std::string::npos) { out = {std::stod(raw), std::stod(raw)}; return; }
        out = {std::stod(raw.substr(0, comma)), std::stod(raw.substr(comma + 1))};
    };
    parseRange("corr_x_range", P.correlation_x_range, 0.001, 0.02);
    parseRange("corr_y_range", P.correlation_y_range, 0.001, 0.02);

    return P;
}

// -----------------------------------------------------------------------
// Build RunOptions from context
// -----------------------------------------------------------------------
RunOptions buildRunOptions(const ScriptContext& ctx)
{
    RunOptions opts;

    auto getb = [&](const std::string& k, bool def) -> bool {
        if (!ctx.has(k)) return def;
        const std::string v = lower(ctx.expand(ctx.get(k)));
        return (v == "true" || v == "yes" || v == "1");
    };
    auto geti = [&](const std::string& k, int def) -> int {
        return ctx.has(k) ? std::stoi(ctx.expand(ctx.get(k))) : def;
    };

    opts.upscale_only              = getb("upscale_only",   false);
    opts.hardcoded_mean            = getb("hardcoded_mean", false);
    opts.solve_fine_scale_transport= getb("fine_transport", true);
    opts.solve_upscale_transport   = getb("upscale_transport", true);
    opts.perform_particle_tracking = getb("particle_tracking", false);
    opts.perform_upscaled_PT       = getb("upscaled_pt", false);
    opts.analyze_qx_ranks          = getb("qx_ranks", true);
    opts.analyze_qx_copula_diagnostics = getb("copula_analysis", false);
    opts.qx_copula_bootstrap       = geti("copula_bootstrap", 5);
    opts.qx_copula_max_points      = geti("copula_max_points", 1200);
    opts.empirical_copula_bins     = geti("copula_bins", 20);
    opts.qx_rank_num_samples       = geti("qx_rank_num_samples", 10000);
    opts.qx_rank_num_delta         = geti("qx_rank_num_delta", 30);

    // Copula dependence model
    if (ctx.has("copula_model")) {
        const std::string m = lower(ctx.expand(ctx.get("copula_model")));
        if (m == "empirical")
            opts.qx_dependence_model = RunOptions::CopulaDependenceModel::Empirical;
        else
            opts.qx_dependence_model = RunOptions::CopulaDependenceModel::GaussianRank;
    }

    // Mixing model
    if (ctx.has("mixing_model")) {
        const std::string m = lower(ctx.expand(ctx.get("mixing_model")));
        if      (m == "exponential")      opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::Exponential;
        else if (m == "exponentialvdep")  opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::ExponentialVDep;
        else if (m == "gaussian")         opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::Gaussian;
        else if (m == "matern")           opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::Matern;
        else if (m == "oneoversum")       opts.upscaled_mixing_model = RunOptions::UpscaledMixingModel::OneOverSum;
    }

    // Mean TS mode
    if (ctx.has("mean_ts_mode")) {
        const std::string m = lower(ctx.expand(ctx.get("mean_ts_mode")));
        if      (m == "first")   opts.mean_ts_mode = RunOptions::MeanTSMode::First;
        else if (m == "longest") opts.mean_ts_mode = RunOptions::MeanTSMode::Longest;
        else if (m == "union")   opts.mean_ts_mode = RunOptions::MeanTSMode::Union;
    }

    if (ctx.has("qx_cdf_path"))
        opts.hardcoded_qx_cdf_path = ctx.expand(ctx.get("qx_cdf_path"));

    return opts;
}

// -----------------------------------------------------------------------
// Build HardcodedMean from context
// -----------------------------------------------------------------------
HardcodedMean buildHardcodedMean(const ScriptContext& ctx)
{
    HardcodedMean H;
    auto getd = [&](const std::string& k, double def) -> double {
        return ctx.has(k) ? std::stod(ctx.expand(ctx.get(k))) : def;
    };
    H.lc_mean  = getd("hm_lc",   0.0);
    H.lx_mean  = getd("hm_lx",   0.0);
    H.ly_mean  = getd("hm_ly",   0.0);
    H.dt_mean  = getd("hm_dt",   0.0);
    H.nu_x     = getd("hm_nu_x", 1.5);
    H.nu_y     = getd("hm_nu_y", 1.5);
    H.qx_const = getd("hm_qx_const", 0.0);
    return H;
}

// -----------------------------------------------------------------------
// Execute a single "run fine" or "run upscale" statement
// -----------------------------------------------------------------------
bool execRun(const std::string& stage,
             ScriptContext& ctx,
             RunOutputs& out,
             int rank)
{
    const SimParams     P    = buildSimParams(ctx);
    RunOptions          opts = buildRunOptions(ctx);
    const HardcodedMean H    = buildHardcodedMean(ctx);

    const std::string output_dir    = ctx.has("output_dir")    ? ctx.expand(ctx.get("output_dir"))    : "./UpscalingResults";
    const std::string resume_run_dir= ctx.has("resume_run_dir")? ctx.expand(ctx.get("resume_run_dir")): "";

    // Override run mode flags based on stage
    if (stage == "fine") {
        opts.upscale_only = false;
        opts.solve_fine_scale_transport = true;
        // upscale stays as configured
    } else if (stage == "upscale") {
        opts.upscale_only = true;
        opts.solve_fine_scale_transport = false;
    } else if (stage == "copula_transport") {
        // handled by execCopulaTransport — should not reach here
        std::cerr << "[warn] run copula_transport: use dedicated path\n";
        return false;
    } else {
        throw std::runtime_error("execRun: unknown stage '" + stage + "'");
    }

    const std::string run_tag = make_run_tag_std_D_aniso_df(P);
    out.run_dir = prepare_run_dir_mpi(output_dir, resume_run_dir, opts, rank, run_tag);

    // Store run_dir back into context for subsequent write/score statements
    ctx.set("run_dir", out.run_dir);

    const bool ok = run_simulation_blocks(P, opts, H, resume_run_dir, out, rank);
    if (!ok && rank == 0)
        std::cerr << "[error] run_simulation_blocks failed for stage '" << stage << "'\n";
    return ok;
}

// -----------------------------------------------------------------------
// Execute "run copula_transport" using saved copula CSV files
// -----------------------------------------------------------------------
bool execCopulaTransport(ScriptContext& ctx, RunOutputs& out, int rank)
{
    // Resolve copula CSV paths
    auto resolvePath = [&](const std::string& key,
                           const std::string& default_filename) -> std::string
    {
        const std::string val = ctx.has(key) ? ctx.expand(ctx.get(key)) : "auto";
        if (lower(val) != "auto") return val;
        // auto: look in run_dir
        if (!ctx.has("run_dir"))
            throw std::runtime_error("execCopulaTransport: 'auto' path requires run_dir to be set");
        return joinPath(ctx.expand(ctx.get("run_dir")), default_filename);
    };

    const std::string adv_path  = resolvePath("copula_adv_path",
                                               "mean_empirical_copula_advection.csv");
    const std::string diff_path = resolvePath("copula_diff_path",
                                               "mean_empirical_copula_velocity.csv");

    if (rank == 0) {
        std::cout << "run copula_transport:\n"
                  << "  theta_adv  = " << adv_path  << "\n"
                  << "  theta_diff = " << diff_path << "\n";
    }

    // We only need this on rank 0 for now (solver is serial in Grid2D)
    if (rank == 0) {
        CMatrix theta_adv  = Grid2D::loadCopulaCSV(adv_path);
        CMatrix theta_diff = Grid2D::loadCopulaCSV(diff_path);

        const SimParams P = buildSimParams(ctx);

        // lambda_a: use fitted lc_mean if "auto"
        double lambda_a = 1.0;
        if (ctx.has("copula_lambda_a")) {
            const std::string v = lower(ctx.expand(ctx.get("copula_lambda_a")));
            if (v != "auto") lambda_a = std::stod(v);
            // If auto: we'd need lc_mean from a prior fine run stored in context
            // For now we leave it as the default and warn
            else if (ctx.has("lc_mean"))
                lambda_a = std::stod(ctx.get("lc_mean"));
            else
                std::cerr << "[warn] copula_lambda_a = auto but lc_mean not set; using 1.0\n";
        }

        const double dt0 = ctx.has("copula_dt0")
                           ? std::stod(ctx.expand(ctx.get("copula_dt0")))
                           : 0.5;

        // Build v_of_u from the mean inverse CDF stored in run_dir
        // (written as mean_qx_inverse_cdf.txt by run_fine_loop_collect)
        TimeSeries<double> invcdf;
        const std::string cdf_path = ctx.has("invcdf_path")
            ? ctx.expand(ctx.get("invcdf_path"))
            : joinPath(ctx.has("run_dir") ? ctx.expand(ctx.get("run_dir")) : ".",
                       "mean_qx_inverse_cdf.txt");

        if (!read_inverse_cdf_any_format(cdf_path, invcdf)) {
            std::cerr << "[error] run copula_transport: could not load inverse CDF from '"
                      << cdf_path << "'\n";
            return false;
        }

        // v(u) = interpolate invcdf at u
        auto v_of_u = [&invcdf](double u) -> double {
            return invcdf.interpol(u);
        };

        // Build a minimal Grid2D just for the upscaled solver
        Grid2D g(P.nx, P.nu);
        g.setDX(P.Lx / P.nx);
        g.SetVal("diffusion", P.Diffusion_coefficient);
        g.SetVal("c_left",    1.0);
        g.setBTCLocations(P.xLocations);

        // dt: use 0.5 * dx / v_max as a safety bound
        double v_max = v_of_u(1.0 - 1.0 / (2.0 * P.nu));
        double dt    = 0.5 * (P.Lx / P.nx) / std::max(v_max, 1e-12);

        TimeSeriesSet<double> btc_set;
        g.SolveCopulaTransport(
            P.t_end_pdf, dt,
            theta_adv, theta_diff,
            lambda_a, dt0,
            v_of_u,
            "copula_",
            &btc_set,
            100,
            ctx.has("run_dir") ? ctx.expand(ctx.get("run_dir")) : ".");

        // Store the BTC output into out for later write statements
        // Use a single-entry vector to match the xLocation convention
        out.Fine_Scale_BTCs_pdf.resize(P.xLocations.size());
        // Write BTC immediately to auto path
        const std::string btc_out = ctx.has("run_dir")
            ? joinPath(ctx.expand(ctx.get("run_dir")), "copula_BTC.csv")
            : "copula_BTC.csv";
        btc_set.write(btc_out);
        if (rank == 0)
            std::cout << "  Wrote copula BTC: " << btc_out << "\n";
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    return true;
}

// -----------------------------------------------------------------------
// Execute a "write" statement
// -----------------------------------------------------------------------
void execWrite(const std::string& kind,
               const std::string& raw_filename,
               ScriptContext& ctx,
               const RunOutputs& out,
               const SimParams& P,
               int rank)
{
    if (rank != 0) return;   // only rank 0 writes

    const std::string filename = ctx.expand(raw_filename);
    const std::string run_dir  = ctx.has("run_dir")
                                 ? ctx.expand(ctx.get("run_dir")) : ".";
    const std::string path     = (filename.front() == '/' || filename.front() == '.')
                                 ? filename
                                 : joinPath(run_dir, filename);

    std::cout << "write " << kind << " -> " << path << "\n";

    if (kind == "btc_mean") {
        out.mean_BTCs.write(path);
    }
    else if (kind == "btc_mean_pdf") {
        out.mean_transport_pdf.write(path);
    }
    else if (kind == "btc_mean_cdf") {
        out.mean_transport_cdf.write(path);
    }
    else if (kind == "pt_mean") {
        out.mean_pt_pdf.write(path);
    }
    else if (kind == "pt_mean_cdf") {
        out.mean_pt_cdf.write(path);
    }
    else if (kind == "btc_compare") {
        for (int i = 0; i < (int)P.xLocations.size(); ++i) {
            std::ostringstream fn;
            fn << path << "_x" << std::fixed << std::setprecision(2) << P.xLocations[i] << ".csv";
            if (i < (int)out.Fine_Scale_BTCs_pdf.size())
                out.Fine_Scale_BTCs_pdf[i].write(fn.str());
        }
    }
    else if (kind == "params") {
        std::ofstream f(path);
        if (!f) { std::cerr << "[warn] write params: cannot open '" << path << "'\n"; return; }
        f << "run_dir=" << run_dir << "\n";
        for (auto& [k, v] : ctx.vars())
            f << k << "=" << v << "\n";
    }
    else {
        std::cerr << "[warn] write: unknown kind '" << kind << "' (ignored)\n";
    }
}

// -----------------------------------------------------------------------
// Execute a "score" statement
// -----------------------------------------------------------------------
void execScore(const std::string& kind,
               const std::string& ref_path_raw,
               ScriptContext& ctx,
               const RunOutputs& out,
               int rank)
{
    if (rank != 0) return;

    const std::string ref_path = ctx.expand(ref_path_raw);
    const std::string run_dir  = ctx.has("run_dir") ? ctx.expand(ctx.get("run_dir")) : ".";

    std::cout << "score " << kind << " against " << ref_path << "\n";

    const double s = score_upscaled_vs_black_mean_from_compare(ref_path, run_dir);
    std::cout << "  score = " << s << "\n";

    // Store score in context so it can be referenced by $score
    ctx.set("score", std::to_string(s));
    ctx.set("last_score", std::to_string(s));
}

// -----------------------------------------------------------------------
// Execute a "plot" statement
// -----------------------------------------------------------------------
void execPlot(const std::string& kind,
              const std::string& raw_filename,
              ScriptContext& ctx,
              int rank)
{
    if (rank != 0) return;

    const std::string filename = ctx.expand(raw_filename);
    const std::string run_dir  = ctx.has("run_dir") ? ctx.expand(ctx.get("run_dir")) : ".";
    const std::string path     = (filename.front() == '/' || filename.front() == '.')
                                 ? filename
                                 : joinPath(run_dir, filename);

    std::cout << "plot " << kind << " -> " << path << "\n";

    if (kind == "btc") {
        // Write a gnuplot script and execute it
        const std::string gp_path = path + ".gp";
        const bool ok = write_btc_compare_plot_gnuplot_by_basename(
            gp_path, path, path, "c/c0", /*skip_base_t=*/true);
        if (ok) run_gnuplot_script(gp_path);
    } else {
        std::cerr << "[warn] plot: unknown kind '" << kind << "' (ignored)\n";
    }
}

// -----------------------------------------------------------------------
// Forward declaration for recursive sweep execution
// -----------------------------------------------------------------------
bool execStatements(const std::vector<Statement>& stmts,
                    ScriptContext& ctx,
                    RunOutputs& out,
                    int rank);

// -----------------------------------------------------------------------
// Execute a sweep block
// -----------------------------------------------------------------------
bool execSweep(const Statement& stmt,
               ScriptContext& ctx,
               RunOutputs& out,
               int rank)
{
    for (const std::string& raw_val : stmt.sweep_values) {
        // Expand value (it might reference other variables)
        const std::string val = ctx.expand(raw_val);

        // Set the sweep variable
        ctx.set(stmt.sweep_var, val);

        if (rank == 0)
            std::cout << "\n=== sweep " << stmt.sweep_var
                      << " = " << val << " ===\n";

        // Execute body
        const bool ok = execStatements(stmt.body, ctx, out, rank);
        if (!ok) {
            std::cerr << "[error] sweep body failed at "
                      << stmt.sweep_var << " = " << val << "\n";
            return false;
        }

        // After each sweep iteration, write mean outputs if run_dir is set
        // (averaging of copulas, correlations, BTCs happened inside run_simulation_blocks)
        if (rank == 0 && ctx.has("run_dir")) {
            // Automatic mean writes that always make sense after a fine run
            const std::string run_dir = ctx.expand(ctx.get("run_dir"));
            if (!out.mean_BTCs.empty())
                out.mean_BTCs.write(joinPath(run_dir, "BTC_mean.csv"));
            if (!out.mean_transport_pdf.empty())
                out.mean_transport_pdf.write(joinPath(run_dir, "BTC_mean_pdf.csv"));
            if (!out.mean_transport_cdf.empty())
                out.mean_transport_cdf.write(joinPath(run_dir, "BTC_mean_cdf.csv"));
            if (!out.mean_pt_pdf.empty())
                out.mean_pt_pdf.write(joinPath(run_dir, "PT_mean_pdf.csv"));
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
    return true;
}

// -----------------------------------------------------------------------
// Execute a list of statements (top-level or sweep body)
// -----------------------------------------------------------------------
bool execStatements(const std::vector<Statement>& stmts,
                    ScriptContext& ctx,
                    RunOutputs& out,
                    int rank)
{
    for (const auto& stmt : stmts) {
        switch (stmt.kind) {

        case StmtKind::Blank:
        case StmtKind::Comment:
            break;

        case StmtKind::Set: {
            // Expand RHS at assignment time so $vars are resolved now
            const std::string val = ctx.expand(stmt.value);
            ctx.set(stmt.key, val);
            if (rank == 0)
                std::cout << "set " << stmt.key << " = " << val << "\n";
            break;
        }

        case StmtKind::Run: {
            const std::string stage = lower(stmt.key);
            if (rank == 0)
                std::cout << "\n>>> run " << stage << "\n";

            bool ok = false;
            if (stage == "copula_transport") {
                ok = execCopulaTransport(ctx, out, rank);
            } else {
                ok = execRun(stage, ctx, out, rank);
                // After fine run, store fitted means in context
                // so "auto" lookups and $lc_mean work in subsequent statements
                if (ok && stage == "fine" && rank == 0) {
                    // These were computed inside run_simulation_blocks via
                    // FineScaleOutputs and written to mean_params.txt
                    // Re-read them from file into context for downstream use
                    const std::string params_path = ctx.has("run_dir")
                        ? joinPath(ctx.expand(ctx.get("run_dir")), "mean_params.txt")
                        : "mean_params.txt";
                    double lc=0, lx=0, ly=0, dt=0, nu_x=0, nu_y=0;
                    if (read_mean_params_txt(params_path, lc, lx, ly, dt, nu_x, nu_y)) {
                        ctx.set("lc_mean",  std::to_string(lc));
                        ctx.set("lx_mean",  std::to_string(lx));
                        ctx.set("ly_mean",  std::to_string(ly));
                        ctx.set("dt_mean",  std::to_string(dt));
                        ctx.set("nu_x_mean",std::to_string(nu_x));
                        ctx.set("nu_y_mean",std::to_string(nu_y));
                    }
                }
            }
            if (!ok) {
                std::cerr << "[error] run " << stage << " failed at line "
                          << stmt.lineno << "\n";
                return false;
            }
            break;
        }

        case StmtKind::Write: {
            const SimParams P = buildSimParams(ctx);
            execWrite(stmt.key, stmt.value, ctx, out, P, rank);
            MPI_Barrier(PETSC_COMM_WORLD);
            break;
        }

        case StmtKind::Score:
            execScore(stmt.value, stmt.key, ctx, out, rank);
            MPI_Barrier(PETSC_COMM_WORLD);
            break;

        case StmtKind::Plot:
            execPlot(stmt.key, stmt.value, ctx, rank);
            break;

        case StmtKind::SweepBegin: {
            const bool ok = execSweep(stmt, ctx, out, rank);
            if (!ok) return false;
            break;
        }

        case StmtKind::SweepEnd:
            // Should never appear in a flat list (parser nests these)
            throw std::runtime_error("execStatements: unexpected SweepEnd at line "
                                     + std::to_string(stmt.lineno));
        }
    }
    return true;
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------
bool executeScript(const ScriptProgram& program, int rank)
{
    ScriptContext ctx;
    RunOutputs    out;

    // Broadcast: all ranks execute the same statement list,
    // but only rank 0 does I/O (write/plot/score).
    // run* calls invoke MPI-aware run_simulation_blocks internally.

    try {
        return execStatements(program, ctx, out, rank);
    } catch (const std::exception& e) {
        if (rank == 0)
            std::cerr << "[fatal] Script execution error: " << e.what() << "\n";
        return false;
    }
}

bool runScriptFile(const std::string& path, int rank)
{
    try {
        if (rank == 0)
            std::cout << "=== Loading script: " << path << " ===\n";

        ScriptProgram program = parseScriptFile(path);

        if (rank == 0)
            std::cout << "=== Parsed " << program.size()
                      << " top-level statements ===\n\n";

        return executeScript(program, rank);

    } catch (const std::exception& e) {
        if (rank == 0)
            std::cerr << "[fatal] " << e.what() << "\n";
        return false;
    }
}
