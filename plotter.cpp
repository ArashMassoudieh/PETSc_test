// plotter.cpp
#include "plotter.h"

#include "sim_helpers.h"   // read_time_series_table_csv, resample_col_linear, list_fine_folders, etc.
#include "TimeSeries.h"
#include "TimeSeriesSet.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <cstdlib>
#include <algorithm>

// -------------------------------------------------------------
// local helper: load CSV time column only
// -------------------------------------------------------------
static inline bool load_csv_time_only_local(const std::string& path, std::vector<double>& t_out)
{
    std::vector<double> t;
    std::vector<std::string> names;
    std::vector<std::vector<double>> cols;
    if (!read_time_series_table_csv(path, t, names, cols)) return false;
    if (t.empty()) return false;
    t_out = std::move(t);
    return true;
}

bool run_final_aggregation_and_plots(
    const std::string& run_dir_in,
    const std::string& up_btc_path,
    const std::string& up_btc_deriv_path,
    TBaseMode tbase_mode,
    AlignMode align_mode,
    bool use_timeseriesset_mean,
    double t_end_cmp,
    double dt_cmp
)
{
    std::string run_dir = run_dir_in;
    if (!run_dir.empty() && run_dir.back() != '/' && run_dir.back() != '\\') run_dir += "/";

    // Need fine folders for --t-fine and for building FineMean stacks
    auto fine_folders = list_fine_folders(run_dir);
    // NOTE: In upscale-only runs you might still have fine_* folders; if not, some parts will warn & skip.

    // -----------------------
    // choose t_base
    // -----------------------
    std::vector<double> t_base;

    auto build_fixed = [&](){
        t_base.clear();
        t_base.reserve((size_t)std::ceil(t_end_cmp / dt_cmp) + 1);
        for (double tt = 0.0; tt <= t_end_cmp + 1e-12; tt += dt_cmp) t_base.push_back(tt);
    };

    if (tbase_mode == TBaseMode::Fixed) {
        build_fixed();
    } else if (tbase_mode == TBaseMode::FromUpscaled) {
        if (!load_csv_time_only_local(up_btc_path, t_base)) {
            std::cerr << "WARNING: --t-upscaled failed; using fixed t_base\n";
            build_fixed();
        }
    } else { // FromFirstFine
        bool ok = false;
        for (auto& pr : fine_folders) {
            int rr = pr.first;
            std::string fine_dir = pr.second;
            if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";
            std::string pfx  = makeRealLabel(rr) + "_";
            std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");
            ok = load_csv_time_only_local(path, t_base);
            if (ok) break;
        }
        if (!ok) {
            std::cerr << "WARNING: --t-fine failed; using fixed t_base\n";
            build_fixed();
        }
    }

    // -----------------------
    // helper: load CSV -> append to TimeSeriesSet on chosen t_base
    // alignment behavior:
    //   Resample: direct interpolation onto t_base
    //   MakeUniform: first uniformize on dt_cmp, then interpolate onto t_base
    // -----------------------
    auto load_csv_as_tsset = [&](const std::string& path,
                                 const std::string& prefix,
                                 TimeSeriesSet<double>& out) -> bool
    {
        std::vector<double> t;
        std::vector<std::string> names;
        std::vector<std::vector<double>> cols;

        if (!read_time_series_table_csv(path, t, names, cols)) {
            std::cerr << "WARNING: could not read CSV: " << path << "\n";
            return false;
        }

        bool any = false;

        for (size_t j = 0; j < names.size(); ++j) {
            if (names[j] == "t" || names[j] == "T" || names[j] == "time") continue;

            std::vector<double> y_aligned;

            if (align_mode == AlignMode::Resample) {
                y_aligned = resample_col_linear(t, cols[j], t_base);
            } else {
                TimeSeries<double> ts0;
                ts0.reserve(t.size());
                for (size_t i = 0; i < t.size(); ++i) ts0.append(t[i], cols[j][i]);

                TimeSeries<double> tsu = ts0.make_uniform(dt_cmp, false);

                std::vector<double> tu, yu;
                tu.reserve(tsu.size());
                yu.reserve(tsu.size());
                for (size_t i = 0; i < tsu.size(); ++i) {
                    tu.push_back(tsu.getTime(i));
                    yu.push_back(tsu.getValue(i));
                }

                y_aligned = resample_col_linear(tu, yu, t_base);
            }

            TimeSeries<double> ts;
            ts.setName(prefix + "_" + names[j]);
            ts.reserve(t_base.size());
            for (size_t i = 0; i < t_base.size(); ++i) ts.append(t_base[i], y_aligned[i]);

            if (!out.append(ts, ts.name())) {
                std::cerr << "WARNING: duplicate series skipped: " << ts.name() << "\n";
            } else {
                any = true;
            }
        }

        return any;
    };

    // ============================================================
    // BTC Compare
    // ============================================================
    TimeSeriesSet<double> BTC_compare;
    TimeSeriesSet<double> BTC_fineMean_only;

    for (auto& pr : fine_folders) {
        int rr = pr.first;
        std::string fine_dir = pr.second;
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        std::string pfx  = makeRealLabel(rr) + "_";
        std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");

        load_csv_as_tsset(path, "Fine_" + makeRealLabel(rr), BTC_compare);
    }

    // FineMean for BTC
    {
        bool initialized = false;
        std::vector<std::string> mean_names;

        if (!use_timeseriesset_mean) {
            std::vector<std::vector<double>> sum_cols;
            std::vector<std::vector<int>>    cnt_cols;
            int used = 0;

            for (auto& pr : fine_folders) {
                int rr = pr.first;
                std::string fine_dir = pr.second;
                if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

                std::string pfx  = makeRealLabel(rr) + "_";
                std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");

                std::vector<double> t;
                std::vector<std::string> names;
                std::vector<std::vector<double>> cols;
                if (!read_time_series_table_csv(path, t, names, cols)) continue;

                std::vector<std::vector<double>> cols_aligned(cols.size());
                for (size_t j = 0; j < cols.size(); ++j) {
                    if (align_mode == AlignMode::Resample) {
                        cols_aligned[j] = resample_col_linear(t, cols[j], t_base);
                    } else {
                        TimeSeries<double> ts0;
                        ts0.reserve(t.size());
                        for (size_t i = 0; i < t.size(); ++i) ts0.append(t[i], cols[j][i]);
                        TimeSeries<double> tsu = ts0.make_uniform(dt_cmp, false);

                        std::vector<double> tu, yu;
                        tu.reserve(tsu.size()); yu.reserve(tsu.size());
                        for (size_t i = 0; i < tsu.size(); ++i) {
                            tu.push_back(tsu.getTime(i));
                            yu.push_back(tsu.getValue(i));
                        }
                        cols_aligned[j] = resample_col_linear(tu, yu, t_base);
                    }
                }

                if (!initialized) {
                    mean_names = names;
                    sum_cols.assign(mean_names.size(), std::vector<double>{});
                    cnt_cols.assign(mean_names.size(), std::vector<int>{});
                    initialized = true;
                }
                if (names.size() != mean_names.size()) continue;

                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;
                    accumulate_sum_count(cols_aligned[j], sum_cols[j], cnt_cols[j]);
                }
                used++;
            }

            if (initialized && used > 0) {
                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;

                    std::vector<double> col_mean = finalize_mean_vec(sum_cols[j], cnt_cols[j]);

                    TimeSeries<double> ts;
                    ts.setName("FineMean_" + mean_names[j]);
                    ts.reserve(t_base.size());
                    for (size_t i = 0; i < t_base.size(); ++i) ts.append(t_base[i], col_mean[i]);

                    BTC_compare.append(ts, ts.name());
                    BTC_fineMean_only.append(ts, ts.name());
                }
                std::cout << "Added FineMean (BTC) by sum/count from " << used << " realizations.\n";
            } else {
                std::cerr << "WARNING: FineMean (BTC) not added.\n";
            }

        } else {
            std::vector<TimeSeriesSet<double>> stacks;
            int used = 0;

            for (auto& pr : fine_folders) {
                int rr = pr.first;
                std::string fine_dir = pr.second;
                if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

                std::string pfx  = makeRealLabel(rr) + "_";
                std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled.csv");

                std::vector<double> t;
                std::vector<std::string> names;
                std::vector<std::vector<double>> cols;
                if (!read_time_series_table_csv(path, t, names, cols)) continue;

                std::vector<std::vector<double>> cols_aligned(cols.size());
                for (size_t j = 0; j < cols.size(); ++j) {
                    if (align_mode == AlignMode::Resample) {
                        cols_aligned[j] = resample_col_linear(t, cols[j], t_base);
                    } else {
                        TimeSeries<double> ts0;
                        ts0.reserve(t.size());
                        for (size_t i = 0; i < t.size(); ++i) ts0.append(t[i], cols[j][i]);
                        TimeSeries<double> tsu = ts0.make_uniform(dt_cmp, false);

                        std::vector<double> tu, yu;
                        tu.reserve(tsu.size()); yu.reserve(tsu.size());
                        for (size_t i = 0; i < tsu.size(); ++i) {
                            tu.push_back(tsu.getTime(i));
                            yu.push_back(tsu.getValue(i));
                        }
                        cols_aligned[j] = resample_col_linear(tu, yu, t_base);
                    }
                }

                if (!initialized) {
                    mean_names = names;
                    stacks.assign(mean_names.size(), TimeSeriesSet<double>{});
                    initialized = true;
                }
                if (names.size() != mean_names.size()) continue;

                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;

                    TimeSeries<double> ts;
                    ts.setName(makeRealLabel(rr));
                    ts.reserve(t_base.size());
                    for (size_t i = 0; i < t_base.size(); ++i) ts.append(t_base[i], cols_aligned[j][i]);
                    stacks[j].append(ts, ts.name());
                }

                used++;
            }

            if (initialized && used > 0) {
                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;

                    TimeSeries<double> m = stacks[j].mean_ts(0);

                    TimeSeries<double> ts;
                    ts.setName("FineMean_" + mean_names[j]);
                    ts.reserve(t_base.size());
                    for (size_t i = 0; i < t_base.size(); ++i) {
                        ts.append(t_base[i], (i < m.size()) ? m.getValue(i) : 0.0);
                    }

                    BTC_compare.append(ts, ts.name());
                    BTC_fineMean_only.append(ts, ts.name());
                }

                std::cout << "Added FineMean (BTC) using TimeSeriesSet::mean_ts from "
                          << used << " realizations.\n";
            } else {
                std::cerr << "WARNING: FineMean (BTC) not added.\n";
            }
        }
    }

    // Add upscaled BTCs
    load_csv_as_tsset(up_btc_path, "Upscaled", BTC_compare);

    // Write compare CSVs
    const std::string out_cmp = joinPath(run_dir, "BTC_Compare.csv");
    BTC_compare.write(out_cmp);
    std::cout << "Wrote BTC WIDE comparison CSV (t_base): " << out_cmp << "\n";

    const std::string fm_csv = joinPath(run_dir, "BTC_FineMean.csv");
    if (BTC_fineMean_only.size() > 0) BTC_fineMean_only.write(fm_csv);

    // GNUPLOT (outputs forced into run_dir by passing absolute fig_prefix)
    {
        const std::string gp1 = joinPath(run_dir, "plot_BTC_compare.gp");
        const std::string fig_prefix = joinPath(run_dir, "BTC_compare");
        if (write_btc_compare_plot_gnuplot_by_basename(gp1, out_cmp, fig_prefix, "Concentration")) {
            int rc = run_gnuplot_script(gp1);
            if (rc != 0) std::cerr << "WARNING: gnuplot failed (rc=" << rc << ") for: " << gp1 << "\n";
        }
    }

    // ============================================================
    // Derivative Compare (kept as in your block; enable/extend as needed)
    // ============================================================
    TimeSeriesSet<double> BTCd_compare;
    TimeSeriesSet<double> BTCd_fineMean_only;

    for (auto& pr : fine_folders) {
        int rr = pr.first;
        std::string fine_dir = pr.second;
        if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

        std::string pfx  = makeRealLabel(rr) + "_";
        std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

        load_csv_as_tsset(path, "FineDeriv_" + makeRealLabel(rr), BTCd_compare);
    }

    // FineDerivMean
    {
        bool initialized = false;
        std::vector<std::string> mean_names;

        if (!use_timeseriesset_mean) {
            std::vector<std::vector<double>> sum_cols;
            std::vector<std::vector<int>>    cnt_cols;
            int used = 0;

            for (auto& pr : fine_folders) {
                int rr = pr.first;
                std::string fine_dir = pr.second;
                if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

                std::string pfx  = makeRealLabel(rr) + "_";
                std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

                std::vector<double> t;
                std::vector<std::string> names;
                std::vector<std::vector<double>> cols;
                if (!read_time_series_table_csv(path, t, names, cols)) continue;

                std::vector<std::vector<double>> cols_aligned(cols.size());
                for (size_t j = 0; j < cols.size(); ++j) {
                    if (align_mode == AlignMode::Resample) {
                        cols_aligned[j] = resample_col_linear(t, cols[j], t_base);
                    } else {
                        TimeSeries<double> ts0;
                        ts0.reserve(t.size());
                        for (size_t i = 0; i < t.size(); ++i) ts0.append(t[i], cols[j][i]);
                        TimeSeries<double> tsu = ts0.make_uniform(dt_cmp, false);

                        std::vector<double> tu, yu;
                        tu.reserve(tsu.size()); yu.reserve(tsu.size());
                        for (size_t i = 0; i < tsu.size(); ++i) {
                            tu.push_back(tsu.getTime(i));
                            yu.push_back(tsu.getValue(i));
                        }
                        cols_aligned[j] = resample_col_linear(tu, yu, t_base);
                    }
                }

                if (!initialized) {
                    mean_names = names;
                    sum_cols.assign(mean_names.size(), std::vector<double>{});
                    cnt_cols.assign(mean_names.size(), std::vector<int>{});
                    initialized = true;
                }
                if (names.size() != mean_names.size()) continue;

                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;
                    accumulate_sum_count(cols_aligned[j], sum_cols[j], cnt_cols[j]);
                }
                used++;
            }

            if (initialized && used > 0) {
                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;

                    std::vector<double> col_mean = finalize_mean_vec(sum_cols[j], cnt_cols[j]);

                    TimeSeries<double> ts;
                    ts.setName("FineDerivMean_" + mean_names[j]);
                    ts.reserve(t_base.size());
                    for (size_t i = 0; i < t_base.size(); ++i) ts.append(t_base[i], col_mean[i]);

                    BTCd_compare.append(ts, ts.name());
                    BTCd_fineMean_only.append(ts, ts.name());
                }
                std::cout << "Added FineDerivMean by sum/count from " << used << " realizations.\n";
            } else {
                std::cerr << "WARNING: FineDerivMean not added.\n";
            }

        } else {
            std::vector<TimeSeriesSet<double>> stacks;
            int used = 0;

            for (auto& pr : fine_folders) {
                int rr = pr.first;
                std::string fine_dir = pr.second;
                if (!fine_dir.empty() && fine_dir.back() != '/' && fine_dir.back() != '\\') fine_dir += "/";

                std::string pfx  = makeRealLabel(rr) + "_";
                std::string path = joinPath(fine_dir, pfx + "BTC_FineScaled_derivative.csv");

                std::vector<double> t;
                std::vector<std::string> names;
                std::vector<std::vector<double>> cols;
                if (!read_time_series_table_csv(path, t, names, cols)) continue;

                std::vector<std::vector<double>> cols_aligned(cols.size());
                for (size_t j = 0; j < cols.size(); ++j) {
                    if (align_mode == AlignMode::Resample) {
                        cols_aligned[j] = resample_col_linear(t, cols[j], t_base);
                    } else {
                        TimeSeries<double> ts0;
                        ts0.reserve(t.size());
                        for (size_t i = 0; i < t.size(); ++i) ts0.append(t[i], cols[j][i]);
                        TimeSeries<double> tsu = ts0.make_uniform(dt_cmp, false);

                        std::vector<double> tu, yu;
                        tu.reserve(tsu.size()); yu.reserve(tsu.size());
                        for (size_t i = 0; i < tsu.size(); ++i) {
                            tu.push_back(tsu.getTime(i));
                            yu.push_back(tsu.getValue(i));
                        }
                        cols_aligned[j] = resample_col_linear(tu, yu, t_base);
                    }
                }

                if (!initialized) {
                    mean_names = names;
                    stacks.assign(mean_names.size(), TimeSeriesSet<double>{});
                    initialized = true;
                }
                if (names.size() != mean_names.size()) continue;

                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;

                    TimeSeries<double> ts;
                    ts.setName(makeRealLabel(rr));
                    ts.reserve(t_base.size());
                    for (size_t i = 0; i < t_base.size(); ++i) ts.append(t_base[i], cols_aligned[j][i]);
                    stacks[j].append(ts, ts.name());
                }

                used++;
            }

            if (initialized && used > 0) {
                for (size_t j = 0; j < mean_names.size(); ++j) {
                    if (mean_names[j] == "t" || mean_names[j] == "T" || mean_names[j] == "time") continue;

                    TimeSeries<double> m = stacks[j].mean_ts(0);

                    TimeSeries<double> ts;
                    ts.setName("FineDerivMean_" + mean_names[j]);
                    ts.reserve(t_base.size());
                    for (size_t i = 0; i < t_base.size(); ++i) {
                        ts.append(t_base[i], (i < m.size()) ? m.getValue(i) : 0.0);
                    }

                    BTCd_compare.append(ts, ts.name());
                    BTCd_fineMean_only.append(ts, ts.name());
                }

                std::cout << "Added FineDerivMean using TimeSeriesSet::mean_ts from "
                          << used << " realizations.\n";
            } else {
                std::cerr << "WARNING: FineDerivMean not added.\n";
            }
        }
    }

    load_csv_as_tsset(up_btc_deriv_path, "UpscaledDeriv", BTCd_compare);

    // If you want derivative compare output (and plot), uncomment these lines:
    /*
    const std::string out_cmp_d = joinPath(run_dir, "BTC_Derivative_Compare.csv");
    BTCd_compare.write(out_cmp_d);
    std::cout << "Wrote BTC DERIVATIVE WIDE comparison CSV (t_base): " << out_cmp_d << "\n";

    const std::string fdm_csv = joinPath(run_dir, "BTC_FineDerivMean.csv");
    if (BTCd_fineMean_only.size() > 0) BTCd_fineMean_only.write(fdm_csv);

    {
        const std::string gp2 = joinPath(run_dir, "plot_BTC_derivative_compare.gp");
        const std::string fig_prefix = joinPath(run_dir, "BTC_derivative_compare");
        if (write_btc_compare_plot_gnuplot_by_basename(gp2, out_cmp_d, fig_prefix, "dC/dt")) {
            int rc = run_gnuplot_script(gp2);
            if (rc != 0) std::cerr << "WARNING: gnuplot failed (rc=" << rc << ") for: " << gp2 << "\n";
        }
    }
    */

    return true;
}
