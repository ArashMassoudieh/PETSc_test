// plotter.h
#pragma once

#include <string>
#include <vector>

// -------------------------------------------------------------
// Compare-time switches
// -------------------------------------------------------------
enum class TBaseMode { Fixed, FromUpscaled, FromFirstFine };
enum class AlignMode { Resample, MakeUniform };

// Run final aggregation + (optional) gnuplot output.
// - Writes BTC_Compare.csv (+ optional BTC_FineMean.csv)
// - Writes derivative compare CSV if you enable that section
// - Writes gnuplot scripts + runs them (as in your block)
//
// Returns true if it ran (even if some series were missing).
// Returns false only on "hard" failures (e.g., cannot list fine folders when required).
bool run_final_aggregation_and_plots(
    const std::string& run_dir,
    const std::string& up_btc_path,
    const std::string& up_btc_deriv_path,
    TBaseMode tbase_mode,
    AlignMode align_mode,
    bool use_timeseriesset_mean,
    double t_end_cmp = 10.0,
    double dt_cmp    = 0.001
);
