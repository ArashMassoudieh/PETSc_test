# File overview: automatic BTC compare plotting for PETSc_test outputs.
# plot_BTC_both_100_auto.gp
#
# Purpose:
#   Auto-find fine_r* folders and generate BTC comparison plots for any
#   x=*BTC_Compare.csv files found inside each folder.
#
# For each compare file:
#   - fine realizations = all column pairs except final pair
#   - upscaled curve    = final column pair
#   - mean curve        = BTC_mean.csv if present (local folder, else parent)
#   - linear main panel + log-y inset
#
# Usage:
#   cd <run_folder_containing_fine_r0001,...>
#   gnuplot /path/to/gnuplotscripts/plot_BTC_both_100_auto.gp
#
# Optional overrides:
#   out_prefix, x_main_max, x_log_max, y_log_min

reset
set datafile separator ','
set datafile missing ''
set datafile columnheaders

if (!exists("out_prefix")) out_prefix = "BTC_compare"
if (!exists("x_main_max")) x_main_max = 10
if (!exists("x_log_max"))  x_log_max  = 20
if (!exists("y_log_min"))  y_log_min  = 1e-6

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'

root_dir = GPVAL_PWD

# Prefer fine_r* folders; if none exist, fall back to current folder.
dir_list = system("find . -maxdepth 1 -type d -name 'fine_r*' | sort")
if (strlen(dir_list) == 0) dir_list = "."

n_dirs = words(dir_list)
print sprintf("Found %d plotting target folder(s).", n_dirs)

# -------- helper-like plotting block via repeated code --------
do for [di=1:n_dirs] {
    rel_dir = word(dir_list, di)

    # Enter folder
    eval sprintf("cd '%s'", rel_dir)
    print sprintf("Processing folder: %s", GPVAL_PWD)

    files = system("ls x=*BTC_Compare.csv 2>/dev/null")
    if (strlen(files) == 0) {
        print "  No x=*BTC_Compare.csv found; trying fine-vs-upscaled fallback."

        fine_file = word(system("ls r*_BTC_FineScaled.csv 2>/dev/null"), 1)
        up_file = "../upscaled_mean/upscaled_BTC_Upscaled.csv"

        if (strlen(fine_file) == 0) {
            print "  No r*_BTC_FineScaled.csv found; skipping."
            eval sprintf("cd '%s'", root_dir)
            continue
        }
        if (system(sprintf("test -f %s", up_file)) ne "") {
            print sprintf("  Missing upscaled reference file: %s ; skipping.", up_file)
            eval sprintf("cd '%s'", root_dir)
            continue
        }

        stats fine_file nooutput
        ncol_f = STATS_columns
        npairs_f = int(ncol_f/2)
        stats up_file nooutput
        ncol_u = STATS_columns
        npairs_u = int(ncol_u/2)
        npairs_fb = (npairs_f < npairs_u) ? npairs_f : npairs_u
        if (npairs_fb < 1) {
            print sprintf("  Skipping fallback: no column pairs in %s", fine_file)
            eval sprintf("cd '%s'", root_dir)
            continue
        }

        do for [p=1:npairs_fb] {
            ft_col = 2*p - 1
            fc_col = 2*p
            ut_col = ft_col
            uc_col = fc_col

            stem = fine_file[1:strlen(fine_file)-4]
            out_png = sprintf("%s_%s_pair%02d_auto.png", out_prefix, stem, p)
            set output out_png
            set multiplot layout 1,1

            set grid
            set key top right
            set xlabel 'Time' font 'Arial,32'
            set ylabel 'c/c_0' font 'Arial,32'
            set format y "%.1f"
            set yrange [0:*]
            set xrange [0:x_main_max]
            unset logscale y
            unset y2tics
            set origin 0,0
            set size 1,1

            plot fine_file using ft_col:fc_col with lines lw 2 lc rgb "#000000" title sprintf("Fine pair %d", p), \
                 up_file   using ut_col:uc_col with lines lw 3 lc rgb "#FF0000" title "Upscaled"

            set origin 0.40, 0.15
            set size 0.55, 0.65
            set key off
            unset title
            unset xlabel
            unset ylabel
            set tics font 'Arial,16'
            set logscale y
            set format y "10^{%T}"
            set yrange [y_log_min:*]
            set xrange [0:x_log_max]
            unset y2tics

            plot fine_file using ft_col:fc_col with lines lw 2 lc rgb "#000000" notitle, \
                 up_file   using ut_col:uc_col with lines lw 3 lc rgb "#FF0000" notitle

            unset multiplot
            print sprintf("  Wrote fallback %s (pair=%d)", out_png, p)
        }

        eval sprintf("cd '%s'", root_dir)
        continue
    }

    # Local mean first, then parent folder fallback
    mean_file = "BTC_mean.csv"
    if (system("test -f BTC_mean.csv") ne "") {
        if (system("test -f ../BTC_mean.csv") eq "") {
            mean_file = "../BTC_mean.csv"
        }
    }

    n_files = words(files)
    do for [i=1:n_files] {
        f = word(files, i)

        stats f nooutput
        ncol = STATS_columns
        if (ncol < 4) {
            print sprintf("  Skipping %s: expected >=4 columns, got %d", f, ncol)
            continue
        }

        npairs = int(ncol/2)
        fine_pairs = npairs - 1
        up_t_col = ncol - 1
        up_c_col = ncol

        stem = f[1:strlen(f)-4]
        out_png = sprintf("%s_%s_auto.png", out_prefix, stem)
        set output out_png
        set multiplot layout 1,1

        # Main (linear)
        set grid
        set key top right
        set xlabel 'Time' font 'Arial,32'
        set ylabel 'c/c_0' font 'Arial,32'
        set format y "%.1f"
        set yrange [0:*]
        set xrange [0:x_main_max]
        unset logscale y
        set origin 0,0
        set size 1,1

        has_mean = (system(sprintf("test -f %s", mean_file)) eq "") ? 1 : 0

        if (fine_pairs >= 1) {
            if (has_mean) {
                plot for [k=1:fine_pairs] f using (column(2*k-1)):(column(2*k)) with lines lw 1 lc rgb "#CCCCCC" notitle, \
                     mean_file using 1:2 with lines lw 3 lc rgb "#000000" title 'Mean', \
                     f using up_t_col:up_c_col with lines lw 3 lc rgb "#FF0000" title 'Upscaled'
            } else {
                plot for [k=1:fine_pairs] f using (column(2*k-1)):(column(2*k)) with lines lw 1 lc rgb "#CCCCCC" notitle, \
                     f using up_t_col:up_c_col with lines lw 3 lc rgb "#FF0000" title 'Upscaled'
            }
        } else {
            plot f using up_t_col:up_c_col with lines lw 3 lc rgb "#FF0000" title 'Upscaled'
        }

        # Inset (log-y)
        set origin 0.40, 0.15
        set size 0.55, 0.65
        set key off
        unset title
        unset xlabel
        unset ylabel
        set tics font 'Arial,16'
        set logscale y
        set format y "10^{%T}"
        set yrange [y_log_min:*]
        set xrange [0:x_log_max]

        if (fine_pairs >= 1) {
            if (has_mean) {
                plot for [k=1:fine_pairs] f using (column(2*k-1)):(column(2*k)) with lines lw 1 lc rgb "#CCCCCC" notitle, \
                     mean_file using 1:2 with lines lw 3 lc rgb "#000000" notitle, \
                     f using up_t_col:up_c_col with lines lw 3 lc rgb "#FF0000" notitle
            } else {
                plot for [k=1:fine_pairs] f using (column(2*k-1)):(column(2*k)) with lines lw 1 lc rgb "#CCCCCC" notitle, \
                     f using up_t_col:up_c_col with lines lw 3 lc rgb "#FF0000" notitle
            }
        } else {
            plot f using up_t_col:up_c_col with lines lw 3 lc rgb "#FF0000" notitle
        }

        unset multiplot
        print sprintf("  Wrote %s (fine_pairs=%d, upscaled cols=%d:%d)", out_png, fine_pairs, up_t_col, up_c_col)
    }

    # Return to root for next folder
    eval sprintf("cd '%s'", root_dir)
}

print "Done."
