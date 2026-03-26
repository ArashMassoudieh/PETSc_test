# File overview: automatic BTC compare plotting for PETSc_test outputs.
# plot_BTC_both_100_auto.gp
#
# Purpose:
#   Auto-plot all x=*BTC_Compare.csv files in current directory.
#   For each file:
#     - Plot all fine-scale realization curves (all column pairs except last pair)
#     - Plot mean curve from BTC_mean.csv (if available)
#     - Plot upscaled curve from the last column pair in x=*BTC_Compare.csv
#     - Add log-scale inset for tail behavior
#
# Usage:
#   cd <run_folder>
#   gnuplot /path/to/gnuplotscripts/plot_BTC_both_100_auto.gp
#
# Optional overrides (gnuplot -e):
#   x_main_max, x_log_max, y_log_min, out_prefix

reset

set datafile separator ','
set datafile missing ''

if (!exists("out_prefix")) out_prefix = "BTC_compare"
if (!exists("x_main_max")) x_main_max = 10
if (!exists("x_log_max"))  x_log_max  = 20
if (!exists("y_log_min"))  y_log_min  = 1e-6

mean_file = "BTC_mean.csv"
files = system("ls x=*BTC_Compare.csv 2>/dev/null")

if (strlen(files) == 0) {
    print "No x=*BTC_Compare.csv files found in current directory."
    exit
}

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'

# Iterate over all compare files
n_files = words(files)
do for [i=1:n_files] {
    f = word(files, i)

    # Read number of columns in compare file
    stats f nooutput
    ncol = STATS_columns

    if (ncol < 4) {
        print sprintf("Skipping %s: expected at least 4 columns, found %d", f, ncol)
        continue
    }

    npairs = int(ncol/2)
    fine_pairs = npairs - 1
    up_t_col = ncol - 1
    up_c_col = ncol

    # Build output name from input file stem
    stem = f[1:strlen(f)-4]
    out_png = sprintf("%s_%s_auto.png", out_prefix, stem)
    set output out_png

    set multiplot layout 1,1

    # ---------- Main plot (linear y) ----------
    set grid
    set key top right
    set xlabel 'Time' font 'Arial,32'
    set ylabel 'c/c_0' font 'Arial,32'
    set format y "%.1f"
    set yrange [0:*]
    set xrange [0:x_main_max]
    unset logscale y

    if (fine_pairs >= 1) {
        if (system(sprintf("test -f %s", mean_file)) eq "") {
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

    # ---------- Inset plot (log y) ----------
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
        if (system(sprintf("test -f %s", mean_file)) eq "") {
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
    print sprintf("Wrote %s from %s (fine pairs=%d, upscaled cols=%d:%d)", out_png, f, fine_pairs, up_t_col, up_c_col)
}
