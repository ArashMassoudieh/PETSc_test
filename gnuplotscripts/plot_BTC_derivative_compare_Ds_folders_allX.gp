#!/usr/bin/gnuplot

# ===============================
# Terminal & general appearance
# ===============================
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set datafile separator ','
set grid
set logscale y
set format y "10^{%T}"
set yrange [1e-4:*]

set xlabel 't' font "Arial,32"
set ylabel 'c/c_0' font "Arial,32"
set key right top

# ===============================
# Line styles
# ===============================
# Mean (black) — different dash types per D
set style line 41 lc rgb '#000000' lw 3.0 dt 1
set style line 42 lc rgb '#000000' lw 3.0 dt 2
set style line 43 lc rgb '#000000' lw 3.0 dt 3
set style line 44 lc rgb '#000000' lw 3.0 dt 4
set style line 45 lc rgb '#000000' lw 3.0 dt 5
set style line 46 lc rgb '#000000' lw 3.0 dt 6
set style line 47 lc rgb '#000000' lw 3.0 dt 7

# Upscaled (red) — different dash types per D
set style line 31 lc rgb '#d60000' lw 3.0 dt 1
set style line 32 lc rgb '#d60000' lw 3.0 dt 2
set style line 33 lc rgb '#d60000' lw 3.0 dt 3
set style line 34 lc rgb '#d60000' lw 3.0 dt 4
set style line 35 lc rgb '#d60000' lw 3.0 dt 5
set style line 36 lc rgb '#d60000' lw 3.0 dt 6
set style line 37 lc rgb '#d60000' lw 3.0 dt 7

# ===============================
# Parameters
# ===============================
Dvals = "0 0.01 0.02 0.05 0.1 0.2 0.5"
Xvals = "0.50 1.50 2.50"

# Helper: strip trailing newline from system() output
stripnl(s) = (strlen(s) > 0 && s[strlen(s):strlen(s)] eq "\n") ? s[1:strlen(s)-1] : s

# Find newest folder matching "*_D<D>"
find_run_dir(D) = stripnl(system(sprintf("bash -lc 'ls -td *_D%s 2>/dev/null | head -n 1'", D)))

# ===============================
# One plot per location X (7 Ds on one plot)
# Plots Mean (black) and Upscaled (red) for each D
# ===============================
do for [X in Xvals] {

    outfile = sprintf("BTC_mean_upscaled_allD_x%s_log.png", X)
    set output outfile
    set title sprintf("BTC Mean vs Upscaled  |  x = %s", X) font "Arial,32"

    plotcmd = ""
    sep = ""
    i = 0

    do for [D in Dvals] {
        i = i + 1
        runDir = find_run_dir(D)

        if (strlen(runDir) == 0) {
            print sprintf("WARNING: No folder found matching '*_D%s' (skipping)", D)
        }

        if (strlen(runDir) > 0) {
            datafile = sprintf("%s/x=%sBTC_Compare.csv", runDir, X)

            stats datafile every ::1 nooutput
            ncol = STATS_columns

            if (ncol >= 4) {
                # Mean (ncol-1)
                plotcmd = plotcmd . sep . sprintf("'%s' every ::1 using 1:%d with lines ls %d title 'Mean D=%s'", datafile, ncol-1, 40+i, D)
                sep = ", "
                # Upscaled (ncol)
                plotcmd = plotcmd . sep . sprintf("'%s' every ::1 using 1:%d with lines ls %d title 'Upscaled D=%s'", datafile, ncol, 30+i, D)
                sep = ", "
            } else {
                print sprintf("WARNING: Not enough columns in %s (ncol=%d)", datafile, ncol)
            }
        }
    }

    if (strlen(plotcmd) > 0) {
        eval("plot " . plotcmd)
    }
}

unset output
