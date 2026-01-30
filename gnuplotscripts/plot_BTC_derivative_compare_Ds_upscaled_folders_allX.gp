#!/usr/bin/gnuplot

# ===============================
# Terminal & appearance
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
# Line styles (all red, different dash types per D)
# ===============================
set style line 31 lc rgb '#000000' lw 3.5 dt 1
set style line 32 lc rgb '#000000' lw 3.5 dt 2
set style line 33 lc rgb '#000000' lw 3.5 dt 3
set style line 34 lc rgb '#000000' lw 3.5 dt 4
set style line 35 lc rgb '#000000' lw 3.5 dt 5
set style line 36 lc rgb '#000000' lw 3.5 dt 6
set style line 37 lc rgb '#000000' lw 3.5 dt 7

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
# ===============================
do for [X in Xvals] {

    outfile = sprintf("BTC_upscaled_allD_x%s.png", X)
    set output outfile
    set title sprintf("Upscaled BTC  |  x = %s", X) font "Arial,32"

    # Build plot command by looping D and appending with commas
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

            # last column = Upscaled; skip header row
            stats datafile every ::1 nooutput
            ncol = STATS_columns

            if (ncol >= 2) {
                plotcmd = plotcmd . sep . sprintf("'%s' every ::1 using 1:%d with lines ls %d title 'D=%s'", datafile, ncol, 30+i, D)
                sep = ", "
            } else {
                print sprintf("WARNING: Bad column count in %s (ncol=%d)", datafile, ncol)
            }
        }
    }

    if (strlen(plotcmd) > 0) {
        eval("plot " . plotcmd)
    }
}

unset output
