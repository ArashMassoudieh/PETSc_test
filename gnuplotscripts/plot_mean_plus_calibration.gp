#!/usr/bin/gnuplot

# ===============================
# Terminal & appearance
# ===============================
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set datafile separator ','
set grid
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]

set xlabel 't' font "Arial,32"
set ylabel 'c/c_0' font "Arial,32"
set key right top
set key font "Arial,18"
set key maxcols 4
set key samplen 2
set key spacing 1.1
set key width -2

# ===============================
# Styles
# ===============================
set style line 1  lc rgb '#000000' lw 4 dt 1  # black mean

# red calibration curves (cycle dash types)
set style line 11 lc rgb '#d60000' lw 3 dt 1
set style line 12 lc rgb '#d60000' lw 3 dt 2
set style line 13 lc rgb '#d60000' lw 3 dt 3
set style line 14 lc rgb '#d60000' lw 3 dt 4
set style line 15 lc rgb '#d60000' lw 3 dt 5
set style line 16 lc rgb '#d60000' lw 3 dt 6
set style line 17 lc rgb '#d60000' lw 3 dt 7
set style line 18 lc rgb '#d60000' lw 3 dt 8
set style line 19 lc rgb '#d60000' lw 3 dt 9

# ===============================
# Parameters (edit this path)
# ===============================
Xvals = "0.50 1.50 2.50"

# Script should be run from: .../UpscalingResults/Calibration
CALIB_ROOT = "."

# Resume folder containing BTC_mean.csv (black mean source)
RESUME_DIR = "/mnt/3rd900/Projects/PETSc_test/UpscalingResults/100Realizations_20260202_003241_std2_D0.1_aniso"
BLACK_MEAN = sprintf("%s/BTC_mean.csv", RESUME_DIR)

# ===============================
# Helpers
# ===============================
stripnl(s) = (strlen(s) > 0 && s[strlen(s):strlen(s)] eq "\n") ? s[1:strlen(s)-1] : s

# Map paired columns in BTC_mean.csv:
#  t,x=0.50,t,x=1.50,t,x=2.50,... so:
#  x=0.50 -> using 1:2
#  x=1.50 -> using 3:4
#  x=2.50 -> using 5:6
col_t(x) = (x eq "0.50") ? 1 : (x eq "1.50") ? 3 : 5
col_c(x) = (x eq "0.50") ? 2 : (x eq "1.50") ? 4 : 6

# Extract df from folder name using bash (no gnuplot regex pain)
# returns e.g. "0.15"
df_from_dir(d) = stripnl(system(sprintf("bash -lc \"basename '%s' | sed -n 's/.*_df\\([0-9.]*\\).*/\\1/p'\"", d)))

# ===============================
# Get calibration run directories (string list)
# ===============================
calib_dirs = stripnl(system(sprintf("bash -lc 'ls -1d %s/run_*_df* 2>/dev/null | sort'", CALIB_ROOT)))

if (strlen(calib_dirs) == 0) {
    print "ERROR: No calibration run folders found (expected run_*_df* under Calibration/)."
    print "       Run this script from the Calibration folder, or fix CALIB_ROOT."
}

# ===============================
# One plot per X
# ===============================
do for [X in Xvals] {

    outfile = sprintf("BTC_mean_plus_calibration_x%s.png", X)
    set output outfile
    set title sprintf("BTC Mean (black) + Calibration (red)  |  x = %s", X) font "Arial,32"

    # ---- start plot with black mean ----
    plotcmd = sprintf("'%s' every ::1 using %d:%d with lines ls 1 title 'Mean'", \
                      BLACK_MEAN, col_t(X), col_c(X))

    # ---- add calibration curves ----
    n = words(calib_dirs)
    i = 0

    do for [k=1:n] {
        i = i + 1
        runDir = word(calib_dirs, k)
        df = df_from_dir(runDir)

        datafile = sprintf("%s/x=%sBTC_Compare.csv", runDir, X)

        # Determine last column (Upscaled) dynamically
        stats datafile every ::1 nooutput
        ncol = STATS_columns

        if (ncol >= 2) {
            ls = 11 + ((i-1) % 9)  # cycle 11..19
            plotcmd = plotcmd . sprintf(", '%s' every ::1 using 1:%d with lines ls %d title 'Upscaled df=%s'", \
                                        datafile, ncol, ls, df)
        } else {
            print sprintf("WARNING: Bad column count in %s (ncol=%d)", datafile, ncol)
        }
    }

    eval("plot " . plotcmd)
}

unset output
