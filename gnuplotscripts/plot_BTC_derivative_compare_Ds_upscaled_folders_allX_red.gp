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
set style line 31 lc rgb '#d60000' lw 3.5 dt 1
set style line 32 lc rgb '#d60000' lw 3.5 dt 2
set style line 33 lc rgb '#d60000' lw 3.5 dt 3
set style line 34 lc rgb '#d60000' lw 3.5 dt 4
set style line 35 lc rgb '#d60000' lw 3.5 dt 5
set style line 36 lc rgb '#d60000' lw 3.5 dt 6
set style line 37 lc rgb '#d60000' lw 3.5 dt 7

# ===============================
# Parameters
# ===============================
Dvals = "0 0.001 0.01 0.05 0.1"
Xvals = "0.50 1.50 2.50"

# Helper: strip trailing newline from system() output
stripnl(s) = (strlen(s) > 0 && s[strlen(s):strlen(s)] eq "\n") ? s[1:strlen(s)-1] : s

# Newest run dir for each D (last after sort)
find_run_dir(Dtok) = stripnl(system(sprintf("bash -lc 'ls -1d run_*_D%s_* 2>/dev/null | sort | tail -n 1'", Dtok)))

# Find BTC file for given runDir and X:
# tries common exact paths, then falls back to find+grep.
find_btc_file(runDir, X) = stripnl(system(sprintf("bash -lc 'rd=\"%s\"; x=\"%s\"; \
c1=\"$rd/x=$x/BTC_Compare.csv\"; \
c2=\"$rd/x=$xBTC_Compare.csv\"; \
c3=\"$rd/x=$x/BTC_compare.csv\"; \
c4=\"$rd/x=$x/BTC.csv\"; \
c5=\"$rd/BTC_Compare_x=$x.csv\"; \
for f in \"$c1\" \"$c2\" \"$c3\" \"$c4\" \"$c5\"; do if [ -f \"$f\" ]; then echo \"$f\"; exit 0; fi; done; \
find \"$rd\" -type f \\( -iname \"*BTC*Compare*.csv\" -o -iname \"*BTC*compare*.csv\" \\) 2>/dev/null | grep -F \"x=$x\" | sort | head -n 1'", runDir, X)))

# ===============================
# One plot per location X (all Ds on one plot)
# ===============================
do for [X in Xvals] {

    outfile = sprintf("BTC_upscaled_allD_x%s.png", X)
    set output outfile
    set title sprintf("Upscaled BTC  |  x = %s", X) font "Arial,32"

    plotcmd = ""
    sep = ""
    i = 0

    do for [D in Dvals] {
        i = i + 1

        runDir = find_run_dir(D)
        if (strlen(runDir) == 0) {
            print sprintf("WARNING: No folder found matching 'run_*_D%s_*' (skipping)", D)
        } else {

            datafile = find_btc_file(runDir, X)
            if (strlen(datafile) == 0) {
                print sprintf("WARNING: No BTC file found for D=%s at x=%s under %s", D, X, runDir)
            } else {

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
    }

    if (strlen(plotcmd) > 0) {
        eval("plot " . plotcmd)
    } else {
        print sprintf("WARNING: No valid curves for x=%s (nothing to plot)", X)
        # still generates an empty image; if you want to skip, comment set output above and handle differently
        plot NaN notitle
    }
}

unset output

