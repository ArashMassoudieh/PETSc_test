# plot_BTC_combined.gp
# Generates plots with linear main plot and log inset

reset

# Common data settings
set datafile separator ','
set datafile missing ''

# ============================================================
# Plot 1: x = 0.50
# ============================================================
set output 'BTC_x0.50_combined.png'

set multiplot layout 1,1

# --- Main plot (linear) ---
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set grid
set key top right
set xlabel 'Time' font 'Arial,32'
set ylabel 'Concentration' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:5]
set title 'Breakthrough Curve at x = 0.50 m' font 'Arial,32'

plot 'x=0.50BTC_Compare.csv' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" title 'Realizations', \
     'BTC_mean.csv' using 1:2 with lines lw 3 lc rgb "#000000" title 'Mean', \
     'x=0.50BTC_Compare.csv' using 41:42 with lines lw 3 lc rgb "#FF0000" title 'Upscaled'

# --- Inset plot (log) ---
set origin 0.55, 0.50     # Position: right side, top half
set size 0.40, 0.45       # Size: 40% width, 45% height
set grid
set key off               # No legend in inset
unset title               # No title in inset
set xlabel 'Time' font 'Arial,18'
set ylabel 'Conc. (log)' font 'Arial,18'
set tics font 'Arial,16'
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

plot 'x=0.50BTC_Compare.csv' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     'BTC_mean.csv' using 1:2 with lines lw 2 lc rgb "#000000" notitle, \
     'x=0.50BTC_Compare.csv' using 41:42 with lines lw 2 lc rgb "#FF0000" notitle

unset multiplot

# ============================================================
# Plot 2: x = 1.50 (similar structure)
# ============================================================
reset
set datafile separator ','
set datafile missing ''
set output 'BTC_x1.50_combined.png'

set multiplot layout 1,1

# Main plot
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set grid
set key top right
set xlabel 'Time' font 'Arial,32'
set ylabel 'Concentration' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:5]
set title 'Breakthrough Curve at x = 1.50 m' font 'Arial,32'

plot 'x=1.50BTC_Compare.csv' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" title 'Realizations', \
     'BTC_mean.csv' using 3:4 with lines lw 3 lc rgb "#000000" title 'Mean', \
     'x=1.50BTC_Compare.csv' using 41:42 with lines lw 3 lc rgb "#FF0000" title 'Upscaled'

# Inset
set origin 0.55, 0.50
set size 0.40, 0.45
set grid
set key off
unset title
set xlabel 'Time' font 'Arial,18'
set ylabel 'Conc. (log)' font 'Arial,18'
set tics font 'Arial,16'
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

plot 'x=1.50BTC_Compare.csv' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     'BTC_mean.csv' using 3:4 with lines lw 2 lc rgb "#000000" notitle, \
     'x=1.50BTC_Compare.csv' using 41:42 with lines lw 2 lc rgb "#FF0000" notitle

unset multiplot

# ============================================================
# Plot 3: x = 2.50
# ============================================================
reset
set datafile separator ','
set datafile missing ''
set output 'BTC_x2.50_combined.png'

set multiplot layout 1,1

# Main plot
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set grid
set key top right
set xlabel 'Time' font 'Arial,32'
set ylabel 'Concentration' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:5]
set title 'Breakthrough Curve at x = 2.50 m' font 'Arial,32'

plot 'x=2.50BTC_Compare.csv' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" title 'Realizations', \
     'BTC_mean.csv' using 5:6 with lines lw 3 lc rgb "#000000" title 'Mean', \
     'x=2.50BTC_Compare.csv' using 41:42 with lines lw 3 lc rgb "#FF0000" title 'Upscaled'

# Inset
set origin 0.55, 0.50
set size 0.40, 0.45
set grid
set key off
unset title
set xlabel 'Time' font 'Arial,18'
set ylabel 'Conc. (log)' font 'Arial,18'
set tics font 'Arial,16'
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

plot 'x=2.50BTC_Compare.csv' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     'BTC_mean.csv' using 5:6 with lines lw 2 lc rgb "#000000" notitle, \
     'x=2.50BTC_Compare.csv' using 41:42 with lines lw 2 lc rgb "#FF0000" notitle

unset multiplot
