# plot_BTC_compare_all.gp
# Generates three separate plots for x=0.50, x=1.50, x=2.50

reset

# Common settings
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set grid
set key top right
set xlabel 'Time' font 'Arial,32'
set ylabel 'Concentration' font 'Arial,32'
set datafile separator ','
set datafile missing ''

# Logarithmic y-axis with superscript notation
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

# ============================================================
# Plot 1: x = 0.50
# ============================================================
set output 'BTC_compare_x0.50.png'
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

# ============================================================
# Plot 2: x = 1.50
# ============================================================
set output 'BTC_compare_x1.50.png'
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

# ============================================================
# Plot 3: x = 2.50
# ============================================================
set output 'BTC_compare_x2.50.png'
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
