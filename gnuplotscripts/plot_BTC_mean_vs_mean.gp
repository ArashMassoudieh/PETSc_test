# plot_BTC_mean_vs_mean.gp
# Mean vs Mean comparison
# Solid  = 300x100
# Dashed = 150x100
# Dotted = 150x50
# Legend only on x = 0.50

reset

set datafile separator ','
set datafile missing ''

# ------------------------------------------------------------
# Three runs (folders) -- all are in CURRENT DIR
# ------------------------------------------------------------
RUN_300x100     = "100Realizations_20260220_092630_std2_D0.01_aniso1&0.1_df0.15_300x100"
RUN_150x100 = "100Realizations_20260220_133516_std2_D0.01_aniso1&0.1_df0.15_150x100"
RUN_150x50  = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

FILE_300     = sprintf("%s/BTC_mean.csv", RUN_300x100)
FILE_150x100 = sprintf("%s/BTC_mean.csv", RUN_150x100)
FILE_150x50  = sprintf("%s/BTC_mean.csv", RUN_150x50)

# ============================================================
# Plot 1: x = 0.50  (cols 1:2)  <-- legend ON
# ============================================================
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_mean_compare_x0.50_combined.png'

set multiplot layout 1,1

set grid
set key top right
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:10]

plot FILE_300     using 1:2 with lines lw 3 lc rgb "#000000" dt 1 title '300x100', \
     FILE_150x100 using 1:2 with lines lw 3 lc rgb "#000000" dt 2 title '150x100', \
     FILE_150x50  using 1:2 with lines lw 3 lc rgb "#000000" dt 3 title '150x50'

set origin 0.40, 0.15
set size 0.55, 0.65
set key off
unset title
unset xlabel
unset ylabel
set tics font 'Arial,16'
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

plot FILE_300     using 1:2 with lines lw 2 lc rgb "#000000" dt 1 notitle, \
     FILE_150x100 using 1:2 with lines lw 2 lc rgb "#000000" dt 2 notitle, \
     FILE_150x50  using 1:2 with lines lw 2 lc rgb "#000000" dt 3 notitle

unset multiplot


# ============================================================
# Plot 2: x = 1.50  (cols 3:4)  <-- legend OFF
# ============================================================
reset
set datafile separator ','
set datafile missing ''

RUN_300x100     = "100Realizations_20260220_092630_std2_D0.01_aniso1&0.1_df0.15_300x100"
RUN_150x100 = "100Realizations_20260220_133516_std2_D0.01_aniso1&0.1_df0.15_150x100"
RUN_150x50  = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

FILE_300     = sprintf("%s/BTC_mean.csv", RUN_300x100)
FILE_150x100 = sprintf("%s/BTC_mean.csv", RUN_150x100)
FILE_150x50  = sprintf("%s/BTC_mean.csv", RUN_150x50)

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_mean_compare_x1.50_combined.png'

set multiplot layout 1,1

set grid
set key off
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:10]

plot FILE_300     using 3:4 with lines lw 3 lc rgb "#000000" dt 1 notitle, \
     FILE_150x100 using 3:4 with lines lw 3 lc rgb "#000000" dt 2 notitle, \
     FILE_150x50  using 3:4 with lines lw 3 lc rgb "#000000" dt 3 notitle

set origin 0.40, 0.30
set size 0.55, 0.65
set key off
unset xlabel
unset ylabel
set tics font 'Arial,16'
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

plot FILE_300     using 3:4 with lines lw 2 lc rgb "#000000" dt 1 notitle, \
     FILE_150x100 using 3:4 with lines lw 2 lc rgb "#000000" dt 2 notitle, \
     FILE_150x50  using 3:4 with lines lw 2 lc rgb "#000000" dt 3 notitle

unset multiplot


# ============================================================
# Plot 3: x = 2.50  (cols 5:6)  <-- legend OFF
# ============================================================
reset
set datafile separator ','
set datafile missing ''

RUN_300x100     = "100Realizations_20260220_092630_std2_D0.01_aniso1&0.1_df0.15_300x100"
RUN_150x100 = "100Realizations_20260220_133516_std2_D0.01_aniso1&0.1_df0.15_150x100"
RUN_150x50  = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

FILE_300     = sprintf("%s/BTC_mean.csv", RUN_300x100)
FILE_150x100 = sprintf("%s/BTC_mean.csv", RUN_150x100)
FILE_150x50  = sprintf("%s/BTC_mean.csv", RUN_150x50)

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_mean_compare_x2.50_combined.png'

set multiplot layout 1,1

set grid
set key off
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:10]

plot FILE_300     using 5:6 with lines lw 3 lc rgb "#000000" dt 1 notitle, \
     FILE_150x100 using 5:6 with lines lw 3 lc rgb "#000000" dt 2 notitle, \
     FILE_150x50  using 5:6 with lines lw 3 lc rgb "#000000" dt 3 notitle

set origin 0.40, 0.30
set size 0.55, 0.65
set key off
unset xlabel
unset ylabel
set tics font 'Arial,16'
set logscale y
set format y "10^{%T}"
set yrange [1e-6:*]
set xrange [0:20]

plot FILE_300     using 5:6 with lines lw 2 lc rgb "#000000" dt 1 notitle, \
     FILE_150x100 using 5:6 with lines lw 2 lc rgb "#000000" dt 2 notitle, \
     FILE_150x50  using 5:6 with lines lw 2 lc rgb "#000000" dt 3 notitle

unset multiplot
