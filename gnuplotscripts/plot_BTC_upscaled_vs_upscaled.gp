# plot_BTC_upscaled_vs_upscaled.gp
# Upscaled vs Upscaled comparison
# Solid red  = RUN_A (e.g., 300x100)
# Dashed red = RUN_B (e.g., 150x50)
# Legend only on x = 0.50
# Uses x=<x>BTC_Compare.csv columns 201:202 (Upscaled)

reset
set datafile separator ','
set datafile missing ''

# ------------------------------------------------------------
# Two runs (folders)
# ------------------------------------------------------------
RUN_A = "100Realizations_std2_D0.01_aniso1&0.1"
RUN_B = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

TAG_A = "300x100"
TAG_B = "150x50"

FILE_A(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_A, x)
FILE_B(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_B, x)

# ============================================================
# Plot 1: x = 0.50  <-- legend ON
# ============================================================
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_upscaled_compare_x0.50_combined.png'

set multiplot layout 1,1

set grid
set key top right
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:10]

plot FILE_A("0.50") using 201:202 with lines lw 3 lc rgb "#FF0000" dt 1 title TAG_A, \
     FILE_B("0.50") using 201:202 with lines lw 3 lc rgb "#FF0000" dt 2 title TAG_B

# --- Inset (log) ---
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

plot FILE_A("0.50") using 201:202 with lines lw 2 lc rgb "#FF0000" dt 1 notitle, \
     FILE_B("0.50") using 201:202 with lines lw 2 lc rgb "#FF0000" dt 2 notitle

unset multiplot


# ============================================================
# Plot 2: x = 1.50  <-- legend OFF
# ============================================================
reset
set datafile separator ','
set datafile missing ''

RUN_A = "100Realizations_std2_D0.01_aniso1&0.1"
RUN_B = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"
TAG_A = "300x100"
TAG_B = "150x50"
FILE_A(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_A, x)
FILE_B(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_B, x)

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_upscaled_compare_x1.50_combined.png'

set multiplot layout 1,1

set grid
set key off
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:10]

plot FILE_A("1.50") using 201:202 with lines lw 3 lc rgb "#FF0000" dt 1 notitle, \
     FILE_B("1.50") using 201:202 with lines lw 3 lc rgb "#FF0000" dt 2 notitle

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

plot FILE_A("1.50") using 201:202 with lines lw 2 lc rgb "#FF0000" dt 1 notitle, \
     FILE_B("1.50") using 201:202 with lines lw 2 lc rgb "#FF0000" dt 2 notitle

unset multiplot


# ============================================================
# Plot 3: x = 2.50  <-- legend OFF
# ============================================================
reset
set datafile separator ','
set datafile missing ''

RUN_A = "100Realizations_std2_D0.01_aniso1&0.1"
RUN_B = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"
TAG_A = "300x100"
TAG_B = "150x50"
FILE_A(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_A, x)
FILE_B(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_B, x)

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_upscaled_compare_x2.50_combined.png'

set multiplot layout 1,1

set grid
set key off
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
set yrange [0:*]
set xrange [0:10]

plot FILE_A("2.50") using 201:202 with lines lw 3 lc rgb "#FF0000" dt 1 notitle, \
     FILE_B("2.50") using 201:202 with lines lw 3 lc rgb "#FF0000" dt 2 notitle

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

plot FILE_A("2.50") using 201:202 with lines lw 2 lc rgb "#FF0000" dt 1 notitle, \
     FILE_B("2.50") using 201:202 with lines lw 2 lc rgb "#FF0000" dt 2 notitle

unset multiplot
