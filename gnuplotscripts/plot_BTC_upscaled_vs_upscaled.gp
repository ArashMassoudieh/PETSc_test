# plot_BTC_upscaled_vs_upscaled.gp
# Upscaled vs Upscaled comparison (fixed columns per folder)
# Solid red  = 300x100  (2-col files => using 1:2)
# Dashed red = 150x100  (wide files => using 201:202)
# Dotted red = 150x50   (wide files => using 201:202)
# Legend only on x = 0.50

reset
set datafile separator ','
set datafile missing ''

# ------------------------------------------------------------
# Three runs (folders)
# ------------------------------------------------------------
RUN_300x100 = "100Realizations_20260220_092630_std2_D0.01_aniso1&0.1_df0.15_300x100"
RUN_150x100 = "100Realizations_20260220_133516_std2_D0.01_aniso1&0.1_df0.15_150x100"
RUN_150x50  = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

TAG_300x100 = "300x100"
TAG_150x100 = "150x100"
TAG_150x50  = "150x50"

FILE_300x100(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_300x100, x)
FILE_150x100(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_150x100, x)
FILE_150x50(x)  = sprintf("%s/x=%sBTC_Compare.csv", RUN_150x50,  x)

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
unset logscale y
set yrange [0:*]
set xrange [0:10]

plot FILE_300x100("0.50") every ::1 using 1:2     with lines lw 3 lc rgb "#FF0000" dt 1 title TAG_300x100, \
     FILE_150x100("0.50") every ::1 using 201:202 with lines lw 3 lc rgb "#FF0000" dt 2 title TAG_150x100, \
     FILE_150x50("0.50")  every ::1 using 201:202 with lines lw 3 lc rgb "#FF0000" dt 3 title TAG_150x50

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

plot FILE_300x100("0.50") every ::1 using 1:2     with lines lw 2 lc rgb "#FF0000" dt 1 notitle, \
     FILE_150x100("0.50") every ::1 using 201:202 with lines lw 2 lc rgb "#FF0000" dt 2 notitle, \
     FILE_150x50("0.50")  every ::1 using 201:202 with lines lw 2 lc rgb "#FF0000" dt 3 notitle

unset multiplot


# ============================================================
# Plot 2: x = 1.50  <-- legend OFF
# ============================================================
reset
set datafile separator ','
set datafile missing ''

RUN_300x100 = "100Realizations_20260220_092630_std2_D0.01_aniso1&0.1_df0.15_300x100"
RUN_150x100 = "100Realizations_20260220_133516_std2_D0.01_aniso1&0.1_df0.15_150x100"
RUN_150x50  = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

FILE_300x100(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_300x100, x)
FILE_150x100(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_150x100, x)
FILE_150x50(x)  = sprintf("%s/x=%sBTC_Compare.csv", RUN_150x50,  x)

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_upscaled_compare_x1.50_combined.png'

set multiplot layout 1,1

set grid
set key off
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
unset logscale y
set yrange [0:*]
set xrange [0:10]

plot FILE_300x100("1.50") every ::1 using 1:2     with lines lw 3 lc rgb "#FF0000" dt 1 notitle, \
     FILE_150x100("1.50") every ::1 using 201:202 with lines lw 3 lc rgb "#FF0000" dt 2 notitle, \
     FILE_150x50("1.50")  every ::1 using 201:202 with lines lw 3 lc rgb "#FF0000" dt 3 notitle

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

plot FILE_300x100("1.50") every ::1 using 1:2     with lines lw 2 lc rgb "#FF0000" dt 1 notitle, \
     FILE_150x100("1.50") every ::1 using 201:202 with lines lw 2 lc rgb "#FF0000" dt 2 notitle, \
     FILE_150x50("1.50")  every ::1 using 201:202 with lines lw 2 lc rgb "#FF0000" dt 3 notitle

unset multiplot


# ============================================================
# Plot 3: x = 2.50  <-- legend OFF
# ============================================================
reset
set datafile separator ','
set datafile missing ''

RUN_300x100 = "100Realizations_20260220_092630_std2_D0.01_aniso1&0.1_df0.15_300x100"
RUN_150x100 = "100Realizations_20260220_133516_std2_D0.01_aniso1&0.1_df0.15_150x100"
RUN_150x50  = "100Realizations_20260218_161839_std2_D0.01_aniso1&0.1_df0.15_150x50"

FILE_300x100(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_300x100, x)
FILE_150x100(x) = sprintf("%s/x=%sBTC_Compare.csv", RUN_150x100, x)
FILE_150x50(x)  = sprintf("%s/x=%sBTC_Compare.csv", RUN_150x50,  x)

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'BTC_upscaled_compare_x2.50_combined.png'

set multiplot layout 1,1

set grid
set key off
set xlabel 'Time' font 'Arial,32'
set ylabel 'c/c_0' font 'Arial,32'
set format y "%.1f"
unset logscale y
set yrange [0:*]
set xrange [0:10]

plot FILE_300x100("2.50") every ::1 using 1:2     with lines lw 3 lc rgb "#FF0000" dt 1 notitle, \
     FILE_150x100("2.50") every ::1 using 201:202 with lines lw 3 lc rgb "#FF0000" dt 2 notitle, \
     FILE_150x50("2.50")  every ::1 using 201:202 with lines lw 3 lc rgb "#FF0000" dt 3 notitle

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

plot FILE_300x100("2.50") every ::1 using 1:2     with lines lw 2 lc rgb "#FF0000" dt 1 notitle, \
     FILE_150x100("2.50") every ::1 using 201:202 with lines lw 2 lc rgb "#FF0000" dt 2 notitle, \
     FILE_150x50("2.50")  every ::1 using 201:202 with lines lw 2 lc rgb "#FF0000" dt 3 notitle

unset multiplot
