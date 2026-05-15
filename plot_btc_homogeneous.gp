# =====================================================================
# plot_btc_homogeneous.gp
#
# Three-panel BTC comparison at x = 0.5, 1.5, 2.5 for the
# near-homogeneous control case (sigma_lnK = 0.05, D = 0).
#
# Fine-scale and upscaled should agree closely. Any visible gap is
# diagnostic of a residual numerical/structural issue.
#
# File layout (TimeSeriesSet::write, 3 locations):
#   col 1, 2 : t, c   for x = 0.5
#   col 3, 4 : t, c   for x = 1.5
#   col 5, 6 : t, c   for x = 2.5
#
# Run from the directory containing 'Homogeneous/':
#     gnuplot plot_btc_homogeneous.gp
# =====================================================================

# ---- (t, c) column pairs ----------------------------------------------
t_05 = 1
c_05 = 2
t_15 = 3
c_15 = 4
t_25 = 5
c_25 = 6

# ---- paths ------------------------------------------------------------
real_dir = "Homogeneous"
mean_f   = "Homogeneous/mean_btc.csv"
ups_f    = "Homogeneous/Upscaled/D_0/btc.csv"
n_real   = 2

# ---- output -----------------------------------------------------------
set terminal pngcairo size 1800,600 enhanced font "Arial,16"
set output "btc_homogeneous.png"

# ---- styling ----------------------------------------------------------
set datafile separator ","
set datafile missing "NaN"

set logscale y
set format y "10^{%L}"
set yrange [1e-4:*]
set xlabel "time"
set ylabel "c(t)"
set grid
set key bottom left

set style line 1 lc rgb "gray70" lw 1     # realizations
set style line 2 lc rgb "black"  lw 3     # fine-scale mean
set style line 3 lc rgb "red"    lw 3     # upscaled

set multiplot layout 1,3 \
    title "Homogeneous control: fine-scale ensemble vs. upscaled (sigma_{lnK}=0.05, D=0)" \
    font "Arial,18"

# ===== panel 1: x = 0.5 ================================================
set title "x = 0.5"
tc = t_05
cc = c_05
plot \
    for [i=1:n_real] sprintf("%s/r%d/btc.csv", real_dir, i) \
        using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 1 notitle, \
    mean_f using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 2 title "fine-scale mean", \
    ups_f  using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 3 title "upscaled"

# ===== panel 2: x = 1.5 ================================================
set title "x = 1.5"
tc = t_15
cc = c_15
plot \
    for [i=1:n_real] sprintf("%s/r%d/btc.csv", real_dir, i) \
        using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 1 notitle, \
    mean_f using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 2 title "fine-scale mean", \
    ups_f  using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 3 title "upscaled"

# ===== panel 3: x = 2.5 ================================================
set title "x = 2.5"
tc = t_25
cc = c_25
plot \
    for [i=1:n_real] sprintf("%s/r%d/btc.csv", real_dir, i) \
        using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 1 notitle, \
    mean_f using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 2 title "fine-scale mean", \
    ups_f  using tc:(column(cc) > 0 ? column(cc) : 1/0) \
        with lines ls 3 title "upscaled"

unset multiplot
set output
