#!/usr/bin/gnuplot

set datafile separator ","
set logscale x
set grid

set style line 1 lc '#cccccc' lw 1.0          # grey: realizations
set style line 2 lc '#000000' lw 3.0          # black thick: mean
set style line 3 lc '#2E73B5' lw 2.5 dt 2    # blue dashed: fit

# ---------------------------------------------------------------------------
# Diffusion X
# ---------------------------------------------------------------------------
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set output 'diffusion_x_correlation.png'

set xrange [0.01 : 0.3]
set yrange [0 : 1.0005]

set xlabel '{/Symbol D}x' font "Arial,32"
set ylabel 'E[{/Symbol w}(x){/Symbol w}(x+{/Symbol D}x)]' font "Arial,32"
set key bottom left font "Arial,24"

plot 'diffusion_x_correlations.txt' every ::1 u 1:2  w l ls 1 notitle, \
     ''                                                    every ::1 u 3:4  w l ls 1 notitle, \
     ''                                                    every ::1 u 5:6  w l ls 1 notitle, \
     ''                                                    every ::1 u 7:8  w l ls 1 notitle, \
     ''                                                    every ::1 u 9:10 w l ls 1 notitle, \
     ''                                                    every ::1 u 11:12 w l ls 1 notitle, \
     ''                                                    every ::1 u 13:14 w l ls 1 notitle, \
     ''                                                    every ::1 u 15:16 w l ls 1 notitle, \
     ''                                                    every ::1 u 17:18 w l ls 1 notitle, \
     ''                                                    every ::1 u 19:20 w l ls 1 notitle, \
     ''                                                    every ::1 u 21:22 w l ls 1 notitle, \
     ''                                                    every ::1 u 23:24 w l ls 1 notitle, \
     ''                                                    every ::1 u 25:26 w l ls 1 notitle, \
     ''                                                    every ::1 u 27:28 w l ls 1 notitle, \
     ''                                                    every ::1 u 29:30 w l ls 1 notitle, \
     ''                                                    every ::1 u 31:32 w l ls 1 notitle, \
     ''                                                    every ::1 u 33:34 w l ls 1 notitle, \
     ''                                                    every ::1 u 35:36 w l ls 1 notitle, \
     ''                                                    every ::1 u 37:38 w l ls 1 notitle, \
     ''                                                    every ::1 u 39:40 w l ls 1 title 'Realizations', \
     'diffusion_x_correlations_mean.txt'        u 1:2 w lp ls 2 pt 6 ps 0.45 title 'Mean', \
     'diffusion_x_correlations_mean_fitted.txt' u 1:2 w l  ls 3          title 'Exponential fit'

# ---------------------------------------------------------------------------
# Diffusion Y
# ---------------------------------------------------------------------------
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set output 'diffusion_y_correlation.png'

set xrange [0.01 : 0.3]
set yrange [0 : 1.0005]

set xlabel '{/Symbol D}y' font "Arial,32"
set ylabel 'E[{/Symbol w}(y){/Symbol w}(y+{/Symbol D}y)]' font "Arial,32"
set key bottom left font "Arial,24"

plot 'diffusion_y_correlations.txt' every ::1 u 1:2  w l ls 1 notitle, \
     ''                                                    every ::1 u 3:4  w l ls 1 notitle, \
     ''                                                    every ::1 u 5:6  w l ls 1 notitle, \
     ''                                                    every ::1 u 7:8  w l ls 1 notitle, \
     ''                                                    every ::1 u 9:10 w l ls 1 notitle, \
     ''                                                    every ::1 u 11:12 w l ls 1 notitle, \
     ''                                                    every ::1 u 13:14 w l ls 1 notitle, \
     ''                                                    every ::1 u 15:16 w l ls 1 notitle, \
     ''                                                    every ::1 u 17:18 w l ls 1 notitle, \
     ''                                                    every ::1 u 19:20 w l ls 1 notitle, \
     ''                                                    every ::1 u 21:22 w l ls 1 notitle, \
     ''                                                    every ::1 u 23:24 w l ls 1 notitle, \
     ''                                                    every ::1 u 25:26 w l ls 1 notitle, \
     ''                                                    every ::1 u 27:28 w l ls 1 notitle, \
     ''                                                    every ::1 u 29:30 w l ls 1 notitle, \
     ''                                                    every ::1 u 31:32 w l ls 1 notitle, \
     ''                                                    every ::1 u 33:34 w l ls 1 notitle, \
     ''                                                    every ::1 u 35:36 w l ls 1 notitle, \
     ''                                                    every ::1 u 37:38 w l ls 1 notitle, \
     ''                                                    every ::1 u 39:40 w l ls 1 title 'Realizations', \
     'diffusion_y_correlations_mean.txt'        u 1:2 w lp ls 2 pt 6 ps 0.45 title 'Mean', \
     'diffusion_y_correlations_mean_fitted.txt' u 1:2 w l  ls 3          title 'Exponential fit'
