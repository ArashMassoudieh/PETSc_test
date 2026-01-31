#!/usr/bin/gnuplot

# Set terminal to PNG with good resolution
set terminal pngcairo enhanced font "Arial,28" size 1200,800

# Output file
set output 'velocity_correlation.png'

# Set plot style
set style line 1 lc rgb '#cccccc' lw 1.5  # Grey for realizations
set style line 2 lc rgb '#000000' lw 3.0  # Black thick for fit

# Axis labels
set xlabel '{/Symbol D}y' font "Arial,32"
set ylabel 'E[{/Symbol w}(y){/Symbol w}(y+{/Symbol D}y)]' font "Arial,32"
set yrange [0.5:1]
# Set logarithmic scale for y-axis
set logscale y

# Grid
set grid

# Key/legend
set key bottom left

# Set datafile separator to comma
set datafile separator ","

# Define the exponential fit function
fit_func(x) = exp(-x/0.354109)

# Plot all 20 realizations in grey and the fit function in black
plot 'fine_r0001/r0001_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0002/r0002_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0003/r0003_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0004/r0004_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0005/r0005_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0006/r0006_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0007/r0007_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0008/r0008_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0009/r0009_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0010/r0010_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0011/r0011_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0012/r0012_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0013/r0013_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0014/r0014_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0015/r0015_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0016/r0016_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0017/r0017_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0018/r0018_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0019/r0019_velocity_correlation_y.txt' u 1:2 w l ls 1 notitle, \
     'fine_r0020/r0020_velocity_correlation_y.txt' u 1:2 w l ls 1 title 'Realizations', \
     fit_func(x) w l ls 2 lw 4 title 'exp(-{/Symbol D}y/{/Symbol l}_y)'
