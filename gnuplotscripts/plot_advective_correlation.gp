#!/usr/bin/gnuplot

# Set terminal to PNG with good resolution
set terminal pngcairo enhanced font "Arial,28" size 1200,800

# Output file
set output 'advective_correlation.png'

# Set plot style
set style line 1 lc rgb '#cccccc' lw 1.5  # Grey for realizations
set style line 2 lc rgb '#000000' lw 3.0  # Black thick for fit

# Axis labels
set xlabel '{/Symbol D}x' font "Arial,32"
set ylabel 'E[{/Symbol w}(x){/Symbol w}(x+{/Symbol D}x)]' font "Arial,32"

# Set logarithmic scale for x-axis
set logscale x

# Grid
set grid

# Key/legend
set key bottom left

# Set datafile separator to comma
set datafile separator ","

# Define the exponential fit function with lambda_a = 0.294484
fit_func(x) = exp(-x/0.294484)

# Plot all 20 realizations in grey (skip header, columns 1-40 contain pairs)
# Using column numbers: 1,2 for first realization, 3,4 for second, etc.
plot 'advective_correlations.txt' every ::1 u 1:2 w l ls 1 notitle, \
     '' every ::1 u 3:4 w l ls 1 notitle, \
     '' every ::1 u 5:6 w l ls 1 notitle, \
     '' every ::1 u 7:8 w l ls 1 notitle, \
     '' every ::1 u 9:10 w l ls 1 notitle, \
     '' every ::1 u 11:12 w l ls 1 notitle, \
     '' every ::1 u 13:14 w l ls 1 notitle, \
     '' every ::1 u 15:16 w l ls 1 notitle, \
     '' every ::1 u 17:18 w l ls 1 notitle, \
     '' every ::1 u 19:20 w l ls 1 notitle, \
     '' every ::1 u 21:22 w l ls 1 notitle, \
     '' every ::1 u 23:24 w l ls 1 notitle, \
     '' every ::1 u 25:26 w l ls 1 notitle, \
     '' every ::1 u 27:28 w l ls 1 notitle, \
     '' every ::1 u 29:30 w l ls 1 notitle, \
     '' every ::1 u 31:32 w l ls 1 notitle, \
     '' every ::1 u 33:34 w l ls 1 notitle, \
     '' every ::1 u 35:36 w l ls 1 notitle, \
     '' every ::1 u 37:38 w l ls 1 notitle, \
     '' every ::1 u 39:40 w l ls 1 title 'Realizations', \
     fit_func(x) w l ls 2 lw 4 title 'exp(-{/Symbol D}x/{/Symbol l}_a)'
