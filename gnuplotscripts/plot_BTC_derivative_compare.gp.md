#!/usr/bin/gnuplot

# Set terminal to PNG with good resolution (same as other plots)
set terminal pngcairo enhanced font "Arial,28" size 1200,800

# Data and grid settings
set datafile separator ','
set grid

# Axis labels with larger fonts
set xlabel 't' font "Arial,32"
set ylabel 'c/c_0' font "Arial,32"

# Set logarithmic scale for y-axis
set yrange [1e-4:*]
set ytics (1e-4, 1e-3, 1e-2, 1e-1, 1)
set mytics 10

# Format y-axis labels with superscript notation
set format y "10^{%T}"

# Set y-axis range with minimum of 1e-6
set autoscale fix y
set yrange [1e-4:*]

# Legend position
set key bottom left

# Define line styles
set style line 1 lc rgb '#cccccc' lw 1.5  # Grey for realizations
set style line 2 lc rgb '#000000' lw 3.0  # Black for mean
set style line 3 lc rgb '#d60000' lw 3.0  # Red for upscaled

# Plot for x=0.50
set output 'BTC_derivative_compare_x_0_50.png'
set title 'x=0.50' font "Arial,32"
plot 'BTC_Derivative_Compare.csv' using 1:2 with lines ls 1 notitle, \
     '' using 1:8 with lines ls 1 notitle, \
     '' using 1:14 with lines ls 1 notitle, \
     '' using 1:20 with lines ls 1 notitle, \
     '' using 1:26 with lines ls 1 notitle, \
     '' using 1:32 with lines ls 1 notitle, \
     '' using 1:38 with lines ls 1 notitle, \
     '' using 1:44 with lines ls 1 notitle, \
     '' using 1:50 with lines ls 1 notitle, \
     '' using 1:56 with lines ls 1 notitle, \
     '' using 1:62 with lines ls 1 notitle, \
     '' using 1:68 with lines ls 1 notitle, \
     '' using 1:74 with lines ls 1 notitle, \
     '' using 1:80 with lines ls 1 notitle, \
     '' using 1:86 with lines ls 1 notitle, \
     '' using 1:92 with lines ls 1 notitle, \
     '' using 1:98 with lines ls 1 notitle, \
     '' using 1:104 with lines ls 1 notitle, \
     '' using 1:110 with lines ls 1 notitle, \
     '' using 1:116 with lines ls 1 title 'Realizations', \
     '' using 1:122 with lines ls 2 title 'Mean', \
     '' using 1:128 with lines ls 3 title 'Upscaled'

# Plot for x=1.50
set output 'BTC_derivative_compare_x_1_50.png'
set title 'x=1.50' font "Arial,32"
plot 'BTC_Derivative_Compare.csv' using 1:4 with lines ls 1 notitle, \
     '' using 1:10 with lines ls 1 notitle, \
     '' using 1:16 with lines ls 1 notitle, \
     '' using 1:22 with lines ls 1 notitle, \
     '' using 1:28 with lines ls 1 notitle, \
     '' using 1:34 with lines ls 1 notitle, \
     '' using 1:40 with lines ls 1 notitle, \
     '' using 1:46 with lines ls 1 notitle, \
     '' using 1:52 with lines ls 1 notitle, \
     '' using 1:58 with lines ls 1 notitle, \
     '' using 1:64 with lines ls 1 notitle, \
     '' using 1:70 with lines ls 1 notitle, \
     '' using 1:76 with lines ls 1 notitle, \
     '' using 1:82 with lines ls 1 notitle, \
     '' using 1:88 with lines ls 1 notitle, \
     '' using 1:94 with lines ls 1 notitle, \
     '' using 1:100 with lines ls 1 notitle, \
     '' using 1:106 with lines ls 1 notitle, \
     '' using 1:112 with lines ls 1 notitle, \
     '' using 1:118 with lines ls 1 title 'Realizations', \
     '' using 1:124 with lines ls 2 title 'Mean', \
     '' using 1:130 with lines ls 3 title 'Upscaled'

# Plot for x=2.50
set output 'BTC_derivative_compare_x_2_50.png'
set title 'x=2.50' font "Arial,32"
plot 'BTC_Derivative_Compare.csv' using 1:6 with lines ls 1 notitle, \
     '' using 1:12 with lines ls 1 notitle, \
     '' using 1:18 with lines ls 1 notitle, \
     '' using 1:24 with lines ls 1 notitle, \
     '' using 1:30 with lines ls 1 notitle, \
     '' using 1:36 with lines ls 1 notitle, \
     '' using 1:42 with lines ls 1 notitle, \
     '' using 1:48 with lines ls 1 notitle, \
     '' using 1:54 with lines ls 1 notitle, \
     '' using 1:60 with lines ls 1 notitle, \
     '' using 1:66 with lines ls 1 notitle, \
     '' using 1:72 with lines ls 1 notitle, \
     '' using 1:78 with lines ls 1 notitle, \
     '' using 1:84 with lines ls 1 notitle, \
     '' using 1:90 with lines ls 1 notitle, \
     '' using 1:96 with lines ls 1 notitle, \
     '' using 1:102 with lines ls 1 notitle, \
     '' using 1:108 with lines ls 1 notitle, \
     '' using 1:114 with lines ls 1 notitle, \
     '' using 1:120 with lines ls 1 title 'Realizations', \
     '' using 1:126 with lines ls 2 title 'Mean', \
     '' using 1:132 with lines ls 3 title 'Upscaled'

unset output
