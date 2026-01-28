#!/usr/bin/gnuplot

# Set terminal and output
set terminal pngcairo enhanced size 1200,800 font "Arial,12"
set output 'x_2_50BTC_Compare.png'

# Optionally, you can use other terminals like:
# set terminal qt size 1200,800  # For interactive window
# set terminal pdf size 12,8     # For PDF output

# Plot styling
set grid
set key top left

# Axis labels
set xlabel 't'
set ylabel 'Value'

# Use logarithmic scale for y-axis since values span many orders of magnitude
set logscale y

# Title
set title 'Comparison of x=2.50 Data Series'

# Data file settings
set datafile separator ","

# Plot all data series
# Using 'using' to specify column pairs and 'every ::1' to skip the header
plot '/mnt/user-data/uploads/x_2_50BTC_Compare.csv' every ::1 using 1:2 with lines title 'r0001\_x=2.50' lw 2, \
     '' every ::1 using 3:4 with lines title 'r0002\_x=2.50' lw 2, \
     '' every ::1 using 5:6 with lines title 'r0003\_x=2.50' lw 2, \
     '' every ::1 using 7:8 with lines title 'r0004\_x=2.50' lw 2, \
     '' every ::1 using 9:10 with lines title 'Upscaled\_x=2.50' lw 2

# Alternative: If you prefer linear scale instead of log scale, comment out the logscale line above
# and uncomment the following:
# unset logscale y
# set format y "%.2e"

print "Plot saved to x_2_50BTC_Compare.png"
