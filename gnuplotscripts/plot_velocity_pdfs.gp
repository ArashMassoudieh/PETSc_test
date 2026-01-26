#!/usr/bin/gnuplot

# Set terminal to PNG with good resolution
set terminal pngcairo enhanced font "Arial,28" size 1200,800

# Output file
set output 'velocity_pdfs.png'

# Set plot style
set style line 1 lc rgb '#cccccc' lw 1.5  # Grey for realizations
set style line 2 lc rgb '#000000' lw 3.0  # Black thick for mean

# Axis labels and titles
set xlabel 'Velocity (v_x)' font "Arial,32"
set ylabel 'PDF' font "Arial,32"

# Set axis ranges
set xrange [0.01:15]

# Set logarithmic scales for both axes
set logscale x
set logscale y

# Format y-axis labels with superscript notation
set format y "10^{%T}"

# Grid
set grid

# Key/legend
set key bottom left

# Set datafile separator to comma
set datafile separator ","

# Plot all 20 realizations in grey (skip first line which is header)
# Using column numbers: 1,2 for first realization, 3,4 for second, etc.
plot 'qx_pdfs.txt' every ::1 u 1:2 w l ls 1 notitle, \
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
     'qx_mean_pdf.txt' u 1:2 w l ls 2 lw 4 title 'Mean PDF'
