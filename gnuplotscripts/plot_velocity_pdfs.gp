reset

# Terminal settings
set terminal pngcairo size 1200,800 enhanced font 'Arial,28'
set output 'velocity_pdfs.png'

# Log-log axes
set logscale x
set logscale y
set format y "10^{%T}"
set xlabel 'Velocity (v_x)' font 'Arial,32'
set ylabel 'PDF' font 'Arial,32'

# Axis ranges (matching your plot)
set xrange [0.001:20]
set yrange [1e-3:10]

# Grid and legend
set grid
set key bottom left

# CSV format
set datafile separator ','

# Plot: 20 realizations in gray + mean in black
plot 'qx_pdfs.txt' using 1:2 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 3:4 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 5:6 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 7:8 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 9:10 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 11:12 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 13:14 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 15:16 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 17:18 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 19:20 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 21:22 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 23:24 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 25:26 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 27:28 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 29:30 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 31:32 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 33:34 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 35:36 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 37:38 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" title 'Realizations', \
     'qx_mean_pdf.txt' using 1:2 with lines lw 3 lc rgb "#000000" title 'Mean PDF'
