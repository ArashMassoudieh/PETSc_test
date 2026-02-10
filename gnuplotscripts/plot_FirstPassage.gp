set terminal pngcairo enhanced size 800,500 font "Arial,14"
set output "FirstPassageTimeDistribution.png"

set xlabel "First Passage Time" font "Arial,16"
set ylabel "Probability Density" font "Arial,16"
set title "First Passage Time Distribution (Ellipse Escape)"

set datafile separator ","
set style fill solid 0.4 border -1
set boxwidth 0.2 relative

set xrange [0:*]
set yrange [0:*]
set grid lc rgb "#cccccc"
set key off

plot "FirstPassageTimeDistribution.csv" using 1:2 with boxes lc rgb "#4472C4" notitle, \
     "" using 1:2 with lines lw 2 lc rgb "#C00000" notitle
