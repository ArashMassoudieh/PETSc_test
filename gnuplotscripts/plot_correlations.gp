set datafile separator ","

# ---------------------------------------------------------------------------
# Diffusion Y
# ---------------------------------------------------------------------------
set terminal png size 900,550 font "DejaVu Sans,13"
set output "diffusion_y_correlation.png"

set logscale x
set xrange [8e-4 : 1.2e-1]
set yrange [0.3 : 1.05]

set xlabel "r" offset 0,-0.5
set ylabel "{ρ}(r)" offset -0.8,0

set key top left box linetype 0 samplen 2.5
set key width 1

set grid xtics ytics lt 1 lw 0.5 lc "#AAAAAA"

plot "diffusion_y_correlations_mean.txt" \
        using 1:2 title "Mean (data)" \
        with linespoints lw 2.5 pt 6 ps 0.7 lc "#000000",\
     "diffusion_y_correlations_mean_fitted.txt" \
        using 1:2 title "Matérn fit" \
        with lines lw 2.5 lc "#2E73B5" dt 2

# ---------------------------------------------------------------------------
# Diffusion X
# ---------------------------------------------------------------------------
set terminal png size 900,550 font "DejaVu Sans,13"
set output "diffusion_x_correlation.png"

set logscale x
set xrange [8e-4 : 1.2e-1]
set yrange [0.85 : 1.005]

set xlabel "r" offset 0,-0.5
set ylabel "{ρ}(r)" offset -0.8,0

set key top left box linetype 0 samplen 2.5
set key width 1

set grid xtics ytics lt 1 lw 0.5 lc "#AAAAAA"

plot "diffusion_x_correlations_mean.txt" \
        using 1:2 title "Mean (data)" \
        with linespoints lw 2.5 pt 6 ps 0.7 lc "#000000",\
     "diffusion_x_correlations_mean_fitted.txt" \
        using 1:2 title "Matérn fit" \
        with lines lw 2.5 lc "#2E73B5" dt 2
