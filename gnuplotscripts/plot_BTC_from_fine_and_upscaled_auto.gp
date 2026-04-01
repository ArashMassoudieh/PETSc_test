# File overview: automatic BTC plotting directly from fine_r* folders and
# upscaled_mean folder, without requiring prebuilt x=*BTC_Compare.csv files.
#
# plot_BTC_from_fine_and_upscaled_auto.gp
#
# Gnuplot 6-safe version:
# - top-level simple one-line if statements are kept
# - inside bracketed blocks (while {...}) use only block-style if {...}
# - no ternary operator
# - no arrays
# - no do for / plot for bracket loops

reset
set datafile separator ','
set datafile missing ''
set datafile columnheaders

if (!exists("out_prefix")) out_prefix = "BTC_compare_from_folders"
if (!exists("mode"))       mode = "cdf"
if (!exists("x_main_max")) x_main_max = 10
if (!exists("x_log_max"))  x_log_max  = 20
if (!exists("y_log_min"))  y_log_min  = 1e-6
if (!exists("fine_glob"))  fine_glob  = "fine_r*"
if (!exists("show_title")) show_title = 1

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'

root_dir = GPVAL_PWD

mode_lc = mode
mode_ok = 0
if (mode_lc eq "cdf") mode_ok = 1
if (mode_lc eq "pdf") mode_ok = 1
if (mode_ok == 0) mode_lc = "cdf"

fine_pattern = "r*_BTC_FineScaled.csv"
up_file      = "upscaled_mean/upscaled_BTC_Upscaled.csv"
mode_tag     = "cdf"
mode_label   = "BTC"
y_label      = "c/c_0"

if (mode_lc eq "pdf") fine_pattern = "r*_BTC_FineScaled_derivative.csv"
if (mode_lc eq "pdf") up_file      = "upscaled_mean/upscaled_BTC_Upscaled_derivative.csv"
if (mode_lc eq "pdf") mode_tag     = "pdf"
if (mode_lc eq "pdf") mode_label   = "BTC PDF"
if (mode_lc eq "pdf") y_label      = "dc/c_0dt"

dir_list = system(sprintf("find . -maxdepth 1 -type d -name '%s' | sort", fine_glob))
if (strlen(dir_list) == 0) exit
if (system(sprintf("test -f '%s'", up_file)) ne "") exit

rep_file = ""
fine_files = ""
n_dirs = words(dir_list)
valid_dirs = 0

di = 1
while (di <= n_dirs) {
    d = word(dir_list, di)
    cand = system(sprintf("ls '%s'/%s 2>/dev/null | head -n 1", d, fine_pattern))
    ff = word(cand, 1)

    ff_ok = 0
    if (strlen(ff) > 0) {
        ff_ok = 1
    }

    rep_empty = 0
    if (strlen(rep_file) == 0) {
        rep_empty = 1
    }

    fine_empty = 0
    if (strlen(fine_files) == 0) {
        fine_empty = 1
    }

    fine_same = 0
    if (fine_files eq ff) {
        fine_same = 1
    }

    if (ff_ok == 1) {
        valid_dirs = valid_dirs + 1
    }

    set_rep = 0
    if (ff_ok == 1) {
        set_rep = set_rep + 1
    }
    if (rep_empty == 1) {
        set_rep = set_rep + 1
    }
    if (set_rep == 2) {
        rep_file = ff
    }

    set_first_fine = 0
    if (ff_ok == 1) {
        set_first_fine = set_first_fine + 1
    }
    if (fine_empty == 1) {
        set_first_fine = set_first_fine + 1
    }
    if (set_first_fine == 2) {
        fine_files = ff
    }

    append_fine = 0
    if (ff_ok == 1) {
        append_fine = append_fine + 1
    }
    if (fine_empty == 0) {
        append_fine = append_fine + 1
    }
    if (fine_same == 0) {
        append_fine = append_fine + 1
    }
    if (append_fine == 3) {
        fine_files = fine_files . " " . ff
    }

    di = di + 1
}

if (strlen(rep_file) == 0) exit

stats rep_file nooutput
ncol_f = STATS_columns
npairs_f = int(ncol_f/2)
if (npairs_f < 1) exit

stats up_file nooutput
ncol_u = STATS_columns
npairs_u = int(ncol_u/2)
if (npairs_u < 1) exit

npairs = npairs_f
if (npairs_u < npairs_f) npairs = npairs_u

print sprintf("Found %d fine folders; usable fine files in %d folder(s); plotting %d BTC pair(s) in mode=%s", n_dirs, valid_dirs, npairs, mode_tag)

mean_file = ""
mean_candidates = "BTC_mean_cdf.csv BTC_mean_transport_cdf.csv BTC_mean.csv"
if (mode_lc eq "pdf") mean_candidates = "BTC_mean_pdf.csv BTC_mean_transport_full.csv BTC_mean.csv"

nmc = words(mean_candidates)
mi = 1
while (mi <= nmc) {
    mf = word(mean_candidates, mi)

    found_mean = 0
    if (system(sprintf("test -f '%s'", mf)) eq "") {
        found_mean = 1
    }

    if (found_mean == 1) {
        mean_file = mf
    }
    if (found_mean == 1) {
        mi = nmc + 1
    }
    if (found_mean == 0) {
        mi = mi + 1
    }
}

if (strlen(mean_file) > 0) print sprintf("Using mean file: %s", mean_file)
if (strlen(mean_file) == 0) print "No root mean file found; plots will show fine realizations + upscaled only."

n_valid = words(fine_files)

p = 1
while (p <= npairs) {
    ft_col = 2*p - 1
    fc_col = 2*p
    ut_col = ft_col
    uc_col = fc_col

    out_png = sprintf("%s_pair%02d_%s_auto.png", out_prefix, p, mode_tag)
    set output out_png
    set multiplot layout 1,1

    set grid
    set key top right
    set xlabel 'Time' font 'Arial,32'
    set ylabel y_label font 'Arial,32'
    if (show_title) {
        set label 1 sprintf('%s from fine folders + upscaled_mean (pair %d)', mode_label, p) at graph 0.02,0.96 front tc rgb '#333333' font 'Arial,18'
    }
    if (!show_title) {
        unset label 1
    }
    set format y "%.1f"
    set yrange [0:*]
    set xrange [0:x_main_max]
    unset logscale y
    unset y2tics
    set origin 0,0
    set size 1,1

    plotcmd = "plot "
    di = 1
    while (di <= n_valid) {
        ff = word(fine_files, di)
        if (di > 1) {
            plotcmd = plotcmd . ", "
        }
        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 1 lc rgb '#CCCCCC' notitle", ff, ft_col, fc_col)
        di = di + 1
    }

    if (strlen(mean_file) > 0) {
        if (n_valid > 0) {
            plotcmd = plotcmd . ", "
        }
        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 3 lc rgb '#000000' title 'Mean'", mean_file, ft_col, fc_col)
        plotcmd = plotcmd . sprintf(", '%s' using %d:%d with lines lw 3 lc rgb '#FF0000' title 'Upscaled'", up_file, ut_col, uc_col)
    } else {
        if (n_valid > 0) {
            plotcmd = plotcmd . ", "
        }
        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 3 lc rgb '#FF0000' title 'Upscaled'", up_file, ut_col, uc_col)
    }

    eval plotcmd

    set origin 0.40, 0.15
    set size 0.55, 0.65
    set key off
    unset xlabel
    unset ylabel
    set tics font 'Arial,16'
    set logscale y
    set format y "10^{%T}"
    set yrange [y_log_min:*]
    set xrange [0:x_log_max]

    plotcmd = "plot "
    di = 1
    while (di <= n_valid) {
        ff = word(fine_files, di)
        if (di > 1) {
            plotcmd = plotcmd . ", "
        }
        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 1 lc rgb '#CCCCCC' notitle", ff, ft_col, fc_col)
        di = di + 1
    }

    if (strlen(mean_file) > 0) {
        if (n_valid > 0) {
            plotcmd = plotcmd . ", "
        }
        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 3 lc rgb '#000000' notitle", mean_file, ft_col, fc_col)
        plotcmd = plotcmd . sprintf(", '%s' using %d:%d with lines lw 3 lc rgb '#FF0000' notitle", up_file, ut_col, uc_col)
    }
    if (strlen(mean_file) == 0) {
        if (n_valid > 0) {
            plotcmd = plotcmd . ", "
        }
        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 3 lc rgb '#FF0000' notitle", up_file, ut_col, uc_col)
    }

    eval plotcmd

    unset multiplot
    unset label 1
    print sprintf("Wrote %s", out_png)

    p = p + 1
}

print "Done."
