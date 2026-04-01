# plot_BTC_from_fine_and_upscaled_auto.gp
# Auto plot BTC CDF + PDF from fine_r* folders and upscaled_mean
# Style matched to old plot_BTC_both_100.gp
#
# Gnuplot 6-safe:
# - no mode needed; does both cdf and pdf
# - no ternary operator
# - no do for / plot for
# - no old-style if/else inside bracket blocks
# - avoids fragile negated string expressions

reset
set datafile separator ','
set datafile missing ''
set datafile columnheaders

if (!exists("out_prefix")) out_prefix = "BTC_compare"
if (!exists("x_main_max")) x_main_max = 10
if (!exists("x_log_max"))  x_log_max  = 20
if (!exists("y_log_min"))  y_log_min  = 1e-6
if (!exists("fine_glob"))  fine_glob  = "fine_r*"

set terminal pngcairo size 1200,800 enhanced font 'Arial,28'

root_dir = GPVAL_PWD

# ------------------------------------------------------------
# discover fine folders once
# ------------------------------------------------------------
dir_list = system(sprintf("find . -maxdepth 1 -type d -name '%s' | sort", fine_glob))
if (strlen(dir_list) == 0) {
    print "No fine_r* folders found."
    exit
}

# ------------------------------------------------------------
# kind = 1 -> cdf
# kind = 2 -> pdf
# ------------------------------------------------------------
kind = 1
while (kind <= 2) {

    fine_pattern    = "r*_BTC_FineScaled.csv"
    up_file         = "upscaled_mean/upscaled_BTC_Upscaled.csv"
    mean_candidates = "BTC_mean_cdf.csv BTC_mean_transport_cdf.csv BTC_mean.csv"
    mode_tag        = "cdf"
    y_label         = "c/c_0"

    if (kind == 2) {
        fine_pattern    = "r*_BTC_FineScaled_derivative.csv"
        up_file         = "upscaled_mean/upscaled_BTC_Upscaled_derivative.csv"
        mean_candidates = "BTC_mean_pdf.csv BTC_mean_transport_full.csv BTC_mean.csv"
        mode_tag        = "pdf"
        y_label         = "dc/c_0dt"
    }

    have_up = 0
    if (system(sprintf("test -f '%s'", up_file)) eq "") {
        have_up = 1
    }

    if (have_up == 1) {

        rep_file   = ""
        fine_files = ""
        n_dirs     = words(dir_list)
        valid_dirs = 0

        di = 1
        while (di <= n_dirs) {
            d    = word(dir_list, di)
            cand = system(sprintf("ls '%s'/%s 2>/dev/null | head -n 1", d, fine_pattern))
            ff   = word(cand, 1)

            ff_ok = 0
            if (strlen(ff) > 0) {
                ff_ok = 1
            }

            if (ff_ok == 1) {
                valid_dirs = valid_dirs + 1
            }

            if ((ff_ok == 1) && (strlen(rep_file) == 0)) {
                rep_file = ff
            }

            if ((ff_ok == 1) && (strlen(fine_files) == 0)) {
                fine_files = ff
            }

            if ((ff_ok == 1) && (strlen(fine_files) > 0)) {
                same_file = 0
                if (fine_files eq ff) {
                    same_file = 1
                }
                if (same_file == 0) {
                    fine_files = fine_files . " " . ff
                }
            }

            di = di + 1
        }

        if (strlen(rep_file) > 0) {

            stats rep_file nooutput
            ncol_f   = STATS_columns
            npairs_f = int(ncol_f/2)

            stats up_file nooutput
            ncol_u   = STATS_columns
            npairs_u = int(ncol_u/2)

            if ((npairs_f >= 1) && (npairs_u >= 1)) {

                npairs = npairs_f
                if (npairs_u < npairs_f) {
                    npairs = npairs_u
                }

                print sprintf("Found %d fine folders; usable fine files in %d folder(s); plotting %d BTC pair(s) for %s", \
                              n_dirs, valid_dirs, npairs, mode_tag)

                # ------------------------------------------------------------
                # mean file detection
                # ------------------------------------------------------------
                mean_file = ""
                nmc = words(mean_candidates)
                mi = 1
                while (mi <= nmc) {
                    mf = word(mean_candidates, mi)
                    found_mean = 0
                    if (system(sprintf("test -f '%s'", mf)) eq "") {
                        found_mean = 1
                    }
                    if ((found_mean == 1) && (strlen(mean_file) == 0)) {
                        mean_file = mf
                    }
                    mi = mi + 1
                }

                if (strlen(mean_file) > 0) {
                    print sprintf("Using mean file for %s: %s", mode_tag, mean_file)
                }
                if (strlen(mean_file) == 0) {
                    print sprintf("No mean file found for %s; plots will show realizations + upscaled only.", mode_tag)
                }

                n_valid = words(fine_files)

                # ------------------------------------------------------------
                # one plot per pair
                # ------------------------------------------------------------
                p = 1
                while (p <= npairs) {
                    ft_col = 2*p - 1
                    fc_col = 2*p
                    ut_col = ft_col
                    uc_col = fc_col

                    # old-style x naming
                    x_guess = sprintf("x%.2f", 0.5 + (p-1)*1.0)
                    if (p == 1) {
                        x_guess = "x0.50"
                    }
                    if (p == 2) {
                        x_guess = "x1.50"
                    }
                    if (p == 3) {
                        x_guess = "x2.50"
                    }

                    out_png = sprintf("%s_%s_%s_combined.png", out_prefix, x_guess, mode_tag)

                    set output out_png
                    set multiplot layout 1,1

                    # ========================================================
                    # Main plot
                    # ========================================================
                    set grid
                    if (p == 1) {
                        set key top right
                    }
                    if (p != 1) {
                        set key off
                    }

                    set xlabel 'Time' font 'Arial,32'
                    set ylabel y_label font 'Arial,32'
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

                        last_and_first = 0
                        if (di == n_valid) {
                            if (p == 1) {
                                last_and_first = 1
                            }
                        }

                        if (last_and_first == 1) {
                            plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 1 lc rgb '#CCCCCC' title 'Realizations'", \
                                                        ff, ft_col, fc_col)
                        }
                        if (last_and_first == 0) {
                            plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 1 lc rgb '#CCCCCC' notitle", \
                                                        ff, ft_col, fc_col)
                        }

                        di = di + 1
                    }

                    if (strlen(mean_file) > 0) {
                        if (n_valid > 0) {
                            plotcmd = plotcmd . ", "
                        }
                        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 3 lc rgb '#000000' title 'Mean'", \
                                                    mean_file, ft_col, fc_col)
                        plotcmd = plotcmd . sprintf(", '%s' using %d:%d with lines lw 3 lc rgb '#FF0000' title 'Upscaled'", \
                                                    up_file, ut_col, uc_col)
                    }

                    if (strlen(mean_file) == 0) {
                        if (n_valid > 0) {
                            plotcmd = plotcmd . ", "
                        }
                        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 3 lc rgb '#FF0000' title 'Upscaled'", \
                                                    up_file, ut_col, uc_col)
                    }

                    eval plotcmd

                    # ========================================================
                    # Inset plot
                    # ========================================================
                    if (p == 1) {
                        set origin 0.40, 0.15
                    }
                    if (p != 1) {
                        set origin 0.40, 0.30
                    }
                    set size 0.55, 0.65
                    set key off
                    unset title
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
                        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 1 lc rgb '#CCCCCC' notitle", \
                                                    ff, ft_col, fc_col)
                        di = di + 1
                    }

                    if (strlen(mean_file) > 0) {
                        if (n_valid > 0) {
                            plotcmd = plotcmd . ", "
                        }
                        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 2 lc rgb '#000000' notitle", \
                                                    mean_file, ft_col, fc_col)
                        plotcmd = plotcmd . sprintf(", '%s' using %d:%d with lines lw 2 lc rgb '#FF0000' notitle", \
                                                    up_file, ut_col, uc_col)
                    }

                    if (strlen(mean_file) == 0) {
                        if (n_valid > 0) {
                            plotcmd = plotcmd . ", "
                        }
                        plotcmd = plotcmd . sprintf("'%s' using %d:%d with lines lw 2 lc rgb '#FF0000' notitle", \
                                                    up_file, ut_col, uc_col)
                    }

                    eval plotcmd

                    unset multiplot
                    print sprintf("Wrote %s", out_png)

                    p = p + 1
                }

            }

            if ((npairs_f < 1) || (npairs_u < 1)) {
                print sprintf("Skipping %s: invalid fine/upscaled column pairs.", mode_tag)
            }

        }

        if (strlen(rep_file) == 0) {
            print sprintf("No usable fine files found for %s.", mode_tag)
        }

    }

    if (have_up == 0) {
        print sprintf("Missing upscaled file for %s: %s", mode_tag, up_file)
    }

    kind = kind + 1
}

print "Done."
