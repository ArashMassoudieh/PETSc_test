#!/usr/bin/gnuplot

set datafile separator ","
unset logscale x
set logscale y
set grid

set style line 1 lc '#cccccc' lw 1.0          # grey: realizations
set style line 2 lc '#000000' lw 3.0          # black thick: mean
set style line 3 lc '#2E73B5' lw 2.5 dt 2    # blue dashed: fit

# ---------------------------------------------------------------------------
# Diffusion X
# ---------------------------------------------------------------------------
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set output 'diffusion_x_correlation.png'

set xrange [0 : 1.0]
set yrange [0.1 : 1.0]

set xlabel '{/Symbol D}x' font "Arial,32"
set ylabel 'E[{/Symbol w}(x){/Symbol w}(x+{/Symbol D}x)]' font "Arial,32"
set key bottom left font "Arial,24"
set format y "10^{%T}"

plot 'diffusion_x_correlations.txt' every ::1 u 1:2 w l ls 1 notitle, \
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
     '' every ::1 u 39:40 w l ls 1 notitle, \
     '' every ::1 u 41:42 w l ls 1 notitle, \
     '' every ::1 u 43:44 w l ls 1 notitle, \
     '' every ::1 u 45:46 w l ls 1 notitle, \
     '' every ::1 u 47:48 w l ls 1 notitle, \
     '' every ::1 u 49:50 w l ls 1 notitle, \
     '' every ::1 u 51:52 w l ls 1 notitle, \
     '' every ::1 u 53:54 w l ls 1 notitle, \
     '' every ::1 u 55:56 w l ls 1 notitle, \
     '' every ::1 u 57:58 w l ls 1 notitle, \
     '' every ::1 u 59:60 w l ls 1 notitle, \
     '' every ::1 u 61:62 w l ls 1 notitle, \
     '' every ::1 u 63:64 w l ls 1 notitle, \
     '' every ::1 u 65:66 w l ls 1 notitle, \
     '' every ::1 u 67:68 w l ls 1 notitle, \
     '' every ::1 u 69:70 w l ls 1 notitle, \
     '' every ::1 u 71:72 w l ls 1 notitle, \
     '' every ::1 u 73:74 w l ls 1 notitle, \
     '' every ::1 u 75:76 w l ls 1 notitle, \
     '' every ::1 u 77:78 w l ls 1 notitle, \
     '' every ::1 u 79:80 w l ls 1 notitle, \
     '' every ::1 u 81:82 w l ls 1 notitle, \
     '' every ::1 u 83:84 w l ls 1 notitle, \
     '' every ::1 u 85:86 w l ls 1 notitle, \
     '' every ::1 u 87:88 w l ls 1 notitle, \
     '' every ::1 u 89:90 w l ls 1 notitle, \
     '' every ::1 u 91:92 w l ls 1 notitle, \
     '' every ::1 u 93:94 w l ls 1 notitle, \
     '' every ::1 u 95:96 w l ls 1 notitle, \
     '' every ::1 u 97:98 w l ls 1 notitle, \
     '' every ::1 u 99:100 w l ls 1 notitle, \
     '' every ::1 u 101:102 w l ls 1 notitle, \
     '' every ::1 u 103:104 w l ls 1 notitle, \
     '' every ::1 u 105:106 w l ls 1 notitle, \
     '' every ::1 u 107:108 w l ls 1 notitle, \
     '' every ::1 u 109:110 w l ls 1 notitle, \
     '' every ::1 u 111:112 w l ls 1 notitle, \
     '' every ::1 u 113:114 w l ls 1 notitle, \
     '' every ::1 u 115:116 w l ls 1 notitle, \
     '' every ::1 u 117:118 w l ls 1 notitle, \
     '' every ::1 u 119:120 w l ls 1 notitle, \
     '' every ::1 u 121:122 w l ls 1 notitle, \
     '' every ::1 u 123:124 w l ls 1 notitle, \
     '' every ::1 u 125:126 w l ls 1 notitle, \
     '' every ::1 u 127:128 w l ls 1 notitle, \
     '' every ::1 u 129:130 w l ls 1 notitle, \
     '' every ::1 u 131:132 w l ls 1 notitle, \
     '' every ::1 u 133:134 w l ls 1 notitle, \
     '' every ::1 u 135:136 w l ls 1 notitle, \
     '' every ::1 u 137:138 w l ls 1 notitle, \
     '' every ::1 u 139:140 w l ls 1 notitle, \
     '' every ::1 u 141:142 w l ls 1 notitle, \
     '' every ::1 u 143:144 w l ls 1 notitle, \
     '' every ::1 u 145:146 w l ls 1 notitle, \
     '' every ::1 u 147:148 w l ls 1 notitle, \
     '' every ::1 u 149:150 w l ls 1 notitle, \
     '' every ::1 u 151:152 w l ls 1 notitle, \
     '' every ::1 u 153:154 w l ls 1 notitle, \
     '' every ::1 u 155:156 w l ls 1 notitle, \
     '' every ::1 u 157:158 w l ls 1 notitle, \
     '' every ::1 u 159:160 w l ls 1 notitle, \
     '' every ::1 u 161:162 w l ls 1 notitle, \
     '' every ::1 u 163:164 w l ls 1 notitle, \
     '' every ::1 u 165:166 w l ls 1 notitle, \
     '' every ::1 u 167:168 w l ls 1 notitle, \
     '' every ::1 u 169:170 w l ls 1 notitle, \
     '' every ::1 u 171:172 w l ls 1 notitle, \
     '' every ::1 u 173:174 w l ls 1 notitle, \
     '' every ::1 u 175:176 w l ls 1 notitle, \
     '' every ::1 u 177:178 w l ls 1 notitle, \
     '' every ::1 u 179:180 w l ls 1 notitle, \
     '' every ::1 u 181:182 w l ls 1 notitle, \
     '' every ::1 u 183:184 w l ls 1 notitle, \
     '' every ::1 u 185:186 w l ls 1 notitle, \
     '' every ::1 u 187:188 w l ls 1 notitle, \
     '' every ::1 u 189:190 w l ls 1 notitle, \
     '' every ::1 u 191:192 w l ls 1 notitle, \
     '' every ::1 u 193:194 w l ls 1 notitle, \
     '' every ::1 u 195:196 w l ls 1 notitle, \
     '' every ::1 u 197:198 w l ls 1 notitle, \
     '' every ::1 u 199:200 w l ls 1 title 'Realizations', \
     'diffusion_x_correlations_mean.txt'        u 1:2 w lp ls 2 pt 6 ps 0.45 title 'Mean', \
     'diffusion_x_correlations_mean_fitted.txt' u 1:2 w l  ls 3          title 'Exponential fit'

# ---------------------------------------------------------------------------
# Diffusion Y
# ---------------------------------------------------------------------------
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set output 'diffusion_y_correlation.png'

set xrange [0 : 0.1]
set yrange [0.1 : 1.0]

set xlabel '{/Symbol D}y' font "Arial,32"
set ylabel 'E[{/Symbol w}(y){/Symbol w}(y+{/Symbol D}y)]' font "Arial,32"
set key bottom left font "Arial,24"
set format y "10^{%T}"

plot 'diffusion_y_correlations.txt' every ::1 u 1:2 w l ls 1 notitle, \
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
     '' every ::1 u 39:40 w l ls 1 notitle, \
     '' every ::1 u 41:42 w l ls 1 notitle, \
     '' every ::1 u 43:44 w l ls 1 notitle, \
     '' every ::1 u 45:46 w l ls 1 notitle, \
     '' every ::1 u 47:48 w l ls 1 notitle, \
     '' every ::1 u 49:50 w l ls 1 notitle, \
     '' every ::1 u 51:52 w l ls 1 notitle, \
     '' every ::1 u 53:54 w l ls 1 notitle, \
     '' every ::1 u 55:56 w l ls 1 notitle, \
     '' every ::1 u 57:58 w l ls 1 notitle, \
     '' every ::1 u 59:60 w l ls 1 notitle, \
     '' every ::1 u 61:62 w l ls 1 notitle, \
     '' every ::1 u 63:64 w l ls 1 notitle, \
     '' every ::1 u 65:66 w l ls 1 notitle, \
     '' every ::1 u 67:68 w l ls 1 notitle, \
     '' every ::1 u 69:70 w l ls 1 notitle, \
     '' every ::1 u 71:72 w l ls 1 notitle, \
     '' every ::1 u 73:74 w l ls 1 notitle, \
     '' every ::1 u 75:76 w l ls 1 notitle, \
     '' every ::1 u 77:78 w l ls 1 notitle, \
     '' every ::1 u 79:80 w l ls 1 notitle, \
     '' every ::1 u 81:82 w l ls 1 notitle, \
     '' every ::1 u 83:84 w l ls 1 notitle, \
     '' every ::1 u 85:86 w l ls 1 notitle, \
     '' every ::1 u 87:88 w l ls 1 notitle, \
     '' every ::1 u 89:90 w l ls 1 notitle, \
     '' every ::1 u 91:92 w l ls 1 notitle, \
     '' every ::1 u 93:94 w l ls 1 notitle, \
     '' every ::1 u 95:96 w l ls 1 notitle, \
     '' every ::1 u 97:98 w l ls 1 notitle, \
     '' every ::1 u 99:100 w l ls 1 notitle, \
     '' every ::1 u 101:102 w l ls 1 notitle, \
     '' every ::1 u 103:104 w l ls 1 notitle, \
     '' every ::1 u 105:106 w l ls 1 notitle, \
     '' every ::1 u 107:108 w l ls 1 notitle, \
     '' every ::1 u 109:110 w l ls 1 notitle, \
     '' every ::1 u 111:112 w l ls 1 notitle, \
     '' every ::1 u 113:114 w l ls 1 notitle, \
     '' every ::1 u 115:116 w l ls 1 notitle, \
     '' every ::1 u 117:118 w l ls 1 notitle, \
     '' every ::1 u 119:120 w l ls 1 notitle, \
     '' every ::1 u 121:122 w l ls 1 notitle, \
     '' every ::1 u 123:124 w l ls 1 notitle, \
     '' every ::1 u 125:126 w l ls 1 notitle, \
     '' every ::1 u 127:128 w l ls 1 notitle, \
     '' every ::1 u 129:130 w l ls 1 notitle, \
     '' every ::1 u 131:132 w l ls 1 notitle, \
     '' every ::1 u 133:134 w l ls 1 notitle, \
     '' every ::1 u 135:136 w l ls 1 notitle, \
     '' every ::1 u 137:138 w l ls 1 notitle, \
     '' every ::1 u 139:140 w l ls 1 notitle, \
     '' every ::1 u 141:142 w l ls 1 notitle, \
     '' every ::1 u 143:144 w l ls 1 notitle, \
     '' every ::1 u 145:146 w l ls 1 notitle, \
     '' every ::1 u 147:148 w l ls 1 notitle, \
     '' every ::1 u 149:150 w l ls 1 notitle, \
     '' every ::1 u 151:152 w l ls 1 notitle, \
     '' every ::1 u 153:154 w l ls 1 notitle, \
     '' every ::1 u 155:156 w l ls 1 notitle, \
     '' every ::1 u 157:158 w l ls 1 notitle, \
     '' every ::1 u 159:160 w l ls 1 notitle, \
     '' every ::1 u 161:162 w l ls 1 notitle, \
     '' every ::1 u 163:164 w l ls 1 notitle, \
     '' every ::1 u 165:166 w l ls 1 notitle, \
     '' every ::1 u 167:168 w l ls 1 notitle, \
     '' every ::1 u 169:170 w l ls 1 notitle, \
     '' every ::1 u 171:172 w l ls 1 notitle, \
     '' every ::1 u 173:174 w l ls 1 notitle, \
     '' every ::1 u 175:176 w l ls 1 notitle, \
     '' every ::1 u 177:178 w l ls 1 notitle, \
     '' every ::1 u 179:180 w l ls 1 notitle, \
     '' every ::1 u 181:182 w l ls 1 notitle, \
     '' every ::1 u 183:184 w l ls 1 notitle, \
     '' every ::1 u 185:186 w l ls 1 notitle, \
     '' every ::1 u 187:188 w l ls 1 notitle, \
     '' every ::1 u 189:190 w l ls 1 notitle, \
     '' every ::1 u 191:192 w l ls 1 notitle, \
     '' every ::1 u 193:194 w l ls 1 notitle, \
     '' every ::1 u 195:196 w l ls 1 notitle, \
     '' every ::1 u 197:198 w l ls 1 notitle, \
     '' every ::1 u 199:200 w l ls 1 title 'Realizations', \
     'diffusion_y_correlations_mean.txt'        u 1:2 w lp ls 2 pt 6 ps 0.45 title 'Mean', \
     'diffusion_y_correlations_mean_fitted.txt' u 1:2 w l  ls 3          title 'Exponential fit'
