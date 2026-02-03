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

# Plot: 100 realizations in gray + mean in black
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
     '' using 39:40 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 41:42 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 43:44 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 45:46 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 47:48 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 49:50 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 51:52 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 53:54 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 55:56 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 57:58 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 59:60 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 61:62 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 63:64 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 65:66 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 67:68 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 69:70 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 71:72 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 73:74 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 75:76 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 77:78 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 79:80 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 81:82 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 83:84 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 85:86 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 87:88 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 89:90 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 91:92 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 93:94 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 95:96 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 97:98 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 99:100 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 101:102 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 103:104 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 105:106 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 107:108 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 109:110 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 111:112 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 113:114 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 115:116 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 117:118 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 119:120 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 121:122 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 123:124 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 125:126 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 127:128 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 129:130 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 131:132 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 133:134 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 135:136 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 137:138 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 139:140 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 141:142 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 143:144 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 145:146 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 147:148 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 149:150 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 151:152 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 153:154 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 155:156 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 157:158 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 159:160 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 161:162 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 163:164 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 165:166 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 167:168 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 169:170 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 171:172 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 173:174 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 175:176 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 177:178 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 179:180 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 181:182 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 183:184 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 185:186 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 187:188 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 189:190 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 191:192 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 193:194 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 195:196 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 197:198 with lines lw 1 lc rgb "#CCCCCC" notitle, \
     '' using 199:200 with lines lw 1 lc rgb "#CCCCCC" title 'Realizations', \
     'qx_mean_pdf.txt' using 1:2 with lines lw 3 lc rgb "#000000" title 'Mean PDF'
