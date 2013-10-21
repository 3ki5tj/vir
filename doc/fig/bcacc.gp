unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript enhanced font "Arial, 20"
                                #size 10, 4 lw 2.0 font "Arial, 14"
set output "bcacc.ps"



tcfont="Arial, 16"
set key spacing 1.5

set style line 1 lc rgb "#222222" lt 1 lw 3 pt 4 ps 1.5
set style line 2 lc rgb "#000000" lt 2 lw 3 pt 6 ps 1.5
set style line 3 lc rgb "#333333" lt 3 lw 3 pt 8 ps 2.0
set style line 4 lc rgb "#111111" lt 4 lw 3 pt 10 ps 2.0



set logscale y

set mxtics 5
set xtics 5 font tcfont offset 0, 0.3
set xlabel "Dimension {/Arial-Italic D}" offset 0, 0.5

#set mytics 5
set ytics font tcfont offset 0., 0
set format y '10^{%T}'
set ylabel "Probability of forming a biconnected cluster" offset .5, 0
#set key left bottom Left reverse width -0 font "Arial, 12"

plot [:45][1e-8:] \
  "bcacc.txt" u ($1):($5)     w lp ls 1 t "n = 5", \
  ""          u ($1):($10)    w lp ls 2 t "n = 10", \
  ""          u ($1):($15)    w lp ls 3 t "n = 15", \
  ""          u ($1):($20)    w lp ls 4 t "n = 20"

unset output

set terminal wxt
reset



