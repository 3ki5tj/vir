unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript eps enhanced font "Arial, 20" size 7, 5
set output "virpy.eps"


set rmargin 2.0

set multiplot

pw = 0.5   # width of a panel
ph = 0.5   # height of a panel
px = -0.22  # x-position of the label (a), (b), ...
py = 0.98  # y-position of the label

tagfont="Arial, 32"
lblfont="Arial, 24"
tcfont="Arial, 18"
set key left top Left reverse spacing 1.5

# black
set style line 1 lc rgb "#aaaaaa" lt 1 lw 1.0
set style line 2 lc rgb "#000000" lt 1 lw 2 pt 6 ps 1.2
set style line 3 lc rgb "#000000" lt 1 lw 2 pt 7 ps 1.2

# red
set style line 4 lc rgb "#ffaaaa" lt 2 lw 1.0
set style line 5 lc rgb "#cc0000" lt 2 lw 2 pt 8 ps 1.4
set style line 6 lc rgb "#cc0000" lt 2 lw 2 pt 9 ps 1.4

# blue
set style line 7 lc rgb "#aaaaff" lt 3 lw 1.0
set style line 8 lc rgb "#0000aa" lt 3 lw 2 pt 4 ps 1.
set style line 9 lc rgb "#0000aa" lt 3 lw 2 pt 5 ps 1.




set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
set xlabel "Order {/Arial-Italic n}" \
  offset 2, 0.5 font lblfont

set logscale y

set mytics 5
set ytics 1e-5, 10, 1e7 font tcfont offset 0., 0
set format y '10^{%T}'
set ylabel '{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}' \
  offset 0.0, 0 font lblfont


#set key bottom right


# top left panel
set size pw, ph
set origin 0, ph
set label 1 "(a)" at graph px, py font tagfont


set label 10 "{/Arial-Italic D} = 13" at graph 0.7, 0.1 font lblfont
plot [2:64][1e-4:1e7] \
  "pyhsdata/BnPY13c.dat"            u ($1):(abs($2))                    w l ls 4 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 5 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 6 notitle, \
  "pyhsdata/BnPY13p.dat"            u ($1):(abs($2))                    w l ls 7 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 8 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 9 notitle, \
  "data/D13r1n64/ZrD13r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  1e-100 w lp ls 2 title "simulation", \
  1e-100 w lp ls 8 title "PY, virial", \
  1e-100 w lp ls 5 title "PY, compressibility", \
  1e-100 lw 0 notitle



# top right panel
set size pw, ph
set origin pw, ph
set label 1 "(b)" at graph px, py font tagfont


set label 10 "{/Arial-Italic D} = 15" at graph 0.7, 0.1 font lblfont
plot [2:64][5e-5:2e6] \
  "pyhsdata/BnPY15c.dat"            u ($1):(abs($2))                    w l ls 4 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 5 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 6 notitle, \
  "pyhsdata/BnPY15p.dat"            u ($1):(abs($2))                    w l ls 7 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 8 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 9 notitle, \
  "data/D15r1n64/ZrD15r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  1e-100 w lp ls 2 title "simulation", \
  1e-100 w lp ls 8 title "PY, virial", \
  1e-100 w lp ls 5 title "PY, compressibility", \
  1e-100 lw 0 notitle



# bottom left panel
set size pw, ph
set origin 0, 0
set label 1 "(c)" at graph px, py font tagfont


set label 10 "{/Arial-Italic D} = 17" at graph 0.7, 0.1 font lblfont
plot [2:64][1e-5:4e5] \
  "pyhsdata/BnPY17c.dat"            u ($1):(abs($2))                    w l ls 4 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 5 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 6 notitle, \
  "pyhsdata/BnPY17p.dat"            u ($1):(abs($2))                    w l ls 7 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 8 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 9 notitle, \
  "data/D17r1n64/ZrD17r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  1e-100 w lp ls 2 title "simulation", \
  1e-100 w lp ls 8 title "PY, virial", \
  1e-100 w lp ls 5 title "PY, compressibility", \
  1e-100 lw 0 notitle



# bottom right panel
set size pw, ph
set origin pw, 0
set label 1 "(d)" at graph px, py font tagfont


set label 10 "{/Arial-Italic D} = 19" at graph 0.7, 0.1 font lblfont
plot [2:64][1e-6:1e5] \
  "pyhsdata/BnPY19c.dat"            u ($1):(abs($2))                    w l ls 4 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 5 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 6 notitle, \
  "pyhsdata/BnPY19p.dat"            u ($1):(abs($2))                    w l ls 7 notitle, \
  ""                                u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 8 notitle, \
  ""                                u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 9 notitle, \
  "data/D19r1n64/ZrD19r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  1e-100 w lp ls 2 title "simulation", \
  1e-100 w lp ls 8 title "PY, virial", \
  1e-100 w lp ls 5 title "PY, compressibility", \
  1e-100 lw 0 notitle


unset multiplot

unset output

set terminal wxt
reset



