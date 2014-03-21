unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript eps enhanced font "Arial, 18"
                                #size 10, 4 lw 2.0 font "Arial, 14"
set output "vir64.eps"



tcfont="Arial, 16"
set key spacing 1.5

set style line 1 lc rgb "#aaaaaa" lt 1 lw 0.5

set style line 2 lt 1 lw 1 pt 6 ps 1.0
set style line 3 lt 1 lw 1 pt 7 ps 1.0

set style line 4 lc rgb "#000080" lt 1 lw 1 pt 6 ps 1.0
set style line 5 lc rgb "#000080" lt 1 lw 1 pt 7 ps 1.0

set style line 6 lc rgb "#008000" lt 1 lw 1 pt 6 ps 1.0
set style line 7 lc rgb "#008000" lt 1 lw 1 pt 7 ps 1.0

set style line 8 lc rgb "#aa0000" lt 1 lw 1 pt 6 ps 1.0
set style line 9 lc rgb "#aa0000" lt 1 lw 1 pt 7 ps 1.0




set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
set xlabel "Order {/Arial-Italic n}" offset 2, 0.5

set logscale y

set ytics font tcfont offset 0., 0
set format y '10^{%T}'
set ylabel '{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}' offset 0.5, 0
#set key left bottom Left reverse width -0 font "Arial, 12"

set rmargin 7

lbfont = "Arial, 16"
xlbl = 64.5
set label "{/Arial-Italic D} = 15"  at xlbl, 2e6   font lbfont
set label "{/Arial-Italic D} = 20"  at xlbl, 3e3   font lbfont
set label "{/Arial-Italic D} = 25"  at xlbl, 3e0   font lbfont
set label "{/Arial-Italic D} = 30"  at xlbl, 3e-3  font lbfont
set label "{/Arial-Italic D} = 35"  at xlbl, 5e-7  font lbfont
set label "{/Arial-Italic D} = 40"  at xlbl, 3e-10 font lbfont
set label "{/Arial-Italic D} = 45"  at xlbl, 1e-13 font lbfont
set label "{/Arial-Italic D} = 50"  at xlbl, 2e-17 font lbfont
set label "{/Arial-Italic D} = 60"  at xlbl, 3e-24 font lbfont
set label "{/Arial-Italic D} = 70"  at xlbl, 1e-30 font lbfont
set label "{/Arial-Italic D} = 80"  at xlbl, 1e-37 font lbfont
set label "{/Arial-Italic D} = 90"  at xlbl, 1e-44 font lbfont
set label "{/Arial-Italic D} = 100" at xlbl, 1e-51 font lbfont

plot [2:64][1e-55:] \
  "data/D15r1n64/ZrD15r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  "data/D20r1n64/ZrD20r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 4 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 5 notitle, \
  "data/D25r1n64/ZrD25r1n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 6 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 7 notitle, \
  "data/D30r2n64/ZrD30r2n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 8 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 9 notitle, \
  "data/D35r2n64/ZrD35r2n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  "data/D40r3n64/ZrD40r3n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 4 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 5 notitle, \
  "data/D45r3n64/ZrD45r3n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 6 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 7 notitle, \
  "data/D50r3n64/ZrD50r3n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 8 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 9 notitle, \
  "data/D60r3n64/ZrD60r3n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  "data/D70r4n64/ZrD70r4n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 4 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 5 notitle, \
  "data/D80r4n64/ZrD80r4n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 6 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 7 notitle, \
  "data/D90r4n64/ZrD90r4n64.data"   u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 8 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 9 notitle, \
  "data/D100r4n64/ZrD100r4n64.data" u ($1):(abs($19))                   w l ls 1 notitle, \
  ""                                u ($1):(($19 > 0) ? abs($19) : 1/0) w p ls 2 notitle, \
  ""                                u ($1):(($19 < 0) ? abs($19) : 1/0) w p ls 3 notitle, \
  1e-100 lw 0 notitle

unset output

set terminal wxt
reset



