#!/usr/bin/env gnuplot
unset multiplot
reset

set encoding cp1250 # make the minus sign longer
##set encoding iso_8859_1
set terminal postscript eps enhanced size 10, 7 font "Arial, 22"
set output "ievirgf.eps"

tlfont="Arial, 24"

tcfont="Arial, 16"

thexlabel='Order, {/Arial-Italic n}'

# height of the bottom panels
bh = 0.5
# height of the top panels
th = 1 - bh

# width of the right panel
rw = 0.52
# width of the left panel
lw = 1 - rw

set logscale y
set ytics 10 font tcfont offset 0.3, 0
set mytics 10
set format y '10^{%T}'

lbfont  = "Arial, 20"

color1a = "#dd0000"
color1b = "#002280"

color2a = "#804000"
color2b = "#000000"

color3a = "#600080"
color3b = "#006000"

color4a = "#606060"
color4b = "#008080"

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 1 lw 1

set style line 2  lc rgb color1a lt 1 pt 4  ps 1.6 # empty square
set style line 3  lc rgb color1a lt 1 pt 5  ps 1.6 # full  square

set style line 4  lc rgb color1b lt 2 pt 12 ps 2.2 # empty diamond
set style line 5  lc rgb color1b lt 2 pt 13 ps 2.2 # full  diamond

set style line 6  lc rgb color2a lt 3 pt 10 ps 2.2 # empty inverted triangle
set style line 7  lc rgb color2a lt 3 pt 11 ps 2.2 # full  inverted triangle

set style line 8  lc rgb color2b lt 4 pt 8  ps 2.2 # empty triangle
set style line 9  lc rgb color2b lt 4 pt 9  ps 2.2 # full  triangle

set style line 10 lc rgb color3a lt 5 pt 6  ps 1.6 # empty circle
set style line 11 lc rgb color3a lt 5 pt 7  ps 1.6 # full  circle

set style line 12 lc rgb color3b lt 6 pt 14 ps 1.8 # empty pentagon
set style line 13 lc rgb color3b lt 6 pt 15 ps 1.8 # full  pentagon

set style line 14 lc rgb color4a lt 7 pt 14 ps 1.6 # empty pentagon
set style line 15 lc rgb color4a lt 7 pt 15 ps 1.6 # full  pentagon

set style line 16 lc rgb color4b lt 8 pt 12 ps 1.8
set style line 17 lc rgb color4b lt 8 pt 13 ps 1.8





set size    lw, th
set origin 0.0, bh

set xtics 1 font tcfont offset 0, 0.5

set xlabel thexlabel font lbfont offset 0, 1

set ylabel \
  '{/Arial-Italic B}_n /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}' \
  font lbfont offset 1.5, -0.5

#set tmargin 1.
#set bmargin 1.5
#set rmargin 0.
#set lmargin 6.0

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
spc = 1.3
set key at 14, 0.9 Left reverse spacing spc font lbfont

plot [3:12][5e-3:1] \
  "gfdata/gaussfD2mpf.dat"                        u ($1):(abs($2))                              w l  ls 2  lw 0.5 notitle, \
  ""                                              u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 2  lw 1.0 notitle, \
  ""                                              u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 3  lw 1.0 notitle, \
  ""                                              u ($1):-1                                     w lp ls 2  lw 1.0 t "{/Arial-Italic D} = 2, Exact", \
  "iedata/gaussf/D2n12R16c-0.28/xBnPYcD2n12.dat"  u ($1):(abs($2))                              w l  ls 4  lw 0.5 notitle, \
  ""                                              u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 4  lw 1.0 notitle, \
  ""                                              u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 5  lw 1.0 notitle, \
  ""                                              u ($1):-1                                     w lp ls 4  lw 1.0 t "{/Arial-Italic D} = 2, Self-consistent, {/Symbol-Oblique k} = 0.28", \
  "iedata/gaussf/D2n12R16/xBnPYcD2n12.dat"        u ($1):(abs($2))                              w l  ls 10 lw 0.5 notitle, \
  ""                                              u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 10 lw 1.0 notitle, \
  ""                                              u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 11 lw 1.0 notitle, \
  ""                                              u ($1):-1                                     w lp ls 10 lw 1.0 t "{/Arial-Italic D} = 2, Self-consistent, {/Symbol-Oblique k} = 0", \
  "gfdata/gaussfD8mpf.dat"                        u ($1):(abs($2))                              w l  ls 6  lw 0.5 notitle, \
  ""                                              u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 6  lw 1.0 notitle, \
  ""                                              u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 7  lw 1.0 notitle, \
  ""                                              u ($1):-1                                     w lp ls 6  lw 1.0 t "{/Arial-Italic D} = 8, Exact", \
  "iedata/gaussf/D8n12R16/xBnPYcD8n12.dat"        u ($1):(abs($2))                              w l  ls 8  lw 0.5 notitle, \
  ""                                              u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 8  lw 1.0 notitle, \
  ""                                              u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 9  lw 1.0 notitle, \
  ""                                              u ($1):-1                                     w lp ls 8  lw 1.0 t "{/Arial-Italic D} = 8, Self-consistent, {/Symbol-Oblique k} = 0", \
  1e-100 lw 0 notitle






unset multiplot
unset output
set terminal wxt
reset



