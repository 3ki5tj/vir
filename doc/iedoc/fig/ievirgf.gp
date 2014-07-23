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

spc = 1.5

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

set style line 2  lc rgb color1a lt 1 pt 4  ps 0.8 # empty square
set style line 3  lc rgb color1a lt 1 pt 5  ps 0.8 # full  square

set style line 4  lc rgb color1b lt 2 pt 12 ps 1.1 # empty diamond
set style line 5  lc rgb color1b lt 2 pt 13 ps 1.1 # full  diamond

set style line 6  lc rgb color2a lt 3 pt 10 ps 1.0 # empty inverted triangle
set style line 7  lc rgb color2a lt 3 pt 11 ps 1.0 # full  inverted triangle

set style line 8  lc rgb color2b lt 4 pt 8  ps 1.0 # empty triangle
set style line 9  lc rgb color2b lt 4 pt 9  ps 1.0 # full  triangle

set style line 10 lc rgb color3a lt 5 pt 6  ps 0.8 # empty circle
set style line 11 lc rgb color3a lt 5 pt 7  ps 0.8 # full  circle

set style line 12 lc rgb color3b lt 6 pt 14 ps 0.9 # empty pentagon
set style line 13 lc rgb color3b lt 6 pt 15 ps 0.9 # full  pentagon

set style line 14 lc rgb color4a lt 7 pt 14 ps 0.7 # empty pentagon
set style line 15 lc rgb color4a lt 7 pt 15 ps 0.7 # full  pentagon

set style line 16 lc rgb color4b lt 8 pt 12 ps 0.8
set style line 17 lc rgb color4b lt 8 pt 13 ps 0.8





set size    lw, th
set origin 0.0, bh

set xtics 10 font tcfont offset 0, 0.5
set mxtics 10
unset xlabel

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
set key Left reverse spacing spc font lbfont

plot [3:12][1e-3:] \
  "gfdata/gaussfD2mpf.dat"                u ($1):(abs($2))                              w l  ls 2  lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 2  lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 3  lw 1.0 notitle, \
  ""                                      u ($1):-1                                     w lp ls 2  lw 1.0 t "{/Arial-Italic D} = 2, Exact", \
  "gfdata/gaussfD4mpf.dat"                u ($1):(abs($2))                              w l  ls 6  lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 6  lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 7  lw 1.0 notitle, \
  ""                                      u ($1):-1                                     w lp ls 6  lw 1.0 t "{/Arial-Italic D} = 4, Exact", \
  "gfdata/gaussfD10mpf.dat"               u ($1):(abs($2))                              w l  ls 10 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 10 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 11 lw 1.0 notitle, \
  ""                                      u ($1):-1                                     w lp ls 10 lw 1.0 t "{/Arial-Italic D} = 10, Exact", \
  1e-100 lw 0 notitle






unset multiplot
unset output
set terminal wxt
reset



