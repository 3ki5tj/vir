#!/usr/bin/env gnuplot



# Comparison of exact formulae for pressure
# in the PY and HNC closures
# This figure is used in the supplemental material, iegraph.tex
# Remember to check the equation numbers



unset multiplot
reset

set encoding cp1250 # make the minus sign longer
##set encoding iso_8859_1
set terminal postscript eps enhanced size 10, 7 font "Arial, 24"
set output "iepres.eps"

tlfont="Arial, 24"

tcfont="Arial, 24"

# height of the bottom panels
bh = 0.5
# height of the top panels
th = 1 - bh

# width of the right panel
rw = 0.52
# width of the left panel
lw = 1 - rw

#set ytics 10 font tcfont offset 0.3, 0
#set mytics 10
#set format y '10^{%T}'

lbfont  = "Arial, 24"

color1a = "#000000"
# "#dd0000"
color1b = "#002280"

color2a = "#000080"
# "#804000"
color2b = "#000000"

color3a = "#6060ff"
# "#600080"
color3b = "#006000"

color4a = "#400000"
color4b = "#008080"

color5a = "#aa2020"

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 1 lw 1

set style line 2  lc rgb color1a lt 1 lw 2 pt 6  ps 1.8 # empty circle
set style line 3  lc rgb color1a lt 1 lw 2 pt 7  ps 1.8 # full  circle

set style line 4  lc rgb color2a lt 1 lw 2 pt 8  ps 1.6 # empty triangle
set style line 5  lc rgb color2a lt 1 lw 2 pt 9  ps 1.6 # full  triangle

set style line 6  lc rgb color3a lt 1 lw 2 pt 10 ps 1.6 # empty inverted triangle
set style line 7  lc rgb color3a lt 1 lw 2 pt 11 ps 1.6 # full  inverted triangle

set style line 8  lc rgb color4a lt 1 lw 2 pt 4  ps 1.5 # empty square
set style line 9  lc rgb color4a lt 1 lw 2 pt 5  ps 1.5 # full  square

set style line 10 lc rgb color5a lt 1 lw 2 pt 12 ps 1.5 # empty diamond
set style line 11 lc rgb color5a lt 1 lw 2 pt 13 ps 1.5 # full  diamond

set style line 12 lc rgb color3b lt 1 lw 2 pt 14 ps 1.6 # empty pentagon
set style line 13 lc rgb color3b lt 1 lw 2 pt 15 ps 1.6 # full  pentagon

set style line 14 lc rgb color4a lt 1 lw 2 pt 14 ps 2.0 # empty pentagon
set style line 15 lc rgb color4a lt 1 lw 2 pt 15 ps 2.0 # full  pentagon

set style line 16 lc rgb color4b lt 1 lw 2 pt 12 ps 1.8
set style line 17 lc rgb color4b lt 1 lw 2 pt 13 ps 1.8





set size    lw, th
set origin 0.0, bh

set xtics 0.2 font tcfont offset 0, 0
set mxtics 2
thexlabel='Density, {/Symbol-Oblique r}'
set xlabel thexlabel font lbfont offset 0, 0.0

set ytics 2.0 font tcfont
set mytics 2
theylabel='{/Symbol-Oblique b}{/Times-Italic P}/{/Symbol-Oblique r} {/Times - 1}'
set ylabel theylabel font lbfont offset 0.5, 1.0

#set tmargin 1.
#set bmargin 1.5
#set rmargin 0.
#set lmargin 6.0

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
spc = 2.0
set key left Left reverse spacing spc font lbfont

plot [:][0:9] \
  "iedata/thermo/PYhs.dat"            u 1:($2/$1-1)             w l  ls 2  lw 1 notitle, \
  ""                                  u 1:($2/$1-1)  every 2::1 w p  ls 2  t "PY, compressibility-route, Eq. (24)", \
  ""                                  u 1:($3/$1-1)  every 2::1 w p  ls 4  t "PY, compressibility-route, Eq. (51)", \
  ""                                  u 1:($4/$1-1)  every 2::1 w p  ls 6  t "PY, compressibility-route, Eq. (52)", \
  "iedata/thermo/HNChs.dat"           u 1:($2/$1-1)             w l  ls 8  lw 1 notitle, \
  ""                                  u 1:($2/$1-1)  every 2::1 w p  ls 8  t "HNC, virial-route, Eq. (42)", \
  ""                                  u 1:($3/$1-1)  every 2::1 w p  ls 10 t "HNC, virial-route, Eq. (43)", \
  1e-100 lw 0 notitle






unset multiplot
unset output
set terminal wxt
reset



