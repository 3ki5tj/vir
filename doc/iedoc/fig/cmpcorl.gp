#!/usr/bin/env gnuplot
unset multiplot
reset

set encoding cp1250 # make the minus sign longer
##set encoding iso_8859_1
set terminal postscript eps enhanced size 10, 7 font "Arial, 22"
set output "cmpcorl.eps"


tlfont="Arial, 24"
tcfont="Arial, 16"
lbfont  = "Arial, 20"
keyfont = "Arial, 20"



thexlabel='{/Arial-Italic r}'

# height of the bottom panels
bh = 0.5
# height of the top panels
th = 1 - bh

# width of the right panel
rw = 0.52
# width of the left panel
lw = 1 - rw

set ytics 0.1 font tcfont offset 0.3, 0
set mytics 10
#set format y '10^{%T}'

spc = 1.2

color1a = "#cc3333"
color1al = "#ffdddd"
color1b = "#000000"

color2a = "#555588"
color2b = "#606060"

color3a = "#404040"
color3b = "#808080"

color4a = "#606060"
color4b = "#448855"

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 1 lw 1

# Mayer sampling
set style line 2  lc rgb color1al lt 1 lw 5 pt 6  ps 0.5 # light circle
set style line 3  lc rgb color1a  lt 1 lw 5 pt 7  ps 0.5 # full  circle

# DSC, kappa = 0
set style line 4  lc rgb color1b lt 2 lw 6 pt 12 ps 1.1 # empty diamond
set style line 5  lc rgb color1b lt 2 lw 6 pt 13 ps 1.1 # full  diamond

# DSC, kappa != 0
set style line 6  lc rgb color2a lt 4 lw 6 pt 10 ps 1.0 # empty inverted triangle
set style line 7  lc rgb color2a lt 4 lw 6 pt 11 ps 1.0 # full  inverted triangle

# lambda-DSC
set style line 16 lc rgb color4b lt 5 lw 6 pt 12 ps 0.8
set style line 17 lc rgb color4b lt 5 lw 6 pt 13 ps 0.8

# PY
set style line 8  lc rgb color2b lt 3 lw 2 pt 8  ps 1.0 # empty triangle
set style line 9  lc rgb color2b lt 3 lw 2 pt 9  ps 1.0 # full  triangle

# HNC
set style line 10 lc rgb color3a lt 7 lw 2 pt 4  ps 0.8 # empty square
set style line 11 lc rgb color3a lt 7 lw 2 pt 5  ps 0.8 # full  square

set style line 12 lc rgb color3b lt 5 lw 3 pt 14 ps 0.9 # empty pentagon
set style line 13 lc rgb color3b lt 5 lw 3 pt 15 ps 0.9 # full  pentagon

set style line 14 lc rgb color4a lt 6 lw 3 pt 14 ps 0.7 # empty pentagon
set style line 15 lc rgb color4a lt 6 lw 3 pt 15 ps 0.7 # full  pentagon

set style line 100 lc rgb "#808080" lt 1 lw 0.5 # thin solid line


tagdx1 = 0.005
tagdx2 = 0.010
tagdy1 = 0.020
tagdy2 = 0.005
tagfont = "Arial, 24"
set label 300 "(a)" at screen       tagdx1,  1 - tagdy1 font tagfont
set label 301 "(b)" at screen  lw + tagdx2,  1 - tagdy1 font tagfont
set label 302 "(c)" at screen       tagdx1, bh - tagdy2 font tagfont
set label 303 "(d)" at screen  lw + tagdx2, bh - tagdy2 font tagfont



set multiplot


title0 = "Mayer sampling"
title1 = "DSC, {/Symbol-Oblique k}_{/Arial-Italic n} = 0"
title2 = "DSC, {/Symbol-Oblique k}_{/Arial-Italic n} = (0.352)_{{/Arial-Italic n} {/Symbol \263} 4}"
title3 = "{/Symbol-Oblique l}-DSC"
titlepy   = "PY"
titlehnc  = "HNC"



# left-top panel

set size    lw, th
set origin 0.0, bh

set ylabel \
  '{/Arial-Italic c}_{&{/=7 i}2}({/Arial-Italic r&{/=7 i}})' \
  font lbfont offset 3.0, 0.0

set tmargin 1.
set bmargin 1.5
set rmargin 4.
set lmargin 6.0

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 1.57, 1.42 Left reverse spacing spc font keyfont

set xtics 0.2 font tcfont offset 0, 0.5
set mxtics 2
unset xlabel

ymax = 1.5
set yrange [-ymax:ymax]
set ytics nomirror 0.5 font tcfont offset 0.3, 0
set mytics 5

y2max = 0.25
set y2range [-y2max:y2max]
set y2tics nomirror 0.1 font tcfont offset -0.3, 0
set my2tics 10

ps1 = 0.7

colorvsp = "#808080"

set arrow 123 from 1, -ymax to 1, ymax lt 1 lw 5 lc rgb colorvsp nohead



plot [0.7:1.5][:] \
  "iedata/cr/crD6n4.dat"        u ($1):($1 < 1 ? $2 : 1/0)              axes x1y1 w l   ls 3          t title0, \
  ""                            u ($1):($1 > 1 ? $2 : 1/0)              axes x1y2 w l   ls 3          notitle, \
  "iedata/cr/crtrD6samp.dat"    u ($1):(($1 > 1 && $4 == 2) ? -$3 : 1/0)axes x1y1 w l   ls 2          notitle, \
  "iedata/cr/yrtrD6n4.dat"      u ($1):($1 < 1 ? $2 : 1/0)              axes x1y2 w l   ls 2          notitle, \
  "iedata/cr/crtrD6.dat"        u ($1):($1 < 1 && $4 == 2 ? $2 : 1/0)   axes x1y1 w l   ls 5          t title1, \
  ""                            u ($1):($1 > 1 && $4 == 2 ? $2 : 1/0)   axes x1y2 w l   ls 5          notitle, \
  "iedata/cr/crtrD6c0.352.dat"  u ($1):($1 < 1 && $4 == 2 ? $2 : 1/0)   axes x1y1 w l   ls 7          t title2, \
  ""                            u ($1):($1 > 1 && $4 == 2 ? $2 : 1/0)   axes x1y2 w l   ls 7          notitle, \
  "iedata/cr/crtrD6lamc.dat"    u ($1):($1 < 1 && $4 == 2 ? $2 : 1/0)   axes x1y1 w l   ls 17         t title3, \
  ""                            u ($1):($1 > 1 && $4 == 2 ? $2 : 1/0)   axes x1y2 w l   ls 17         notitle, \
  "iedata/cr/crtrD6py.dat"      u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 9          t titlepy, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 9          notitle, \
  "iedata/cr/crtrD6hnc.dat"     u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 11         t titlehnc, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 11         notitle, \
  0 ls 100 notitle





# right-top panel

set size    rw, th
set origin  lw, bh

set rmargin 5.0
set lmargin 6.0

set ylabel \
  '{/Arial-Italic c}_{3&{/=7 i}}({/Arial-Italic r}&{/=7 i})' \
  font lbfont offset 3.0, 0.0

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 1.90, 0.65 Left reverse spacing spc font keyfont

ymax = .7
set yrange [-ymax:ymax]
set ytics nomirror 0.2 font tcfont offset 0.3, 0
set mytics 2

y2max = 0.15
set y2range [-y2max:y2max]
set y2tics nomirror 0.1 font tcfont offset -0.3, 0
set my2tics 10

set arrow 123 from 1, -ymax to 1, ymax lt 1 lw 5 lc rgb colorvsp nohead



plot [0.6:1.8][:] \
  "iedata/cr/crD6n5.dat"        u ($1):($1 < 1 ? $2 : 1/0)              axes x1y1 w l   ls 3          t title0, \
  ""                            u ($1):($1 > 1 ? $2 : 1/0)              axes x1y2 w l   ls 3          notitle, \
  "iedata/cr/crtrD6samp.dat"    u ($1):(($1 > 1 && $4 == 3) ? -$3 : 1/0)axes x1y1 w l   ls 2          notitle, \
  "iedata/cr/yrtrD6n5.dat"      u ($1):($1 < 1 ? $2 : 1/0)              axes x1y2 w l   ls 2          notitle, \
  "iedata/cr/crtrD6.dat"        u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 5          t title1, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 5          notitle, \
  "iedata/cr/crtrD6c0.352.dat"  u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 7          t title2, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 7          notitle, \
  "iedata/cr/crtrD6lamc.dat"    u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 17         t title3, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 17         notitle, \
  "iedata/cr/crtrD6py.dat"      u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 9          t titlepy, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 9          notitle, \
  "iedata/cr/crtrD6hnc.dat"     u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 11         t titlehnc, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 11         notitle, \
  0 ls 100 notitle





# left-bottom panel

set size    lw, bh
set origin 0.0, 0.0

set tmargin 0.
set bmargin 2.5
set xlabel thexlabel font lbfont offset 0, 1.0

set lmargin 6.0
set ylabel \
  '{/Arial-Italic c}_{4&{/=7 i}}({/Arial-Italic r}&{/=7 i})' \
  font lbfont offset 3.5, 0.0

set rmargin 4.

# Left: align text to the left
# reverse: symbol first, texticc -DD=6 -DN=3 ../../mcrat0.c && time ./a.out --nstcr=20 --nstcrrep=100000 -11e10 next
# invert: first drawn shown last in the legend
set key at 2.10, 0.48 Left reverse spacing spc font keyfont

#
ymax = 0.5
set yrange [-ymax:ymax]
set ytics nomirror 0.2 font tcfont offset 0.3, 0
set mytics 2

y2max = 0.25
set y2range [-y2max:y2max]
set y2tics nomirror 0.1 font tcfont offset -0.3, 0
set my2tics 10

set arrow 123 from 1, -ymax to 1, ymax lt 1 lw 5 lc rgb colorvsp nohead



plot [0.4:2.0][:] \
  "iedata/cr/crD6n6.dat"        u ($1):($1 < 1 ? $2 : 1/0)              axes x1y1 w l   ls 3          t title0, \
  ""                            u ($1):($1 > 1 ? $2 : 1/0)              axes x1y2 w l   ls 3          notitle, \
  "iedata/cr/crtrD6samp.dat"    u ($1):(($1 > 1 && $4 == 4) ? -$3 : 1/0)axes x1y1 w l   ls 2          notitle, \
  "iedata/cr/yrtrD6n6.dat"      u ($1):($1 < 1 ? $2 : 1/0)              axes x1y2 w l   ls 2          notitle, \
  "iedata/cr/crtrD6.dat"        u ($1):($1 < 1 && $4 == 4 ? $2 : 1/0)   axes x1y1 w l   ls 5          t title1, \
  ""                            u ($1):($1 > 1 && $4 == 4 ? $2 : 1/0)   axes x1y2 w l   ls 5          notitle, \
  "iedata/cr/crtrD6c0.352.dat"  u ($1):($1 < 1 && $4 == 4 ? $2 : 1/0)   axes x1y1 w l   ls 7          t title2, \
  ""                            u ($1):($1 > 1 && $4 == 4 ? $2 : 1/0)   axes x1y2 w l   ls 7          notitle, \
  "iedata/cr/crtrD6lamc.dat"    u ($1):($1 < 1 && $4 == 4 ? $2 : 1/0)   axes x1y1 w l   ls 17         t title3, \
  ""                            u ($1):($1 > 1 && $4 == 4 ? $2 : 1/0)   axes x1y2 w l   ls 17         notitle, \
  "iedata/cr/crtrD6py.dat"      u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 9          t titlepy, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 9          notitle, \
  "iedata/cr/crtrD6hnc.dat"     u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 11         t titlehnc, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 11         notitle, \
  0 ls 100 notitle





unset arrow




# right-bottom panel

set size    rw, bh
set origin  lw, 0.0

set lmargin 6.0
set rmargin 5.0

set ylabel \
  '{/Arial-Italic c}_{5&{/=7 i}}({/Arial-Italic r}&{/=7 i})' \
  font lbfont offset 3.0, 0.0

set ytics 10 font tcfont offset 0.3, 0
set mytics 10

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 2.32, 0.47 Left reverse spacing spc font keyfont

ymax = 0.50
set yrange [-ymax:ymax]
set ytics nomirror 0.2 font tcfont offset 0.3, 0
set mytics 2

y2max = 0.25
set y2range [-y2max:y2max]
set y2tics nomirror 0.1 font tcfont offset -0.3, 0
set my2tics 10

set arrow 123 from 1, -ymax to 1, ymax lt 1 lw 5 lc rgb colorvsp nohead



plot [0.4:2.2][:] \
  "iedata/cr/crD6n7.dat"        u ($1):($1 < 1 ? $2 : 1/0)              axes x1y1 w l   ls 3          t title0, \
  ""                            u ($1):($1 > 1 ? $2 : 1/0)              axes x1y2 w l   ls 3          notitle, \
  "iedata/cr/crtrD6samp.dat"    u ($1):(($1 > 1 && $4 == 5) ? -$3 : 1/0)axes x1y1 w l   ls 2          notitle, \
  "iedata/cr/yrtrD6n7.dat"      u ($1):($1 < 1 ? $2 : 1/0)              axes x1y2 w l   ls 2          notitle, \
  "iedata/cr/crtrD6.dat"        u ($1):($1 < 1 && $4 == 5 ? $2 : 1/0)   axes x1y1 w l   ls 5          t title1, \
  ""                            u ($1):($1 > 1 && $4 == 5 ? $2 : 1/0)   axes x1y2 w l   ls 5          notitle, \
  "iedata/cr/crtrD6c0.352.dat"  u ($1):($1 < 1 && $4 == 5 ? $2 : 1/0)   axes x1y1 w l   ls 7          t title2, \
  ""                            u ($1):($1 > 1 && $4 == 5 ? $2 : 1/0)   axes x1y2 w l   ls 7          notitle, \
  "iedata/cr/crtrD6lamc.dat"    u ($1):($1 < 1 && $4 == 5 ? $2 : 1/0)   axes x1y1 w l   ls 17         t title3, \
  ""                            u ($1):($1 > 1 && $4 == 5 ? $2 : 1/0)   axes x1y2 w l   ls 17         notitle, \
  "iedata/cr/crtrD6py.dat"      u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 9          t titlepy, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 9          notitle, \
  "iedata/cr/crtrD6hnc.dat"     u ($1):($1 < 1 && $4 == 3 ? $2 : 1/0)   axes x1y1 w l   ls 11         t titlehnc, \
  ""                            u ($1):($1 > 1 && $4 == 3 ? $2 : 1/0)   axes x1y2 w l   ls 11         notitle, \
  0 ls 100 notitle





unset multiplot
unset output
set terminal wxt
reset



