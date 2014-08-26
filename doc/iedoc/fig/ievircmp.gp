unset multiplot
reset

set encoding cp1250 # make minus sign longer
##set encoding iso_8859_1
set terminal postscript eps enhanced size 7, 9 font "Arial, 20"
set output "ievircmp.eps"

tlfont="Arial, 24"

tcfont="Arial, 16"
thexlabel='Order {/Arial-Italic n}'
theylabel='{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}'


titledsc = "DSC"

# height of the bottom panels
bh = 0.5
# height of the top panels
th = 1 - bh

# width of the right panel
rw = 0.52
# width of the left panel
lw = 1 - rw

spc = 1.2

set logscale y
set format y '10^{%T}'
set ytics 1e-8,10,0.1 font tcfont offset 0.3, 0
set ytics add ("1 " 1.0)
set mytics 10

lbfont  = "Arial, 20"

color1a = "#cc2222"
color1b = "#000000"

color2a = "#666666"  # "#804000"
color2b = "#666666"

color3a = "#666666"  # "#600080"
color3b = "#666666"  # "#006000"

color4a = "#666666"  # "#806000"
color4b = "#666666"  # "#ccaa44"

color5a = "#666666"  # "#608080"
color5b = "#666666"  # "#44cccc"

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 1 lw 1

set style line 2  lc rgb color1a lt 1 lw 4.0 pt 4  ps 1.4 # empty square
set style line 3  lc rgb color1a lt 1 lw 4.0 pt 5  ps 1.4 # full  square

set style line 4  lc rgb color1b lt 2 lw 4.0 pt 12 ps 2.0 # empty diamond
set style line 5  lc rgb color1b lt 2 lw 4.0 pt 13 ps 2.0 # full  diamond

set style line 6  lc rgb color2a lt 4 lw 2.0 pt 10 ps 1.7 # empty inverted triangle
set style line 7  lc rgb color2a lt 4 lw 2.0 pt 11 ps 1.7 # full  inverted triangle

set style line 8  lc rgb color2b lt 4 lw 2.0 pt 8  ps 1.7 # empty triangle
set style line 9  lc rgb color2b lt 4 lw 2.0 pt 9  ps 1.7 # full  triangle

set style line 10 lc rgb color3a lt 4 lw 2.0 pt 6  ps 1.4 # empty circle
set style line 11 lc rgb color3a lt 4 lw 2.0 pt 7  ps 1.4 # full  circle

set style line 12 lc rgb color3b lt 4 lw 2.0 pt 14 ps 1.6 # empty pentagon
set style line 13 lc rgb color3b lt 4 lw 2.0 pt 15 ps 1.6 # full  pentagon

set style line 14 lc rgb color4a lt 4 lw 2.0 pt 3 ps 1.3 # empty pentagon
set style line 15 lc rgb color4a lt 4 lw 6.0 pt 3 ps 1.3 # full  pentagon

set style line 16 lc rgb color4b lt 4 lw 2.0 pt 2 ps 1.3 # empty diamond
set style line 17 lc rgb color4b lt 4 lw 9.0 pt 2 ps 1.3 # full  diamond

set style line 18 lc rgb color5a lt 4 lw 3.0 pt 1 ps 1.6 # empty diamond
set style line 19 lc rgb color5a lt 4 lw 9.0 pt 1 ps 1.6 # full  diamond

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




# left-top panel

set size    lw, th
set origin 0.0, bh

set xtics 4 font tcfont offset 0, 0.5
set mxtics 4
unset xlabel

set ylabel theylabel font lbfont offset 1.5, -1.5

set tmargin 1.
set bmargin 1.5
set lmargin 6.0
set rmargin 0.5

set label 100 "{/Arial-Italic D} = 2" at 13.5, 7e-5 font tlfont

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 13.5, 10e-4 Left reverse spacing spc font lbfont

plot [2:16][5e-5:1] \
  "data/D2/BnD2n14.dat"                   u ($1):(abs($2))                              w l  ls 2         notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 2         notitle, \
  ""                                      u ($1):-1                                     w lp ls 2         t "Mayer sampling", \
  "iedata/xBnPYcD2n32.dat"                u ($1):(abs($2))                              w l  ls 4         notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 4         notitle, \
  ""                                      u ($1):-1                                     w lp ls 4         t titledsc, \
  "iedata/py/xBnPYD2n32.dat"              u ($1):(abs($3))                              w l  ls 10 lw 0.5 notitle, \
  ""                                      u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 10        notitle, \
  ""                                      u ($1):-1                                     w lp ls 10        t "PY, virial", \
  ""                                      u ($1):(abs($2))                              w l  ls 12 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 12        notitle, \
  ""                                      u ($1):-1                                     w lp ls 12        t "PY, compressibility", \
  ""                                      u ($1):(abs($4))                              w l  ls 18 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 18        notitle, \
  ""                                      u ($1):-1                                     w lp ls 18        t "PY, {/Symbol-Oblique c}", \
  "iedata/hnc/xBnHNCD2n32.dat"            u ($1):(abs($3))                              w l  ls 6  lw 0.5 notitle, \
  ""                                      u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 6         notitle, \
  ""                                      u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 7         notitle, \
  ""                                      u ($1):-1                                     w lp ls 6         t "HNC, virial", \
  ""                                      u ($1):(abs($2))                              w l  ls 8  lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 8         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 9         notitle, \
  ""                                      u ($1):-1                                     w lp ls 8         t "HNC, compressibility", \
  1e-100 lw 0 notitle






# right-top panel

set size    rw, th
set origin  lw, bh

set lmargin 7.0
set rmargin 1.5

set ylabel theylabel font lbfont offset 1.5, 0.0


set label 100 "{/Arial-Italic D} = 7" at 26, 8e-6 font tlfont

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 27, 3.5 Left reverse spacing spc font lbfont

plot [2:32][4e-6:4] \
  "data/D7/BnD7n20.dat"                   u ($1):(abs($2))                              w l  ls 2         notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 2         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 3         notitle, \
  ""                                      u ($1):-1                                     w lp ls 2         t "Mayer sampling", \
  "iedata/xBnPYcD7n128.dat"               u ($1):(abs($2))                              w l  ls 4         notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 4         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 5         notitle, \
  ""                                      u ($1):-1                                     w lp ls 4         t titledsc, \
  "iedata/py/BnPYD7n32R34M65536f128.dat"  u ($1):(abs($3))                              w l  ls 10 lw 0.5 notitle, \
  ""                                      u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 10        notitle, \
  ""                                      u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 11        notitle, \
  ""                                      u ($1):-1                                     w lp ls 10        t "PY, virial", \
  ""                                      u ($1):(abs($2))                              w l  ls 12 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 12        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 13        notitle, \
  ""                                      u ($1):-1                                     w lp ls 12        t "PY, compressibility", \
  ""                                      u ($1):(abs($4))                              w l  ls 18 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 18        notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)             w p  ls 19        notitle, \
  ""                                      u ($1):-1                                     w lp ls 18        t "PY, {/Symbol-Oblique c}", \
  "iedata/pyhnc/xBnHNCD7n36.dat"          u ($1):(abs($3))                              w l  ls 6  lw 0.5 notitle, \
  ""                                      u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 6         notitle, \
  ""                                      u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 7         notitle, \
  ""                                      u ($1):-1                                     w lp ls 6         t "HNC, virial", \
  ""                                      u ($1):(abs($2))                              w l  ls 8  lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 8         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 9         notitle, \
  ""                                      u ($1):-1                                     w lp ls 8         t "HNC, compressibility", \
  ""                                      u ($1):(abs($4))                              w l  ls 16 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 16        notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)             w p  ls 17        notitle, \
  ""                                      u ($1):-1                                     w lp ls 16        t "HNC, cavity", \
  1e-100 lw 0 notitle





# left-bottom panel

set size    lw, bh
set origin 0.0, 0.0

set tmargin 0.
set bmargin 2.5
set xlabel thexlabel font lbfont offset 2, 1.0

set lmargin 6.0
set rmargin 0.5
set ylabel theylabel font lbfont offset 1.5, -0.5

set label 100 "{/Arial-Italic D} = 10" at 26.0, 8.0e-4 font tlfont

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 26.5, 3.5 Left reverse spacing spc font lbfont

plot [2:32][5e-4:4] \
  "data/D10r1n32/BnD10n32.dat"              u ($1):(abs($2))                              w l  ls 2         notitle, \
  ""                                        u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 2         notitle, \
  ""                                        u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 3         notitle, \
  ""                                        u ($1):-1                                     w lp ls 2         t "Mayer sampling", \
  "iedata/xBnPYcD10n128.dat"                u ($1):(abs($2))                              w l  ls 4         notitle, \
  ""                                        u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 4         notitle, \
  ""                                        u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 5         notitle, \
  ""                                        u ($1):-1                                     w lp ls 4         t titledsc, \
  "iedata/py/xBnPYD10n64.dat"               u ($1):(abs($3))                              w l  ls 10 lw 0.5 notitle, \
  ""                                        u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 10        notitle, \
  ""                                        u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 11        notitle, \
  ""                                        u ($1):-1                                     w lp ls 10        t "PY, virial", \
  ""                                        u ($1):(abs($2))                              w l  ls 12 lw 0.5 notitle, \
  ""                                        u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 12        notitle, \
  ""                                        u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 13        notitle, \
  ""                                        u ($1):-1                                     w lp ls 12        t "PY, compressibility", \
  ""                                        u ($1):(abs($4))                              w l  ls 18 lw 0.5 notitle, \
  ""                                        u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 18        notitle, \
  ""                                        u ($1):(($4 < 0) ? abs($4) : 1/0)             w p  ls 19        notitle, \
  ""                                        u ($1):-1                                     w lp ls 18 lw 3.0 t "PY, {/Symbol-Oblique c}", \
  "iedata/pyhnc/xBnHNCD10n36.dat"           u ($1):(abs($3))                              w l  ls 6         notitle, \
  ""                                        u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 6         notitle, \
  ""                                        u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 7         notitle, \
  ""                                        u ($1):-1                                     w lp ls 6         t "HNC, virial", \
  ""                                        u ($1):(abs($2))                              w l  ls 8  lw 0.5 notitle, \
  ""                                        u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 8         notitle, \
  ""                                        u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 9         notitle, \
  ""                                        u ($1):-1                                     w lp ls 8         t "HNC, compressibility", \
  ""                                        u ($1):(abs($4))                              w l  ls 16 lw 0.5 notitle, \
  ""                                        u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 16        notitle, \
  ""                                        u ($1):(($4 < 0) ? abs($4) : 1/0)             w p  ls 17        notitle, \
  ""                                        u ($1):-1                                     w lp ls 16        t "HNC, cavity", \
  1e-100 lw 0 notitle




unset arrow




# right-bottom panel

set size    rw, bh
set origin  lw, 0.0

set lmargin 7.0
set rmargin 1.5


set ylabel theylabel font lbfont offset 1.5, 4.5

set label 100 "{/Arial-Italic D} = 15" at 26, 1.5e-4 font tlfont

# Left: align text to the left
# reverse: symbol first, text next
# invert: first drawn shown last in the legend
set key at 26.5, 0.9 Left reverse spacing spc font lbfont

plot [2:32][1e-4:1] \
  "data/D15r1n64/BnD15n64.dat"            u ($1):(abs($2))                              w l  ls 2         notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 2         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 3         notitle, \
  ""                                      u ($1):-1                                     w lp ls 2         t "Mayer sampling", \
  "iedata/xBnPYcD15n128.dat"              u ($1):(abs($2))                              w l  ls 4         notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 4         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 5         notitle, \
  ""                                      u ($1):-1                                     w lp ls 4         t titledsc, \
  "iedata/pyhnc/xBnPYD15n36.dat"          u ($1):(abs($3))                              w l  ls 10 lw 0.5 notitle, \
  ""                                      u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 10        notitle, \
  ""                                      u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 11        notitle, \
  ""                                      u ($1):-1                                     w lp ls 10        t "PY, virial", \
  ""                                      u ($1):(abs($2))                              w l  ls 12 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 12        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 13        notitle, \
  ""                                      u ($1):-1                                     w lp ls 12        t "PY, compressibility", \
  ""                                      u ($1):(abs($4))                              w l  ls 18 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 18        notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)             w p  ls 19        notitle, \
  ""                                      u ($1):-1                                     w lp ls 18        t "PY, {/Symbol-Oblique c}", \
  "iedata/pyhnc/xBnHNCD15n36.dat"         u ($1):(abs($3))                              w l  ls 6  lw 0.5 notitle, \
  ""                                      u ($1):(($3 > 0) ? abs($3) : 1/0)             w p  ls 6         notitle, \
  ""                                      u ($1):(($3 < 0) ? abs($3) : 1/0)             w p  ls 7         notitle, \
  ""                                      u ($1):-1                                     w lp ls 6         t "HNC, virial", \
  ""                                      u ($1):(abs($2))                              w l  ls 8  lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)             w p  ls 8         notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)             w p  ls 9         notitle, \
  ""                                      u ($1):-1                                     w lp ls 8         t "HNC, compressibility", \
  ""                                      u ($1):(abs($4))                              w l  ls 16 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)             w p  ls 16        notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)             w p  ls 17        notitle, \
  ""                                      u ($1):-1                                     w lp ls 16        t "HNC, cavity", \
  1e-100 lw 0 notitle

unset multiplot
unset output
set terminal wxt
reset



