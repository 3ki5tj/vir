unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript eps enhanced size 7, 10 font "Arial, 20"
set output "ievirlow.eps"

tcfont="Arial, 16"
thexlabel='Order {/Arial-Italic n}'
theylabel='{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}'

# height of the bottom panels
bh = 0.55
# height of the top panels
th = 1 - bh

# width of the right panel
rw = 0.52
# width of the left panel
lw = 1 - rw

set logscale y
set ytics font tcfont offset 0.3, 0
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

set style line 2  lc rgb color1a lt 1 pt 4  ps 2.0 # empty square
set style line 3  lc rgb color1a lt 1 pt 5  ps 2.0 # full  square

set style line 4  lc rgb color1b lt 1 pt 12 ps 2.6 # empty diamond
set style line 5  lc rgb color1b lt 1 pt 13 ps 2.6 # full  diamond

set style line 6  lc rgb color2a lt 1 pt 10 ps 2.4 # empty inverted triangle
set style line 7  lc rgb color2a lt 1 pt 11 ps 2.4 # full  inverted triangle

set style line 8  lc rgb color2b lt 1 pt 8  ps 2.4 # empty triangle
set style line 9  lc rgb color2b lt 1 pt 9  ps 2.4 # full  triangle

set style line 10 lc rgb color3a lt 1 pt 6  ps 2.0 # empty circle
set style line 11 lc rgb color3a lt 1 pt 7  ps 2.0 # full  circle

set style line 12 lc rgb color3b lt 1 pt 14 ps 2.2 # empty pentagon
set style line 13 lc rgb color3b lt 1 pt 15 ps 2.2 # full  pentagon

set style line 14 lc rgb color4a lt 1 pt 4  ps 2.0
set style line 15 lc rgb color4a lt 1 pt 5  ps 2.0

set style line 16 lc rgb color4b lt 1 pt 12 ps 2.4
set style line 17 lc rgb color4b lt 1 pt 13 ps 2.4



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

set ylabel theylabel offset 1.3, 0.7

set tmargin 1.
set bmargin 1.5
set rmargin 0.
set lmargin 6.0

set label 101 "{/Arial-Italic D} = 2" at 18.0, 3.0e-4 rotate by -55  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 3" at 12.0, 1.5e-4 rotate by -72  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 4" at  5.0, 5.5e-3 rotate by -78  textcolor rgb color3b font lbfont

plot [2:22][8e-6:1] \
  "data/D2/BnD2n14.dat"                   u ($1):(abs($2)):3                            w l  ls 2  lt 1 lw 0.5 notitle, \
  ""                                      u ($1):(abs($2)):3                            w e  ls 2       lw 3.0 notitle, \
  "data/D3/BnD3n12.dat"                   u ($1):(abs($2)):3                            w l  ls 6  lt 2 lw 0.5 notitle, \
  ""                                      u ($1):(abs($2)):3                            w e  ls 6       lw 3.0 notitle, \
  "data/D4/BnD4n11.dat"                   u ($1):(abs($2)):3                            w l  ls 10 lt 4 lw 0.5 notitle, \
  ""                                      u ($1):(abs($2)):3                            w e  ls 10      lw 3.0 notitle, \
  "iedata/hBnPYcD2n32R34M262144.dat"      u ($1):(abs($4))                              w l  ls 4  lt 1 lw 0.5 notitle, \
  ""                                      u ($1):(abs($4))                              w p  ls 4       lw 3.0 notitle, \
  "iedata/BnPYcD3n16R18M4194304f128.dat"  u ($1):(($1 <= 12) ? abs($4) : 1/0)           w l  ls 8  lt 2 lw 0.5 notitle, \
  ""                                      u ($1):(($1 <= 12) ? abs($4) : 1/0)           w p  ls 8       lw 3.0 notitle, \
  "iedata/hBnPYcD4n16R18M262144.dat"      u ($1):(($1 <= 8)  ? abs($4) : 1/0)           w l  ls 12 lt 4 lw 0.5 notitle, \
  ""                                      u ($1):(($1 <= 8 && $4 > 0) ? abs($4) : 1/0)  w p  ls 12      lw 3.0 notitle, \
  ""                                      u ($1):(($1 <= 8 && $4 < 0) ? abs($4) : 1/0)  w p  ls 13      lw 3.0 notitle, \
  1e-100 lw 0 notitle





# right-top panel

set size    rw, th
set origin  lw, bh

set rmargin 1.0
set lmargin 7.0
unset ylabel

set label 101 "{/Arial-Italic D} = 5" at  13, 2e-5   rotate by 0  textcolor rgb color1a font lbfont
set label 102 "{/Arial-Italic D} = 6" at  19, 2.5e-4 rotate by 0  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 7" at  22, 5e-3   rotate by 38 textcolor rgb color3b font lbfont
set label 104 "{/Arial-Italic D} = 8" at  23, 3e-2   rotate by 45 textcolor rgb color4b font lbfont

plot [2:28][1e-5:1] \
  "data/D5/BnD5n12.dat"                   u ($1):(abs($2))                                  w l ls 2  lt 1 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w e ls 2       lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w e ls 3       lw 1.0 notitle, \
  "data/D6/BnD6n16.dat"                   u ($1):(abs($2))                                  w l ls 6  lt 2 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w e ls 6       lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w e ls 7       lw 1.0 notitle, \
  "data/D7/BnD7n20.dat"                   u ($1):(abs($2))                                  w l ls 10 lt 4 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w e ls 10      lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w e ls 11      lw 1.0 notitle, \
  "data/D8/BnD8n24.dat"                   u ($1):(abs($2))                                  w l ls 14 lt 5 lw 0.5 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w e ls 14      lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w e ls 15      lw 1.0 notitle, \
  "iedata/BnPYcD5n16R16M4194304.dat"      u ($1):(($1 <= 5) ? abs($4) : 1/0)                w l ls 4  lt 1 lw 0.5 notitle, \
  ""                                      u ($1):(($1 <= 5 && $4 > 0) ? abs($4) : 1/0)      w p ls 4       lw 1.0 notitle, \
  "iedata/hBnPYcD6n32R34M262144.dat"      u ($1):(abs($4))                                  w l ls 8  lt 2 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)                 w p ls 8       lw 1.0 notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)                 w p ls 9       lw 1.0 notitle, \
  "iedata/BnPYcD7n32R34M4194304.dat"      u ($1):(abs($4))                                  w l ls 12 lt 4 lw 0.5 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)                 w p ls 12      lw 1.0 notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)                 w p ls 13      lw 1.0 notitle, \
  "iedata/hBnPYcD8n64R66M262144.dat"      u ($1):(abs($4))                                  w l ls 16 lt 5 lw 0.5 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0)                   w p ls 16      lw 2.0 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0)                   w p ls 17      lw 2.0 notitle, \
  1e-100 lw 0 notitle

unset label 104




# left-bottom panel

set size    lw, bh
set origin 0.0, 0.0

set tmargin 0.
set bmargin 2.5
set xlabel thexlabel font lbfont offset 2, 1.0

set lmargin 6.0
set format y '10^{%T}'
set ylabel theylabel font lbfont offset 1.5, 0.0

set rmargin 0.

set label 101 "{/Arial-Italic D} = 9"   at  23.5, 2.4e-2  rotate by 0  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 10"  at  10.0, 1.0e-2  rotate by 0  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 11"  at  14.0, 2.6e-3  rotate by 0  textcolor rgb color3b font lbfont

set arrow from 25.0, 2.6e-2 to 23.2, 3.6e-2 ls 4  lt 1 nohead
set arrow from 13.0, 9.0e-3 to 16.0, 7.5e-3 ls 8  lt 2 nohead
set arrow from 14.0, 2.8e-3 to 12.1, 3.1e-3 ls 12 lt 4 nohead

plot [2:28][2e-3:4e-1] \
  "data/D9r1n20/BnD9n20.dat"              u ($1):(abs($2))                      w l ls 2  lt 1   lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3   w e ls 2  ps 1.5 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3   w e ls 3  ps 1.5 lw 1.0 notitle, \
  "data/D10r1n32/BnD10n32.dat"            u ($1):(abs($2))                      w l ls 6  lt 2   lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3   w e ls 6  ps 1.8 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3   w e ls 7  ps 1.8 lw 1.0 notitle, \
  "data/D11r1n32/BnD11n32.dat"            u ($1):(abs($2))                      w l ls 10 lt 4   lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3   w e ls 10 ps 1.5 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3   w e ls 11 ps 1.5 lw 1.0 notitle, \
  "iedata/BnPYcD9n32R34M4194304f128.dat"  u ($1):(abs($4))                      w l ls 4  lt 1   lw 0.3 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0)       w p ls 4  ps 1.8 lw 1.0 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0)       w p ls 5  ps 1.8 lw 1.0 notitle, \
  "iedata/hBnPYcD10n128R130M262144.dat"   u ($1):(abs($4))                      w l ls 8  lt 2   lw 0.3 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0)       w p ls 8  ps 1.5 lw 1.0 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0)       w p ls 9  ps 1.5 lw 1.0 notitle, \
  "iedata/BnPYcD11n32R34M4194304f128.dat" u ($1):(abs($4))                      w l ls 12 lt 4   lw 0.3 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0)       w p ls 12 ps 1.6 lw 1.0 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0)       w p ls 13 ps 1.6 lw 1.0 notitle, \
  1e-100 lw 0 notitle

unset arrow




# right-bottom panel

set size    rw, bh
set origin  lw, 0.0

set lmargin 7.0
set rmargin 1.0

set xtics 10 font tcfont offset 0, 0.5
set mxtics 10
unset ylabel
set ytics 1e-4, 10

set label 101 "{/Arial-Italic D} = 12"  at  20.0, 3.0e0   rotate by 0  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 13"  at   6.0, 2.0e-2  rotate by 0  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 14"  at  30.0, 1.0e-2  rotate by 0  textcolor rgb color3b font lbfont

set arrow from 26.0, 2.2e0  to 28.7, 3.8e-1 ls 4  lt 1 nohead
set arrow from 10.0, 1.5e-2 to 11.0, 1.5e-3 ls 8  lt 2 nohead
set arrow from 30.0, 1.2e-2 to 25.3, 3.0e-2 ls 12 lt 4 nohead

plot [2:64][1e-4:1e7] \
  "data/D12r1n64/BnD12n64.dat"            u ($1):(abs($2))                    w l ls 2  lt 1   lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 2  ps 1.0 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 3  ps 1.0 lw 1.0 notitle, \
  "data/D13r1n64/BnD13n64.dat"            u ($1):(abs($2))                    w l ls 6  lt 2   lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 6  ps 1.2 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 7  ps 1.2 lw 1.0 notitle, \
  "data/D14r1n64/BnD14n64.dat"            u ($1):(abs($2))                    w l ls 10 lt 4   lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 10 ps 1.0 lw 1.0 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 11 ps 1.0 lw 1.0 notitle, \
  "iedata/hBnPYcD12n128R130M262144.dat"   u ($1):(abs($4))                    w l ls 4  lt 1   lw 0.3 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)   w p ls 4  ps 1.2 lw 1.0 notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)   w p ls 5  ps 1.2 lw 1.0 notitle, \
  "iedata/BnPYcD13n64R66M4194304f128.dat" u ($1):(abs($4))                    w l ls 8  lt 2   lw 0.3 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0)     w p ls 8  ps 1.2 lw 1.0 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0)     w p ls 9  ps 1.2 lw 1.0 notitle, \
  "iedata/hBnPYcD14n128R130M262144.dat"   u ($1):(abs($4))                    w l ls 12 lt 4   lw 0.3 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)   w p ls 12 ps 1.1 lw 1.0 notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)   w p ls 13 ps 1.1 lw 1.0 notitle, \
  1e-100 lw 0 notitle

unset arrow



unset multiplot
unset output
set terminal wxt
reset



