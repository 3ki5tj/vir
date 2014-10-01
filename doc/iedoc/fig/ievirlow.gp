#!/usr/bin/env gnuplot
unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript eps enhanced size 7, 7 font "Arial, 20"
set output "ievirlow.eps"

tcfont="Arial, 16"
thexlabel='Order {/Arial-Italic n}'
theylabel='{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}'

# height of the bottom panels
bh = 0.4
# height of the top panels
th = 1 - bh

# width of the right panel
rw = 0.52
# width of the left panel
lw = 1 - rw

set logscale y
set ytics ("1 " 1.0, "10" 10.0)
set ytics add 1e-9,10,1e9 font tcfont offset 0.3, 0
set mytics 10
set format y '10^{%T}'

lbfont  = "Arial, 18"
lbfontsm  = "Arial, 12"

color1a   = "#cc3333"
color1b   = "#222222"
color1b1  = "#448855"
color1b2  = "#555588"

color2a   = color1a   # "#000000"
color2b   = color1b   # "#804000"
color2b1  = color1b1  # "#c07030"
color2b2  = color1b2  # "#dd9900"

color3a   = color1a   # "#600080"
color3b   = color1b   # "#006000"
color3b1  = color1b1
color3b2  = color1b2

color4a   = color1a   # "#a0a0a0"
color4b   = color1b   # "#008080"
color4b1  = color1b1
color4b2  = color1b2

colorarr  = "#777777"

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 1 lw 1

# Mayer sampling
set style line 2  lc rgb color1a  lt 1 lw 1.5 pt 4  ps 1.4 # empty square
set style line 3  lc rgb color1a  lt 1 lw 1.5 pt 5  ps 1.4 # full  square

# Detailed self-consistent (DSC)
set style line 4  lc rgb color1b  lt 2 lw 1.5 pt 12 ps 2.0 # empty diamond
set style line 5  lc rgb color1b  lt 2 lw 1.5 pt 13 ps 2.0 # full  diamond

# DSC_lambda
set style line 24 lc rgb color1b1 lt 5 lw 1.5 pt 12 ps 2.0 # empty diamond
set style line 25 lc rgb color1b1 lt 5 lw 1.5 pt 13 ps 2.0 # full  diamond

# modified DSC
set style line 34 lc rgb color1b2 lt 4 lw 1.5 pt 12 ps 2.0 # empty diamond
set style line 35 lc rgb color1b2 lt 4 lw 1.5 pt 13 ps 2.0 # full  diamond



# Mayer sampling
set style line 6  lc rgb color2a  lt 1 lw 1.5 pt 10 ps 2.0 # empty inverted triangle
set style line 7  lc rgb color2a  lt 1 lw 1.5 pt 11 ps 2.0 # full  inverted triangle

# DSC
set style line 8  lc rgb color2b  lt 2 lw 1.5 pt 8  ps 2.0 # empty triangle
set style line 9  lc rgb color2b  lt 2 lw 1.5 pt 9  ps 2.0 # full  triangle

# DSC_lambda
set style line 28 lc rgb color2b1 lt 5 lw 1.5 pt 8  ps 2.0 # empty triangle
set style line 29 lc rgb color2b1 lt 5 lw 1.5 pt 9  ps 2.0 # full  triangle

# modified DSC
set style line 38 lc rgb color2b2 lt 4 lw 1.5 pt 8  ps 2.0 # empty triangle
set style line 39 lc rgb color2b2 lt 4 lw 1.5 pt 9  ps 2.0 # full  triangle



# Mayer sampling
set style line 10 lc rgb color3a  lt 1 lw 1.5 pt 6  ps 1.4 # empty circle
set style line 11 lc rgb color3a  lt 1 lw 1.5 pt 7  ps 1.4 # full  circle

# DSC
set style line 12 lc rgb color3b  lt 2 lw 1.5 pt 14 ps 1.8 # empty pentagon
set style line 13 lc rgb color3b  lt 2 lw 1.5 pt 15 ps 1.8 # full  pentagon

# DSC_lambda
set style line 22 lc rgb color2b1 lt 5 lw 1.5 pt 14 ps 1.8 # empty pentagon
set style line 23 lc rgb color2b1 lt 5 lw 1.5 pt 15 ps 1.8 # full  pentagon



# Mayer sampling
set style line 14 lc rgb color4a  lt 1 lw 2.0 pt 1  ps 2.0 # add +, thin
set style line 15 lc rgb color4a  lt 1 lw 7.0 pt 1  ps 2.0 # add +, thick

# DSC
set style line 16 lc rgb color4b  lt 2 lw 2.0 pt 2  ps 1.2 # cross X, thin
set style line 17 lc rgb color4b  lt 2 lw 9.0 pt 2  ps 1.2 # cross X, thick

# DSC_lambda
set style line 26 lc rgb color4b1 lt 5 lw 2.0 pt 2  ps 1.2 # cross X, thin
set style line 27 lc rgb color4b1 lt 5 lw 9.0 pt 2  ps 1.2 # cross X, thick



tagdx1 = 0.005
tagdx2 = 0.015
tagdy1 = 0.020
tagdy2 = 0.005
tagfont = "Arial, 24"
set label 1000 "(a)" at screen       tagdx1,  1 - tagdy1 font tagfont
set label 1001 "(b)" at screen  lw + tagdx2,  1 - tagdy1 font tagfont
set label 1002 "(c)" at screen       tagdx1, bh - tagdy2 font tagfont
set label 1003 "(d)" at screen  lw + tagdx2, bh - tagdy2 font tagfont



set multiplot




# left-top panel

set size    lw, th
set origin 0.0, bh

set xtics 4 font tcfont offset 0, 0.5
set mxtics 4
unset xlabel

set ylabel theylabel offset 1.3, 0.0

set tmargin 1.
set bmargin 1.5
set rmargin 0.
set lmargin 6.0

set label 101 "{/Arial-Italic D} = 2, {/Symbol-Oblique k} = 0" \
  at 16.0, 1.3e-4 rotate by -53  textcolor rgb color1b font lbfont

set label 301 "{/Arial-Italic D} = 2, {/Symbol-Oblique k} = (0.194)_{{/Arial-Italic n} {/Symbol \263} 6}" \
  at 14.5, 3.0e-3 rotate by -52  textcolor rgb color1b2 font lbfont

set label 102 "{/Arial-Italic D} = 3, {/Symbol-Oblique k} = 0" \
  at  9.5, 2.4e-4 rotate by -72  textcolor rgb color2b font lbfont

set label 302 "{/Arial-Italic D} = 3, {/Symbol-Oblique k} = (0.335)_{{/Arial-Italic n} {/Symbol \263} 5}" \
  at 10.5, 8.5e-4 rotate by -68  textcolor rgb color2b2 font lbfont

set label 103 "{/Arial-Italic D} = 4" \
  at  5.0, 5.5e-3 rotate by -78  textcolor rgb color3b font lbfont

#set label 203 "{/Arial-Italic D} = 4, DSC_{/Symbol-Oblique l}" \
#  at  7.5, 0.5e-4 rotate by -85  textcolor rgb color3b1 font lbfont

plot [2:22][8e-6:2] \
  "data/D2/BnD2n14.dat"                   u ($1):(abs($2)):3                            w l  ls 2              notitle, \
  ""                                      u ($1):(abs($2)):3                            w p  ls 2              notitle, \
  "data/D3/BnD3n12.dat"                   u ($1):(abs($2)):3                            w l  ls 6              notitle, \
  ""                                      u ($1):(abs($2)):3                            w p  ls 6              notitle, \
  "data/D4/BnD4n11.dat"                   u ($1):(($1 <= 10) ? abs($2) : 1/0):3         w l  ls 10             notitle, \
  ""                                      u ($1):(($1 <= 10) ? abs($2) : 1/0):3         w p  ls 10             notitle, \
  "iedata/xBnPYcD2n32.dat"                u ($1):(abs($2))                              w l  ls 4              notitle, \
  ""                                      u ($1):(abs($2))                              w p  ls 4              notitle, \
  "iedata/kappa/hBnPYcD2n36R38M32768c0.194L4.dat"  \
                                          u ($1):(abs($4))                              w l  ls 34             notitle, \
  ""                                      u ($1):(abs($4))                              w p  ls 34             notitle, \
  "iedata/xBnPYcD3n16.dat"                u ($1):(($1 <= 12) ? abs($2) : 1/0)           w l  ls 8              notitle, \
  ""                                      u ($1):(($1 <= 12) ? abs($2) : 1/0)           w p  ls 8              notitle, \
  "iedata/kappa/BnPYcD3n36R38M32768c0.335L3f128.dat"  \
                                          u ($1):(($1 <= 12) ? abs($2) : 1/0)           w l  ls 38             notitle, \
  ""                                      u ($1):(($1 <= 12) ? abs($2) : 1/0)           w p  ls 38             notitle, \
  "iedata/xBnPYcD4n16.dat"                u ($1):(($1 <= 8)  ? abs($2) : 1/0)           w l  ls 12             notitle, \
  ""                                      u ($1):(($1 <= 8 && $2 > 0) ? abs($2) : 1/0)  w p  ls 12             notitle, \
  ""                                      u ($1):(($1 <= 8 && $2 < 0) ? abs($2) : 1/0)  w p  ls 13             notitle, \
  1e-100 lw 0 notitle



#  "iedata/lamc/xBnPYlD4n12.dat"           u ($1):(($1 <= 10) ? abs($2) : 1/0)           w l  ls 22             notitle, \
#  ""                                      u ($1):(($1 <= 10 && $2 > 0) ? abs($2) : 1/0) w p  ls 22             notitle, \
#  ""                                      u ($1):(($1 <= 10 && $2 < 0) ? abs($2) : 1/0) w p  ls 23             notitle, \



unset label





# right-top panel

set size    rw, th
set origin  lw, bh

set rmargin 1.5
set lmargin 5.0
unset ylabel

set label 101 "{/Arial-Italic D} = 5, DSC_{/Symbol-Oblique l}" \
    at   9.5, 2.5e-4   rotate by -55  textcolor rgb color1b1 font lbfont

set label 102 "{/Arial-Italic D} = 6, {/Symbol-Oblique k} = 0" \
    at  20, 5.0e-4 rotate by 8  textcolor rgb color2b font lbfont

set label 202 "{/Arial-Italic D} = 6, {/Symbol-Oblique k} = (0.352)_{{/Arial-Italic n} {/Symbol \263} 4}" \
    at  24.0, 2.5e-4 rotate by 10  textcolor rgb color2b2 font lbfontsm

set label 302 "{/Arial-Italic D} = 6, DSC_{/Symbol-Oblique l}" \
    at  24.5, 1.2e-4 rotate by -20 textcolor rgb color2b1 font lbfont

set label 103 "{/Arial-Italic D} = 7" at  22, 5e-3   rotate by 30 textcolor rgb color3b font lbfont
set label 104 "{/Arial-Italic D} = 8" at  23, 3e-2   rotate by 36 textcolor rgb color4b font lbfont



plot [2:32][1e-5:2] \
  "data/D5/BnD5n12.dat"                   u ($1):(abs($2))                                  w l ls 2              notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w p ls 2              notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w p ls 3              notitle, \
  "data/D6/BnD6n16.dat"                   u ($1):(abs($2))                                  w l ls 6              notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w p ls 6              notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w p ls 7              notitle, \
  "data/D7/BnD7n20.dat"                   u ($1):(abs($2))                                  w l ls 10             notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w p ls 10             notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w p ls 11             notitle, \
  "data/D8/BnD8n24.dat"                   u ($1):(abs($2))                                  w l ls 14             notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3               w p ls 14             notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3               w p ls 15             notitle, \
  "iedata/xBnPYcD5n16.dat"                u ($1):(($1 <= 5) ? abs($2) : 1/0)                w l ls 4              notitle, \
  ""                                      u ($1):(($1 <= 5 && $2 > 0) ? abs($2) : 1/0)      w p ls 4              notitle, \
  "iedata/lamc/xBnPYlD5n16.dat"           u ($1):(($1 <= 12) ? abs($2) : 1/0)               w l ls 24             notitle, \
  ""                                      u ($1):(($1 <= 12 && $2 > 0) ? abs($2) : 1/0)     w p ls 24             notitle, \
  ""                                      u ($1):(($1 <= 12 && $2 < 0) ? abs($2) : 1/0)     w p ls 25             notitle, \
  "iedata/xBnPYcD6n128.dat"               u ($1):(abs($2))                                  w l ls 8              notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)                 w p ls 8              notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)                 w p ls 9              notitle, \
  "iedata/kappa/hBnPYcD6n36R38M32768c0.352L2.dat"  \
                                          u ($1):(abs($4))                                  w l ls 38             notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0)                 w p ls 38             notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0)                 w p ls 39             notitle, \
  "iedata/lamc/xBnPYlD6n32.dat"           u ($1):(abs($2))                                  w l ls 28             notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)                 w p ls 28             notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)                 w p ls 29             notitle, \
  "iedata/xBnPYcD7n128.dat"               u ($1):(abs($2))                                  w l ls 12             notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)                 w p ls 12             notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)                 w p ls 13             notitle, \
  "iedata/xBnPYcD8n128.dat"               u ($1):(abs($2))                                  w l ls 16             notitle, \
  ""                                      u ($1):($2 > 0 ? abs($2) : 1/0)                   w p ls 16             notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0)                   w p ls 17             notitle, \
  1e-100 lw 0 notitle





#  "iedata/lamc/xBnPYlD7n36.dat"           u ($1):(abs($2))                                  w l ls 22             notitle, \
#  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)                 w p ls 22             notitle, \
#  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)                 w p ls 23             notitle, \

#  "iedata/lamc/xBnPYlD8n36.dat"           u ($1):(abs($2))                                  w l ls 26             notitle, \
#  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)                 w p ls 26             notitle, \
#  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)                 w p ls 27             notitle, \



unset label




# left-bottom panel

set size    lw, bh
set origin 0.0, 0.0

set tmargin 0.
set bmargin 2.5
set xlabel thexlabel font lbfont offset 2, 1.0

set lmargin 6.0
set format y '10^{%T}'
set ylabel theylabel font lbfont offset 1.5, 0.7

set rmargin 0.

set label 101 "{/Arial-Italic D} = 9"   at  19.5, 6.0e-3  rotate by 0  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 10"  at  11.0, 1.5e-2  rotate by 0  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 11"  at  14.0, 2.6e-3  rotate by 0  textcolor rgb color3b font lbfont

set arrow from 20.5, 6.8e-3 to 19.0, 1.2e-2 filled ls 4  lt 1 lc rgb colorarr head
set arrow from 13.0, 1.3e-2 to 16.0, 7.7e-3 filled ls 8  lt 1 lc rgb colorarr head
set arrow from 14.0, 2.8e-3 to 12.1, 3.2e-3 filled ls 12 lt 1 lc rgb colorarr head

plot [2:24][2e-3:3e-1] \
  "data/D9r1n20/BnD9n20.dat"              u ($1):(abs($2))                      w l ls 2         lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3   w p ls 2  ps 1.0        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3   w p ls 3  ps 1.0        notitle, \
  "data/D10r1n32/BnD10n32.dat"            u ($1):(abs($2))                      w l ls 6         lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3   w p ls 6  ps 1.2        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3   w p ls 7  ps 1.2        notitle, \
  "data/D11r1n32/BnD11n32.dat"            u ($1):(abs($2))                      w l ls 10        lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0):3   w p ls 10 ps 1.0        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3   w p ls 11 ps 1.0        notitle, \
  "iedata/xBnPYcD9n128.dat"               u ($1):(abs($2))                      w l ls 4         lw 0.3 notitle, \
  ""                                      u ($1):($2 > 0 ? abs($2) : 1/0)       w p ls 4  ps 1.4        notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0)       w p ls 5  ps 1.4        notitle, \
  "iedata/xBnPYcD10n128.dat"              u ($1):(abs($2))                      w l ls 8         lw 0.3 notitle, \
  ""                                      u ($1):($2 > 0 ? abs($2) : 1/0)       w p ls 8  ps 1.2        notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0)       w p ls 9  ps 1.2        notitle, \
  "iedata/xBnPYcD11n128.dat"              u ($1):(abs($2))                      w l ls 12        lw 0.3 notitle, \
  ""                                      u ($1):($2 > 0 ? abs($2) : 1/0)       w p ls 12 ps 1.0        notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0)       w p ls 13 ps 1.0        notitle, \
  1e-100 lw 0 notitle

unset label
unset arrow




# right-bottom panel

set size    rw, bh
set origin  lw, 0.0

set lmargin 5.0
set rmargin 1.5

set xtics 10 font tcfont offset 0, 0.5
set mxtics 10
unset ylabel
#set ytics 1e-4, 10

set label 101 "{/Arial-Italic D} = 12"  at  23.0, 3.0e0   rotate by 0  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 13"  at   6.0, 2.0e-2  rotate by 0  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 14"  at  30.0, 1.0e-2  rotate by 0  textcolor rgb color3b font lbfont

set arrow from 26.0, 2.2e0  to 29.0, 3.9e-1 filled ls 4  lt 1 lc rgb colorarr head  # D = 12
set arrow from  8.0, 1.5e-2 to 11.0, 1.7e-3 filled ls 8  lt 1 lc rgb colorarr head  # D = 13
set arrow from 30.0, 1.2e-2 to 25.4, 3.0e-2 filled ls 12 lt 1 lc rgb colorarr head  # D = 14

plot [2:40][5e-4:5e1] \
  "data/D12r1n64/BnD12n64.dat"            u ($1):(abs($2))                    w l ls 2         lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 2  ps 1.0        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 3  ps 1.0        notitle, \
  "data/D13r1n64/BnD13n64.dat"            u ($1):(abs($2))                    w l ls 6         lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 6  ps 1.5        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 7  ps 1.5        notitle, \
  "data/D14r1n64/BnD14n64.dat"            u ($1):(abs($2))                    w l ls 10        lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 10 ps 1.2        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 11 ps 1.2        notitle, \
  "iedata/xBnPYcD12n128.dat"              u ($1):(abs($2))                    w l ls 4         lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 4  ps 1.4        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 5  ps 1.4        notitle, \
  "iedata/xBnPYcD13n128.dat"              u ($1):(abs($2))                    w l ls 8         lw 0.3 notitle, \
  ""                                      u ($1):($2 > 0 ? abs($2) : 1/0)     w p ls 8  ps 1.5        notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0)     w p ls 9  ps 1.5        notitle, \
  "iedata/xBnPYcD14n128.dat"              u ($1):(abs($2))                    w l ls 12        lw 0.3 notitle, \
  ""                                      u ($1):(($2 > 0) ? abs($2) : 1/0)   w p ls 12 ps 1.2        notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0)   w p ls 13 ps 1.2        notitle, \
  1e-100 lw 0 notitle

unset arrow
unset label



unset multiplot
unset output
set terminal wxt
reset



