unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript eps enhanced size 10, 7 font "Arial, 18"
set output "ievir.eps"

tcfont="Arial, 14"
thexlabel='Order {/Arial-Italic n}'
theylabel='{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}'

# height of the bottom panels
bh = 0.55
# height of the top panels
th = 1 - bh

# width of the central panel
cw = 0.225
# width of the left panel
lw = 0.26

rx = lw + cw
# width of the right panel
rw = 1 - lw - cw

set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
unset xlabel

set logscale y
set ytics font tcfont offset 0.3, 0
set format y '10^{%T}'
set ylabel theylabel offset 2.3, 0.7

lbfont  = "Arial, 16"
lbfont2 = "Arial, 12"


set multiplot




# left-top panel

set size    lw, th
set origin 0.0, bh

set tmargin 1.
set bmargin 1.5
set rmargin 0.
set lmargin 5.

color1a = "#cc2020"
color1b = "#2050aa"

color2a = "#cc2080"
color2b = "#208020"

color3a = "#aa7720"
color3b = "#40aaaa"

color4a = "#804020"
color4b = "#808080"

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 3 lw 1

set style line 2  lc rgb color1a lt 1 lw 2 pt 6  ps 1.0
set style line 3  lc rgb color1a lt 1 lw 2 pt 7  ps 1.0

set style line 4  lc rgb color1b lt 1 lw 2 pt 4  ps 1.5
set style line 5  lc rgb color1b lt 1 lw 2 pt 5  ps 1.5

set style line 6  lc rgb color2a lt 1 lw 2 pt 6  ps 1.0
set style line 7  lc rgb color2a lt 1 lw 2 pt 7  ps 1.0

set style line 8  lc rgb color2b lt 1 lw 2 pt 8  ps 1.5
set style line 9  lc rgb color2b lt 1 lw 2 pt 9  ps 1.5

set style line 10 lc rgb color3a lt 1 lw 2 pt 6  ps 1.0
set style line 11 lc rgb color3a lt 1 lw 2 pt 7  ps 1.0

set style line 12 lc rgb color3b lt 1 lw 2 pt 10 ps 1.5
set style line 13 lc rgb color3b lt 1 lw 2 pt 11 ps 1.5

set style line 14 lc rgb color4a lt 1 lw 2 pt 6  ps 1.0
set style line 15 lc rgb color4a lt 1 lw 2 pt 7  ps 1.0

set style line 16 lc rgb color4b lt 1 lw 2 pt 12 ps 1.5
set style line 17 lc rgb color4b lt 1 lw 2 pt 13 ps 1.5




set label 101 "{/Arial-Italic D} = 2" at 18.5, 3e-4  rotate by -65  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 3" at 13, 1.5e-4  rotate by -75  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 4" at  5, 4e-3    rotate by -80  textcolor rgb color3b font lbfont

plot [2:24][8e-6:] \
  "data/D2/BnD2n14.dat"                   u ($1):(abs($2)):3 w errorbars ls 2  notitle, \
  "data/D3/BnD3n12.dat"                   u ($1):(abs($2)):3 w errorbars ls 6  notitle, \
  "data/D4/BnD4n11.dat"                   u ($1):(abs($2)):3 w errorbars ls 10 notitle, \
  "iedata/hBnPYcD2n32R34M32768.dat"       u ($1):(abs($4))                    w l ls 1 notitle, \
  ""                                      u ($1):(abs($4))                    w p ls 4 notitle, \
  "iedata/BnPYcD3n16R18M1048576.dat"      u ($1):(($1 <= 12) ? abs($4) : 1/0)   w l ls 1 notitle, \
  ""                                      u ($1):(($1 <= 12) ? abs($4) : 1/0) w p ls 8 notitle, \
  "iedata/hBnPYcD4n10R12M32768.dat"       u ($1):(($1 <= 8) ? abs($4) : 1/0)    w l ls 1 notitle, \
  ""                                      u ($1):(($1 <= 8 && $4 > 0) ? abs($4) : 1/0)  w p ls 12 notitle, \
  ""                                      u ($1):(($1 <= 8 && $4 < 0) ? abs($4) : 1/0)  w p ls 13 notitle, \
  1e-100 lw 0 notitle





# central-top panel

set size    cw, th
set origin  lw, bh

set lmargin 0
unset ylabel
set format y ''

set label 101 "{/Arial-Italic D} = 5" at  13, 2e-5   rotate by 0  textcolor rgb color1a font lbfont
set label 102 "{/Arial-Italic D} = 6" at  17, 2.5e-4 rotate by 0  textcolor rgb color2a font lbfont
set label 103 "{/Arial-Italic D} = 7" at  23, 6e-3   rotate by 45 textcolor rgb color3b font lbfont
set label 104 "{/Arial-Italic D} = 8" at  24, 4e-2   rotate by 50 textcolor rgb color4b font lbfont

plot [2:32][8e-6:] \
  "data/D5/BnD5n12.dat"                   u ($1):($2 > 0 ? abs($2) : 1/0):3 w errorbars ls 2  notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0):3 w errorbars ls 3  notitle, \
  "data/D6/BnD6n16.dat"                   u ($1):($2 > 0 ? abs($2) : 1/0):3 w errorbars ls 6  notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0):3 w errorbars ls 7  notitle, \
  "data/D7/BnD7n20.dat"                   u ($1):($2 > 0 ? abs($2) : 1/0):3 w errorbars ls 10 notitle, \
  ""                                      u ($1):($2 < 0 ? abs($2) : 1/0):3 w errorbars ls 11 notitle, \
  "data/D8/BnD8n24.dat"                   u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 14 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 15 notitle, \
  "iedata/BnPYcD5n16R16M4194304.dat"      u ($1):(($1 <= 5) ? abs($4) : 1/0)            w l ls 1 notitle, \
  ""                                      u ($1):(($1 <= 5 && $4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  "iedata/hBnPYcD6n30R33M32768.dat"       u ($1):(($1 <= 4) ? abs($4) : 1/0)   w l ls 1 notitle, \
  ""                                      u ($1):(($1 <= 4 && $4 > 0) ? abs($4) : 1/0) w p ls 8 notitle, \
  "iedata/BnPYcD7n32R34M4194304.dat"      u ($1):(($1 <= 1000) ? abs($4) : 1/0)    w l ls 1 notitle, \
  ""                                      u ($1):(($1 <= 1000 && $4 > 0) ? abs($4) : 1/0)  w p ls 12 notitle, \
  ""                                      u ($1):(($1 <= 1000 && $4 < 0) ? abs($4) : 1/0)  w p ls 13 notitle, \
  "iedata/hBnPYcD8n128R130M32768.dat"     u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0) w p ls 16 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0) w p ls 17 notitle, \
  1e-100 lw 0 notitle

unset label 104




# left-bottom panel

set size    lw, bh
set origin 0.0, 0.0

set tmargin 0.
set bmargin 2.5
set xlabel thexlabel offset 2, 1.0

set lmargin 5.
set format y '10^{%T}'
set ylabel theylabel offset 2.3, 4.0

# line styles for the small panels
set style line 1  lc rgb "#aaaaaa" lt 3 lw 0.2

set style line 2  lc rgb color1a lt 1 lw 1 pt 6  ps 0.7
set style line 3  lc rgb color1a lt 1 lw 1 pt 7  ps 0.7

set style line 4  lc rgb color1b lt 1 lw 1 pt 4  ps 0.7
set style line 5  lc rgb color1b lt 1 lw 1 pt 5  ps 0.7

set style line 6  lc rgb color2a lt 1 lw 1 pt 6  ps 0.8
set style line 7  lc rgb color2a lt 1 lw 1 pt 7  ps 0.8

set style line 8  lc rgb color2b lt 1 lw 1 pt 8  ps 0.8
set style line 9  lc rgb color2b lt 1 lw 1 pt 9  ps 0.8

set style line 10 lc rgb color3a lt 1 lw 1 pt 6  ps 0.8
set style line 11 lc rgb color3a lt 1 lw 1 pt 7  ps 0.8

set style line 12 lc rgb color3b lt 1 lw 1 pt 10 ps 0.8
set style line 13 lc rgb color3b lt 1 lw 1 pt 11 ps 0.8


set label 101 "{/Arial-Italic D} = 9"   at  52, 2e2   rotate by 56  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 10"  at  59, 2.8e4 rotate by 60  textcolor rgb color2b font lbfont2
set label 103 "{/Arial-Italic D} = 11"  at  56, 1.5e5 rotate by 65  textcolor rgb color3b font lbfont

plot [2:64][5e-4:1e7] \
  "data/D9r1n20/BnD9n20.dat"              u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 2 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 3 notitle, \
  "data/D10r1n32/BnD10n32.dat"            u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 6 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 7 notitle, \
  "data/D11r1n32/BnD11n32.dat"            u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 10 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 11 notitle, \
  "iedata/BnPYcD9n64R66M65536.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD10n128R130M32768.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0) w p ls 8 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0) w p ls 9 notitle, \
  "iedata/BnPYcD11n64R66M2097152.dat"     u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0) w p ls 12 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0) w p ls 13 notitle, \
  1e-100 lw 0 notitle





# central-bottom panel

set size    cw, bh
set origin  lw, 0.0

set lmargin 0
set format y ''
unset ylabel

set label 101 "{/Arial-Italic D} = 12"  at  20, 4e0   rotate by 0  textcolor rgb color1b font lbfont
set label 102 "{/Arial-Italic D} = 13"  at  6, 3e-2   rotate by 0  textcolor rgb color2b font lbfont
set label 103 "{/Arial-Italic D} = 14"  at  30, 1e-2  rotate by 0  textcolor rgb color3b font lbfont

set arrow from 24, 3 to 28.7, 0.38 ls 4 nohead
set arrow from 9, 2e-2 to 11.0, 1.5e-3 ls 8 nohead
set arrow from 30, 1.2e-2 to 25.3, 3e-2 ls 12 nohead

plot [2:64][5e-4:1e7] \
  "data/D12r1n64/BnD12n64.dat"            u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 2 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 3 notitle, \
  "data/D13r1n64/BnD13n64.dat"            u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 6 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 7 notitle, \
  "data/D14r1n64/BnD14n64.dat"            u ($1):(($2 > 0) ? abs($2) : 1/0):3 w errorbars ls 10 notitle, \
  ""                                      u ($1):(($2 < 0) ? abs($2) : 1/0):3 w errorbars ls 11 notitle, \
  "iedata/hBnPYcD12n128R130M32768.dat"    u ($1):(abs($4))  w l ls 1 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD13n128R131M1048576.dat"   u ($1):(abs($4))                w l ls 1 notitle, \
  ""                                      u ($1):($4 > 0 ? abs($4) : 1/0) w p ls 8 notitle, \
  ""                                      u ($1):($4 < 0 ? abs($4) : 1/0) w p ls 9 notitle, \
  "iedata/hBnPYcD14n128R129M32768.dat"    u ($1):(abs($4))  w l ls 1 notitle, \
  ""                                      u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 12 notitle, \
  ""                                      u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 13 notitle, \
  1e-100 lw 0 notitle

unset arrow



# right panel




set size    rw, 1
set origin  rx, 0

set tmargin 1.
set lmargin 6.
set rmargin 5.
set format y '10^{%T}'
set ylabel theylabel offset 2.5, 2.0


set style line 1 lc rgb "#aaaaaa" lt 1 lw 0.5

set style line 2 lc rgb color1a lt 1 lw 1 pt 6 ps 0.7
set style line 3 lc rgb color1a lt 1 lw 1 pt 7 ps 0.7

set style line 4 lc rgb color1b lt 1 lw 1 pt 4 ps 0.6
set style line 5 lc rgb color1b lt 1 lw 1 pt 5 ps 0.6

unset label

lbfont = "Arial, 14"
xlbl = 128 + 1

D13nmax = 64
D13xmax = D13nmax - 10
D13ymax = 3e6

D14nmax = 96
D14xmax = D14nmax - 10
D14ymax = 5e13


set label "{/Arial-Italic D} = 15"  at xlbl, 3e21   font lbfont
set label "{/Arial-Italic D} = 16"  at xlbl, 1e21   font lbfont
set label "{/Arial-Italic D} = 17"  at xlbl, 4e20   font lbfont
set label "{/Arial-Italic D} = 18"  at xlbl, 1.5e20 font lbfont
set label "{/Arial-Italic D} = 19"  at xlbl, 5e19   font lbfont
set label "{/Arial-Italic D} = 20"  at xlbl, 1.5e19 font lbfont
set label "{/Arial-Italic D} = 21"  at xlbl, 4e18   font lbfont
set label "{/Arial-Italic D} = 22"  at xlbl, 1e18   font lbfont
set label "{/Arial-Italic D} = 23"  at xlbl, 1.5e17 font lbfont
set label "{/Arial-Italic D} = 24"  at xlbl, 3e16   font lbfont
set label "{/Arial-Italic D} = 25"  at xlbl, 5e15   font lbfont
set label "{/Arial-Italic D} = 26"  at xlbl, 1e15   font lbfont
set label "{/Arial-Italic D} = 27"  at xlbl, 2e14   font lbfont
set label "{/Arial-Italic D} = 28"  at xlbl, 3e13   font lbfont
set label "{/Arial-Italic D} = 29"  at xlbl, 5e12   font lbfont
set label "{/Arial-Italic D} = 30"  at xlbl, 1e12   font lbfont

plot [2:128][1e-9:1e22] \
  "data/D15r1n128/BnD15n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D16r1n128/BnD16n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D17r1n128/BnD17n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D18r1n128/BnD18n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D19r1n128/BnD19n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D20r1n128/BnD20n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D25r1n128/BnD25n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "data/D30r2n128/BnD30n128.dat"              u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 notitle, \
  "iedata/BnPYcD15n128R131M262144p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD16n128R129M32768.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD17n128R131M262144p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD18n128R129M32768.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD19n128R131M262144p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD20n128R130M65536.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD21n128R131M262144p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD22n128R130M32768.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD23n128R131M262144p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD24n128R130M32768.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD25n128R131M32768p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD26n128R130M32768.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD27n128R131M32768p256.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD28n128R130M32768.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/BnPYcD29n128R131M32768p512.dat"    u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  "iedata/hBnPYcD30n128R130M65536.dat"        u ($1):(abs($4))                  w l ls 1 notitle, \
  ""                                          u ($1):(($4 > 0) ? abs($4) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($4 < 0) ? abs($4) : 1/0) w p ls 5 notitle, \
  1e-100 lw 0 notitle


unset multiplot
unset output
set terminal wxt
reset



