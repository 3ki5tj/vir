#!/usr/bin/env gnuplot
unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript eps enhanced size 7, 9 font "Arial, 20"
set output "ievirhigh.eps"




set multiplot


# the main plot
set origin 0, 0
set size 1, 1


tcfont="Arial, 16"
lbfont = "Arial, 20"
thexlabel='Order {/Arial-Italic n}'
theylabel='{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}'


set xtics font tcfont offset 0, 0.3
set logscale y
set ytics ("1 " 1.0)
set ytics add 1e-100,1e5,1e100 font tcfont offset 0.3, 0
set format y '10^{%T}'

lbfont2 = "Arial, 12"

color1a = "#ee2222"
color1b = "#000000"

color2a = "#cc3333"
color2b = "#444444"




set tmargin 1.
set lmargin 6.
set rmargin 6.
set format y '10^{%T}'
set ylabel theylabel font lbfont offset 1.5, 2.0

set xtics 10
set mxtics 10
set xlabel thexlabel font lbfont offset 0.0, 0.8


set style line 1  lc rgb "#aaaaaa" lt 1 lw 0.2

# Mayer sampling
set style line 2  lc rgb color1a lt 1 lw 2.0 pt 4  ps 1.0 # empty square
set style line 3  lc rgb color1a lt 1 lw 2.0 pt 5  ps 1.0 # full  square

# DSC
set style line 4  lc rgb color1b lt 2 lw 2.0 pt 12 ps 1.4 # empty diamond
set style line 5  lc rgb color1b lt 2 lw 2.0 pt 13 ps 1.4 # full  diamond

# Mayer sampling
set style line 6  lc rgb color2a lt 1 lw 2.0 pt 10 ps 1.4 # empty inverted triangle
set style line 7  lc rgb color2a lt 1 lw 2.0 pt 11 ps 1.4 # full  inverted triangle

# DSC
set style line 8  lc rgb color2b lt 2 lw 2.0 pt 8  ps 1.4 # empty triangle
set style line 9  lc rgb color2b lt 2 lw 2.0 pt 9  ps 1.4 # full  triangle



xlbl = 128 + 1

D13nmax = 64
D13xmax = D13nmax - 10
D13ymax = 3e6

D14nmax = 96
D14xmax = D14nmax - 10
D14ymax = 5e13


set label 15 "{/Arial-Italic D} = 15"  at xlbl, 2.5e21 font lbfont textcolor rgb color1b
set label 16 "{/Arial-Italic D} = 16"  at xlbl, 1e21   font lbfont textcolor rgb color2b
set label 17 "{/Arial-Italic D} = 17"  at xlbl, 4.5e20 font lbfont textcolor rgb color1b
set label 18 "{/Arial-Italic D} = 18"  at xlbl, 1.8e20 font lbfont textcolor rgb color2b
set label 19 "{/Arial-Italic D} = 19"  at xlbl, 5.5e19 font lbfont textcolor rgb color1b
set label 20 "{/Arial-Italic D} = 20"  at xlbl, 1.5e19 font lbfont textcolor rgb color2b
set label 21 "{/Arial-Italic D} = 21"  at xlbl, 4e18   font lbfont textcolor rgb color1b
set label 22 "{/Arial-Italic D} = 22"  at xlbl, 9e17   font lbfont textcolor rgb color2b
set label 23 "{/Arial-Italic D} = 23"  at xlbl, 1.8e17 font lbfont textcolor rgb color1b
set label 24 "{/Arial-Italic D} = 24"  at xlbl, 3.6e16 font lbfont textcolor rgb color2b
set label 25 "{/Arial-Italic D} = 25"  at xlbl, 5.9e15 font lbfont textcolor rgb color1b
set label 26 "{/Arial-Italic D} = 26"  at xlbl, 1e15   font lbfont textcolor rgb color2b
set label 27 "{/Arial-Italic D} = 27"  at xlbl, 2e14   font lbfont textcolor rgb color1b
set label 28 "{/Arial-Italic D} = 28"  at xlbl, 3e13   font lbfont textcolor rgb color2b
set label 29 "{/Arial-Italic D} = 29"  at xlbl, 5e12   font lbfont textcolor rgb color1b
set label 30 "{/Arial-Italic D} = 30"  at xlbl, 1e12   font lbfont textcolor rgb color2b

plot [2:128][1e-9:1e22] \
  "data/D15r1n128/BnD15n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D16r1n128/BnD16n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D17r1n128/BnD17n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D18r1n128/BnD18n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D19r1n128/BnD19n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D20r1n128/BnD20n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D21r1n128/BnD21n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D22r1n128/BnD22n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D23r1n128/BnD23n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D24r1n128/BnD24n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D25r1n128/BnD25n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D26r1n128/BnD26n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D27r1n128/BnD27n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D28r1n128/BnD28n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "data/D29r1n128/BnD29n128.dat"              u ($1):(abs($2))                  w l ls 2 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3        notitle, \
  "data/D30r2n128/BnD30n128.dat"              u ($1):(abs($2))                  w l ls 6 lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6        notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7        notitle, \
  "iedata/xBnPYcD15n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD16n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD17n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD18n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD19n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD20n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD21n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD22n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD23n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD24n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD25n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD26n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD27n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD28n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  "iedata/xBnPYcD29n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 notitle, \
  "iedata/xBnPYcD30n128.dat"                  u ($1):(abs($2))                  w l ls 1 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 notitle, \
  1e-100 lw 0 notitle



# magnification of the top-right corner
set origin 0.13, 0.53
set size 0.37, 0.43

# remove unnecessary decorations
unset ylabel
unset xlabel
unset label 26
unset label 27
unset label 28
unset label 29
unset label 30
set tmargin 0
set lmargin 0
set rmargin 0
set bmargin 0

ps1 = 1.4
ps2 = 2.0
ps3 = 2.0

plot [100:128][3e14:3e21] \
  "data/D15r1n128/BnD15n128.dat"              u ($1):(abs($2))                  w l ls 2         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 ps ps1         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 ps ps1         notitle, \
  "data/D16r1n128/BnD16n128.dat"              u ($1):(abs($2))                  w l ls 6         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7 ps ps3         notitle, \
  "data/D17r1n128/BnD17n128.dat"              u ($1):(abs($2))                  w l ls 2         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 ps ps1         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 ps ps1         notitle, \
  "data/D18r1n128/BnD18n128.dat"              u ($1):(abs($2))                  w l ls 6         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7 ps ps3         notitle, \
  "data/D19r1n128/BnD19n128.dat"              u ($1):(abs($2))                  w l ls 2         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 ps ps1         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 ps ps1         notitle, \
  "data/D20r1n128/BnD20n128.dat"              u ($1):(abs($2))                  w l ls 6         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7 ps ps3         notitle, \
  "data/D21r1n128/BnD21n128.dat"              u ($1):(abs($2))                  w l ls 2         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 ps ps1         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 ps ps1         notitle, \
  "data/D22r1n128/BnD22n128.dat"              u ($1):(abs($2))                  w l ls 6         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7 ps ps3         notitle, \
  "data/D23r1n128/BnD23n128.dat"              u ($1):(abs($2))                  w l ls 2         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 ps ps1         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 ps ps1         notitle, \
  "data/D24r1n128/BnD24n128.dat"              u ($1):(abs($2))                  w l ls 6         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 6 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 7 ps ps3         notitle, \
  "data/D25r1n128/BnD25n128.dat"              u ($1):(abs($2))                  w l ls 2         lw 0.3 notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 2 ps ps1         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 3 ps ps1         notitle, \
  "iedata/xBnPYcD15n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 ps ps2         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 ps ps2         notitle, \
  "iedata/xBnPYcD16n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 ps ps3         notitle, \
  "iedata/xBnPYcD17n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 ps ps2         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 ps ps2         notitle, \
  "iedata/xBnPYcD18n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 ps ps3         notitle, \
  "iedata/xBnPYcD19n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 ps ps2         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 ps ps2         notitle, \
  "iedata/xBnPYcD20n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 ps ps3         notitle, \
  "iedata/xBnPYcD21n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 ps ps2         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 ps ps2         notitle, \
  "iedata/xBnPYcD22n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 ps ps2         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 ps ps2         notitle, \
  "iedata/xBnPYcD23n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 ps ps3         notitle, \
  "iedata/xBnPYcD24n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 8 ps ps2         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 9 ps ps2         notitle, \
  "iedata/xBnPYcD25n128.dat"                  u ($1):(abs($2))                  w l ls 1                notitle, \
  ""                                          u ($1):(($2 > 0) ? abs($2) : 1/0) w p ls 4 ps ps3         notitle, \
  ""                                          u ($1):(($2 < 0) ? abs($2) : 1/0) w p ls 5 ps ps3         notitle, \
  1e-100 lw 0 notitle



unset multiplot

unset output
set terminal wxt
reset



