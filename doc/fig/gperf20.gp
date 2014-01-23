# performance of the grand ensemble technique

unset multiplot
reset

set encoding cp1250 # make minus sign longer
#set encoding iso_8859_1
set terminal postscript landscape enhanced font "Arial, 24" size 16, 11
set output "gperf20.ps"



set multiplot

pw = 0.5  # width of each panel
ph = 0.5  # height of each panel
px = -0.22  # x-position of the label
py = 0.98   # y-position of the label
iw = 0.31  # width of the inset
ih = 0.28  # height of the inset
ix = 0.04  # x-coordinate of the inset (relative to the panel)
iy = 0.20  # y-coordinate of the inset (relative to the panel)
ilx = 0.08  # inset label x
ily = 0.84  # inset label y
iwr = 0.30
ihr = 0.26
ixr = 0.03
iyr = 0.22


fnG2  = "data/benchD20n64/refine.data"
fnG3  = "data/benchD20r1n64/ZrD20r1n64.data"
fnG3r = "data/benchD20r1n64/ZrrD20r1n64.data"

tagfont="Arial, 32"
tcfont="Arial, 18"
smtcfont="Arial, 14"
set key spacing 1.5

set style line 1 lc rgb "#aaaaaa" lt 1 lw 0.5
set style line 2 lc rgb "#000000" lt 1 lw 1 pt 6 ps 1.2
set style line 3 lc rgb "#000000" lt 1 lw 1 pt 7 ps 1.2
set style line 4 lc rgb "#000000" lt 1 lw 1

set style line 5 lc rgb "#ff8080" lt 2 lw 0.5
set style line 6 lc rgb "#ff0000" lt 2 lw 1 pt 8 ps 1.2
set style line 7 lc rgb "#ff0000" lt 2 lw 1 pt 9 ps 1.2
set style line 8 lc rgb "#ff0000" lt 2 lw 1



# top left panel
set size pw, ph
set origin 0, ph

set label 10 "(a)" at graph px, py font tagfont

set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
set xlabel "Order {/Arial-Italic n}" offset 2, 0.5

set logscale y
set mytics 5
set ytics 1e-5, 10, 1e7 font tcfont offset 0., 0
set format y '10^{%T}'
set ylabel '{/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1}' offset 0.0, 0

set key at graph 1.06, 0.02 right bottom Left reverse

#set title "Algorithm G3 ({/Arial-Italic M} = 1) versus algorithm G2, {/Arial-Italic D} = 20"
plot [2:64][1e-5:1e4] \
  fnG3      u 1:(abs($19))                        w l ls 1 notitle, \
  ""        u ($1):(($19 > 0) ? abs($19) : 1/0)   w p ls 2 notitle, \
  ""        u ($1):(($19 < 0) ? abs($19) : 1/0)   w p ls 3 notitle, \
  fnG2      u 1:(abs($23))                        w l ls 5 notitle, \
  ""        u ($1):(($23 > 0) ? abs($23) : 1/0)   w p ls 6 notitle, \
  ""        u ($1):(($23 < 0) ? abs($23) : 1/0)   w p ls 7 notitle, \
  0 w lp ls 2 t "algorithm G3", \
  0 w lp ls 6 t "algorithm G2", \
  0 w l notitle


# inset of the top left panel
set size iw, ih
set origin ix, ph + iy

unset label 10

set mxtics 5
set xtics 10 font smtcfont offset 0, 0.5
unset xlabel

unset logscale y
set mytics 5
set ytics 0.001 font smtcfont offset 0.5, 0
set format y "%g"
unset ylabel

set label 3 "{/Symbolic-Oblique e} ({/Arial-Italic B_n} /{/Arial-Italic B}_2^{{/Arial-Italic n}-1})" \
  at graph ilx, ily font tcfont

unset title
plot [2:64][:0.0034] \
  fnG3   u 1:($20/abs($19))                    w lp ls 2 ps 0.6 notitle, \
  fnG2   u 1:($24/abs($23))                    w lp ls 6 ps 0.7 notitle, \
  0 w l notitle

unset label 3



# top right panel
set size pw, ph
set origin pw, ph

set label 1 "(b)" at graph px, py font tagfont

set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
set xlabel "Order {/Arial-Italic n}" offset 2, 0.5

unset logscale y
set mytics 5
set ytics 10 font tcfont offset 0., 0
set format y '%g'
set ylabel '{/Arial-Italic Z_n} / ({/Arial-Italic Z}_{{/Arial-Italic n}-1}{/Arial-Italic V_{D }})' offset 0.0, 0

set key right bottom

#set title "Algorithm G3 ({/Arial-Italic M} = 1) versus algorithm G2, {/Arial-Italic D} = 20"
plot [2:64][0:64] \
  fnG3      u 1:($2)                        w p ls 2 notitle, \
  fnG2      u 1:($2)                        w p ls 6 notitle, \
  -10 w lp ls 2 t "algorithm G3", \
  -10 w lp ls 6 t "algorithm G2", \
  -10 w l notitle

unset xlabel
unset ylabel


# inset of the top right panel
set size iwr, ihr
set origin pw + ixr, ph + iyr

unset label 1

set mxtics 5
set xtics 10 font smtcfont offset 0, 0.5
unset xlabel

unset logscale y
set mytics 5
set ytics 0.001 font smtcfont offset 0.5, 0
set format y "%g"
unset ylabel

set label 3 "{/Symbolic-Oblique e} [{/Arial-Italic Z_n} / ({/Arial-Italic Z}_{{/Arial-Italic n}-1}{/Arial-Italic V_{D }})]" \
  at graph ilx, ily font tcfont

unset title
plot [2:64][:0.0020] \
  fnG3   u 1:($22/abs($2))                    w lp ls 2 ps 0.6 notitle, \
  fnG2   u 1:($26/abs($2))                    w lp ls 6 ps 0.7 notitle, \
  0 w l notitle

unset label 3





# bottom left panel
set size pw, ph
set origin 0, 0

set label 1 "(c)" at graph px, py font tagfont

set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
set xlabel "Order {/Arial-Italic n}" offset 2, 0.5

set ytics 1e9 font tcfont
set format y '%.0t{/Symbol \264}10^{%T}'
set ylabel "Histograms" offset 1.5, 0
set lmargin 7

plot [2:64][0:] \
  fnG3r    u ($1*.5):8    w lp ls 2 ps 0.6 notitle, \
  fnG2     u 1:6          w lp ls 6        notitle, \
  -1 w lp ls 2 t "algorithm G3", \
  -1 w lp ls 6 t "algorithm G2", \
  -1 w l notitle


# bottom right panel
set size pw, ph
set origin pw, 0

set label 1 "(d)" at graph px, py font tagfont

set mxtics 10
set xtics 10 font tcfont offset 0, 0.3
set xlabel "Order {/Arial-Italic n}" offset 2, 0.5

set ytics 0.1 font tcfont
set format y '%g'
set ylabel "Acceptance ratios of\nensemble transitions" \
  offset 1.5, 0
set lmargin 6.3

plot [2:64][0:0.4] \
  fnG3r    u ($1*.5):(($11 + $12)*.5)   w lp ls 2 ps 0.6 notitle, \
  fnG2     u 1:(($9 + $10)*.5)          w lp ls 6        notitle, \
  -1 w lp ls 2 t "algorithm G3", \
  -1 w lp ls 6 t "algorithm G2", \
  -1 w l notitle



unset multiplot

unset output

set terminal wxt
reset



