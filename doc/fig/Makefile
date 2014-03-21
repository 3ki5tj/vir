GP = $(wildcard *.gp)
EPS = $(patsubst %.gp,%.eps, $(GP))
PDF = $(patsubst %.gp,%.pdf, $(GP))
GP_PNG = $(patsubst %.gp,%.png, $(GP))
#PNG = $(wildcard *.png)
PNG = reprg.png ergo.png
TIF = $(patsubst %.png,%.tif, $(PNG))
deps = Makefile bcacc.txt

default: $(PDF) $(TIF)

all: $(PDF) $(GP_PNG) $(TIF)

$(GP_PNG) : %.png : %.pdf $(deps)
	convert -rotate 90 -density 144 $< $@

# the converted PNG is kind of small, so we set the density 144 dots per inch
# see more in http://www.imagemagick.org/script/command-line-options.php

$(PDF) : %.pdf : %.eps $(deps)
	epstopdf $<

$(EPS) : %.eps : %.gp $(deps)
	gpload $<

$(TIF) : %.tif : %.png $(deps)
	convert $< $@
clean:
	rm -f $(PDF) $(PS) $(GP_PNG) *~
	rstrip.py -Rv

