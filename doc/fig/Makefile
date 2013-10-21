
PS = $(patsubst %.gp,%.ps, $(wildcard *.gp))
PDF = $(patsubst %.gp,%.pdf, $(wildcard *.gp))
PNG = $(patsubst %.gp,%.png, $(wildcard *.gp))
deps = Makefile bcacc.txt

default: $(PDF)

all: $(PDF) $(PNG)

$(PNG) : %.png : %.pdf $(deps)
	convert -rotate 90 -density 144 $< $@

# the converted PNG is kind of small, so we set the density 144 dots per inch
# see more in http://www.imagemagick.org/script/command-line-options.php

$(PDF) : %.pdf : %.ps $(deps)
	epstopdf $<

$(PS) : %.ps : %.gp $(deps)
	gpload $<

clean:
	rm -f $(PDF) $(PS) $(PNG) *~
	rstrip.py -Rv

