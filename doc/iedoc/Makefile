prj=ievir

figs=ievircmp ievirhigh ievirlow cmpcorl iepres
pdffigs=$(patsubst %,fig/%.pdf, $(figs))
epsfigs=$(patsubst %,fig/%.eps, $(figs))
figgps=$(patsubst %,fig/%.gp,  $(figs))

latexpdfs = iegraph.pdf $(prj)notes.pdf

all: $(prj)doc.pdf $(latexpdfs)

$(prj)doc.pdf: $(prj)doc.tex $(prj)text.pdf \
  		$(prj)SM1.pdf \
  		iegraph.pdf $(pdffigs) Makefile
	make -C fig
	pdflatex $<

$(latexpdfs): %.pdf : %.tex liquid.bib Makefile
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

iegraph.pdf: fig/iepres.pdf

fig/iepres.pdf:
	$(MAKE) -C fig iegraph.pdf

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi *.out \
	      *Notes.bib *~
	rstrip.py

$(prj)doc.zip: $(prj)doc.tex $(prj)doc.pdf $(prj)text.doc \
               $(prj)SM1.doc \
  	       iegraph.tex iegraph.pdf liquid.bib \
	       $(prj)notes.tex $(prj)notes.pdf \
	       ChangeLog \
	       tables/TableII.html \
	       tables/TableIII.html \
	       tables/TableIV.html \
	       tables/TableV.html \
               $(prj).enl $(prj).Data \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       README Makefile
	zip -r $@ $^

zip: $(prj)doc.zip

mobi: $(prj)text.mobi
	cp $< /media/Kindle/documents

$(prj)text.mobi: $(prj)text.htm
	kindlegen -c1 -verbose $<

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/iedoc
