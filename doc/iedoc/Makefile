prj=ievir

figs=ievircmp ievirhigh ievirlow cmpcorl
pdffigs=$(patsubst %,fig/%.pdf, $(figs))
epsfigs=$(patsubst %,fig/%.eps, $(figs))

$(prj)doc.pdf: $(prj)doc.tex $(prj)text.pdf \
  		iegraph.pdf $(pdffigs) Makefile
	make -C fig
	pdflatex $<

iegraph.pdf: iegraph.tex liquid.bib Makefile
	latex $<
	bibtex iegraph
	latex $<
	pdflatex $<

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi \
	      *Notes.bib *~

$(prj)doc.zip: $(prj)doc.tex $(prj)doc.pdf $(prj)text.doc \
  	       iegraph.tex iegraph.pdf liquid.bib \
               $(prj).enl $(prj).Data \
	       $(epsfigs) $(pdffigs) \
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
