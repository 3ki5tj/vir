prj=ievir

figs=ievircmp.pdf ievirhigh.pdf ievirlow.pdf ljievir.pdf

$(prj)doc.pdf: $(prj)doc.tex $(prj)text.pdf iegraph.pdf $(figs)
	make -C fig
	pdflatex $<

iegraph.pdf: iegraph.tex $(prj).bib
	latex $<
	bibtex $*
	latex $<
	pdflatex $<

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi \
	      *Notes.bib *~

$(prj)doc.zip: $(prj)doc.pdf $(prj)text.doc iegraph.pdf \
               $(prj).enl $(prj).Data $(figs) README
	zip -r $@ $^

zip: $(prj)doc.zip

mobi: $(prj)text.mobi
	cp $< /media/Kindle/documents

$(prj)text.mobi: $(prj)text.htm
	kindlegen -c1 -verbose $<

