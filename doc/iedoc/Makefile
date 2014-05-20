prj=ie

$(prj).pdf : $(prj).tex $(prj).bib $(prj).bbl
	latex $(prj)
	bibtex $(prj)
	latex $(prj)
	pdflatex $(prj)

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi \
	      *Notes.bib *~
