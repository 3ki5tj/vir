prj=vir

figs=bcacc.pdf vir64.pdf ergo.png

$(prj).pdf: $(prj).tex  $(prj)text.pdf Makefile $(figs)
	pdflatex $(prj).tex

clean:
	rm -rf $(prj).log *.aux *~

