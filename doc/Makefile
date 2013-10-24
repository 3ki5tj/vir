prj=vir

figs=bcacc.pdf vir64.pdf ergo.png

$(prj)doc.pdf: $(prj)doc.tex  $(prj)text.pdf Makefile $(figs)
	pdflatex $<

clean:
	rm -rf $(prj).log *.aux *~

$(prj)doc.zip: $(prj)doc.pdf $(prj)text.doc $(prj).enl $(prj).Data $(figs) README
	zip -r $@ $^

zip: $(prj)doc.zip
