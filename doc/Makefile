prj=vir

figs=bcacc.pdf vir64.pdf virpy.pdf gperf20.pdf ergo.png reprg.png

$(prj)doc.pdf: $(prj)doc.tex $(prj)text.pdf $(prj)supp.pdf Makefile $(figs)
	pdflatex $<

clean:
	rm -rf $(prj).log $(prj)doc.pdf *.aux *~

$(prj)doc.zip: $(prj)doc.tex $(prj)doc.pdf $(prj)text.doc \
  	       $(prj)supp.doc \
	       paper.pdf \
	       $(prj).enl $(prj).Data $(figs) README Makefile
	zip -r $@ $^

zip: $(prj)doc.zip

mobi: $(prj)text.mobi
	cp $< /media/Kindle/documents

$(prj)text.mobi: $(prj)text.htm
	kindlegen -c1 -verbose $<

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/doc
