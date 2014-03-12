prj=biconnected

mobi: $(prj).mobi
	cp $< /media/Kindle/documents

$(prj).mobi: $(prj).htm
	kindlegen -c1 -verbose $<
