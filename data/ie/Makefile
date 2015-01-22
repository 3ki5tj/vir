prj = ie
zipfile=$(prj)data.zip

zip: $(zipfile)

$(zipfile): *.dat Makefile
	zip -9 $@ *.dat \
	  TableII/*.dat TableII/*.sh TabelII/*.py TableII/README* \
	  TabelII/*.html TableII/*.c TableII/*.h \
	  ybg/*.dat kirk/*.dat \
	  hnc/*.dat py/*.dat pyhnc/*.dat \
	  cr/*.dat err/*.dat iebenchmark/*.dat \
	  thermo/*.dat lamc/*.dat \
	  hurst/*.dat ir/*.dat rowlinson/*.dat \
	  bpgg/*.dat mp/*.dat ry/*.dat \
	  Makefile

#host=cz1@he.host.utmb.edu
host=cz1@129.109.88.204

Bossman: $(zipfile)
	rsync -avz $< $(host):/Bossman/cz1/vir/data/
