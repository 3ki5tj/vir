prj = gaussf

orders =  1  2  3  4  5  6  7  8  9 10 \
         11 12 13 14 15 16 17 18 19 20 \
         21 22 23 24

organize:
	$(foreach ord, $(orders), ./gaussfmerge $(ord);)
	python gfcollect.py

$(prj)data.zip: gaussfD*mpf.dat gaussfn*ldblmpf.dat \
 		gaussfmerge.c gfcollect.py \
		GNUmakefile NOTES
	zip -r $@ $^

zip: $(prj)data.zip

#host=cz1@he.host.utmb.edu
host=cz1@129.109.88.204

syncin:
	rsync -avz $(host):lwork/vir/gaussf/gaussf*.dat .

Bossman: $(prj)data.zip
	rsync -avz $< $(host):/Bossman/cz1/vir/data/

