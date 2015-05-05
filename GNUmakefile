prj = vir

# portable source code package
# includes zcom.h etc.
$(prj)pack.zip::
	zip $@ GNUmakefile \
	  	prog/*.[ch] prog/README prog/Makefile \
	  	prog/intg/*.c prog/intg/README prog/intg/Makefile \
		prog/samp/*.c prog/samp/README* prog/samp/Makefile \
	  	prog/samp/*.py prog/samp/bak/*.c \
		prog/java/*.java prog/java/*.html prog/java/Makefile \
		prog/java/default* \
		prog/ie/*.[ch] prog/ie/*.py \
		prog/ie/README prog/ie/NOTES \
		data/Bring.dat data/Z.dat data/*.py \
		data/bcacc.txt \
		deploy/stampede/* \
		doc/README doc/Makefile \
		doc/vir.enl doc/vir.Data \
		doc/*.tex doc/*.doc \
		doc/fig/Makefile doc/fig/ergo.cdr \
		doc/fig/*.gp \
		iedoc/README.md iedoc/Makefile \
		iedoc/ievir.enl iedoc/ievir.Data \
		iedoc/*.tex iedoc/*.doc \
		iedoc/fig/Makefile \
		iedoc/fig/*.gp \
		--exclude="*~"
	zip --symlinks $@ \
	  	prog/samp/*.h prog/intg/*.h \
	  	doc/fig/bcacc.txt doc/fig/data


pack: $(prj)pack.zip

$(prj)src.zip::
	git archive --format=zip -9 HEAD > $@

zip: $(prj)src.zip

usbdir = /media/`whoami`/C3/code

usb: $(prj)src.zip $(prj)pack.zip
	mv $(prj)src.zip $(usbdir)/
	mv $(prj)pack.zip $(usbdir)/

#directories with `make clean'
subdirs = prog prog/samp prog/intg

clean:
	$(RM) -f *~ $(prj).o $(prj).zip */*~ */*/*~ */a.out *.tmp
	$(RM) -f $(prj)src.zip $(prj)pack.zip
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

