prj = vir
prjsrc = $(prj)src

#directories with `make clean'
subdirs = prog prog/samp prog/intg

clean:
	$(RM) -f *~ $(prj).o $(prj).zip */*~ */*/*~ */a.out *.tmp
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

$(prjsrc).zip::
	git archive --format=zip -9 HEAD > $@

usbdir = /media/C3/code

usb: $(prjsrc).zip
	mv $(prjsrc).zip $(usbdir)/

