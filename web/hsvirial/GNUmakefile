prj=~/work/vir/web/hsvirial

local: prep
	cd ~/software/google_appengine && \
	echo "the website can be viewed as localhost:8080" && \
	./dev_appserver.py $(prj)

upload: prep
	cd ~/software/google_appengine && \
	./appcfg.py update $(prj)

prep::
	make -C prog/java jar
	make -C prog/c zip

htmls = $(wildcard *.html)

check: $(htmls)
	for d in $(htmls); do (aspell -H check $$d ); done

clean:
	find . -name "*~" | xargs rm -f
	rstrip.py -Rv

