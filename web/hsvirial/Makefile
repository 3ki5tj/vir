upload::
	gcloud app deploy --project=hsvirial app.yaml

local::
	dev_appserver.py .

prep::
	make -C prog/java jar
	make -C prog/c zip

htmls = $(wildcard *.html)

remote = cz1@he.host.utmb.edu

Bossman: $(htmls) prep VirSampApp.jar
	rsync -avzL *.html *.jar css js images data \
	  prog/java/*.java \
	  --exclude=*~ \
	  $(remote):/Bossman/cz1/vir/web/

check: $(htmls)
	for d in $(htmls); do (aspell -H check $$d ); done

clean:
	find . -name "*~" | xargs rm -f
	rstrip.py -Rv

