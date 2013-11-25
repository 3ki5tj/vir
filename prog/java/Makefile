VirSampApp.class: VirSampApp.java
	javac $<

classes=$(wildcard *.class)

run: $(classes)
	appletviewer virsamp.html

jar: VirSampApp.jar

VirSampApp.jar: VirSampApp.class $(classes) Makefile
	jar cf $@ $< $(classes)

clean:
	rm -f *.class *~
