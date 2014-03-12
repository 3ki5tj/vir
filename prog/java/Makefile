VirSampApp.class: VirSampApp.java
	javac $<

VirSampApp.java: VirSampApp.0.java MyScheme.java MyCanvas.java Diagram.java \
                 Bits.java Ave.java XYZModel.java Atom.java
	python comb.py

classes=$(wildcard *.class)

run: $(classes)
	appletviewer virsamp.html

jar: VirSampApp.jar

VirSampApp.jar: VirSampApp.class $(classes) Makefile
	jar cf $@ $< $(classes)

clean:
	rm -f *.class *~
