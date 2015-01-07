prj = VirSamp

javas = $(prj)App.0.java MCSamp.java \
	Diagram.java DiagramMap.java Bits.java \
  	XYScheme.java XYZCanvas.java XYZModel.java \
	XYZCanvasVir.java XYZModelMC.java \
	Atom.java \
        Ave.java MDS.java

all: jar

$(prj)App.class: $(prj)App.java
	javac $<

$(prj)App.java: $(javas)
	python comb.py

classes=$(wildcard *.class)

run: $(classes) $(prj)App.class
	appletviewer virsamp.html

jar: $(prj)App.jar

$(prj)App.jar: $(prj)App.class $(classes) $(javas) Makefile
	jar cf $@ $< $(classes)

clean:
	rm -f *.class *~
