JCC = javac
JFLAGS = -g

default: RandMST.class EuclideanDistExp.class Experiment.class RandDistExp.class Vertex.class Script.class

randmst: RandMST.class EuclideanDistExp.class Experiment.class RandDistExp.class Vertex.class

RandMST.class: RandMST.java
	$(JCC) $(JFLAGS) RandMST.java

EuclideanDistExp.class: EuclideanDistExp.java
	$(JCC) $(JFLAGS) EuclideanDistExp.java

Experiment.class: Experiment.java
	$(JCC) $(JFLAGS) Experiment.java

RandDistExp.class: RandDistExp.java
	$(JCC) $(JFLAGS) RandDistExp.java

Vertex.class: Vertex.java
	$(JCC) $(JFLAGS) Vertex.java

Script.class: Script.java
	$(JCC) $(JFLAGS) Script.java

clean: 
	$(RM) *.class
