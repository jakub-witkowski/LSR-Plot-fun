CC=g++
CFLAGS=-std=c++11 `root-config --cflags`
LDFLAGS=`root-config --libs`
SOURCES=plot.cpp

all: compile run

compile: $(SOURCES)
	$(CC) $(SOURCES) $(LDFLAGS) -o plot $(CFLAGS)
run:
	./plot 1050.csv 4