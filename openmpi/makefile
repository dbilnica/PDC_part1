CXX=mpicc
CXXFLAGS=-g -o

all: tsp

tsp: main.o tsp.o
	$(CXX) $(CXXFLAGS) -o tsp main.o tsp.o -lstdc++

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

tsp.o: tsp.cpp tsp.h
	$(CXX) $(CXXFLAGS) -c tsp.cpp

clean:
	rm -f tsp main.o tsp.o