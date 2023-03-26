CC = gcc
CFLAGS = -O3 -g -fopenmp
LIBS = -lm

all: tsp

tsp: main.o tsp.o queue.o
	$(CC) -o $@ main.o tsp.o queue.o $(CFLAGS) $(LIBS)

main.o: main.cpp
	$(CC) -c main.cpp $(CFLAGS)

tsp.o: tsp.cpp
	$(CC) -c tsp.cpp $(CFLAGS)

queue.o: queue.hpp
	$(CC) -c queue.hpp $(CFLAGS)

clean:
	rm -f tsp main.o tsp.o queue.o
