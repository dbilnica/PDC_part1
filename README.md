# Part 1 - Serial implementation

Write a serial implementation of the algorithm in C or C++. Name the source file of this implementation tsp.c. As stated, your program should expect two input parameters.
Make sure to include a Makefile, the simple command

`$ make`

Should generate an executable with the same name as the source file (minus the extension). Also, to
be uniform across groups, in this makefile please use the gcc compiler with optimization flag -O3.
This version will serve as your base for comparisons and must be as efficient as possible.


# Part 2 - OpenMP implementation
Write an OpenMP implementation of the algorithm, with the same rules and input/output descriptions.
Name this source code tsp-omp.c. You can start by simply adding OpenMP directives, but you are
free, and encouraged, to modify the code in order to make the parallelization as effective and as
scalable as possible. Be careful about synchronization and load balancing!
Important note: in order to test for scalability, we will run this program assigning different values
to the shell variable OMP NUM THREADS. Please do not override this variable or set the number of
threads in your program, otherwise we will not be able to properly evaluate it