#include <iostream>
#include <vector>
#include "tsp.h"
#include <mpi.h>

using namespace std;

void printResults(const pair<vector<int>,double>& v) {
    for (int i : v.first) {
        cout << i << " ";
    }
    cout << endl;
    cout << v.second << endl;
}

int main(int argc, char* argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double exec_time;
    exec_time = -omp_get_wtime();
    double maxValue = stod(argv[2]); //string to double

    pair<vector<vector<double>>,vector<pair<double,double>>>  inputs;

    if(rank == 0) {
        inputs = parse_inputs(argv[1]);
    }

    // Broadcast inputs to all processors
    MPI_Bcast(&inputs, sizeof(inputs), MPI_BYTE, 0, MPI_COMM_WORLD);

    pair<vector<int>,double> results = tsp(inputs, maxValue);

    // Collect results on the root processor
    if(rank == 0) {
        // Print results here
        if(!results.first.empty()  && results.second != 0){
            printResults(results);
        }

        exec_time += omp_get_wtime();
        fprintf(stderr, "%.1fs\n", exec_time);
    }

    MPI_Finalize();
    return 0;
}