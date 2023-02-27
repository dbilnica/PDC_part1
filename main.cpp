#include <iostream>
#include <omp.h>
#include "tsp.cpp"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_file" << std::endl;
        return 1;
    }
    double exec_time;
    exec_time = -omp_get_wtime();
    double maxValue = stod(argv[2]); //string to double

    vector<vector<int>> distances = parse_inputs(argv[1]);
    tsp(distances, maxValue);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    
    return 0;
}
