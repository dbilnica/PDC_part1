#include <iostream>
#include <omp.h>
#include "tsp.cpp"

using namespace std;

void print_roads(const vector<Road>& roads) {
    for (const Road& road : roads) {
        cout << road.firstCity << " " << road.secondCity << " " << road.cost << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_file" << std::endl;
        return 1;
    }
    double exec_time;
    exec_time = -omp_get_wtime();
    double maxValue = stod(argv[2]); //string to double

    vector<Road> roads = parse_inputs(argv[1]);
    print_roads(roads);

    tsp(roads, maxValue);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);


    //print_result(); // to the stdout!

    return 0;
}
