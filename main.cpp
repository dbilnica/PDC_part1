#include <iostream>
#include <omp.h>
#include <vector>
#include "tsp.cpp"

using namespace std;

void printResults(const pair<vector<int>,double>& v) {
    cout << "\n\nBest tour: ";
    for (int i : v.first) {
        cout << i << " ";
    }
    cout << "Best tour cost: " << v.second << endl;
    cout << endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " input_file" << endl;
        return 1;
    }
    double exec_time;
    exec_time = -omp_get_wtime();
    double maxValue = stod(argv[2]); //string to double

    pair<vector<vector<double>>,vector<pair<double,double>>>  inputs = parse_inputs(argv[1]);
    pair<vector<int>,double> results = tsp(inputs, maxValue);


    // Print results here
    if(!results.first.empty()  && results.second != 0){
        printResults(results);
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);

    return 0;
}