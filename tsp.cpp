#include <string>
#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class Node {
public:
    vector<int> path;
    double cost;
    double lower_bound;
    int nodes;
    int currentCity;

    Node(vector<int> path, double cost, double lower_bound, int nodes, int currentCity) {
        path = path,
        cost = cost,
        lower_bound = lower_bound,
        nodes = nodes,
        currentCity = currentCity;
    }

    // Compare roads based on their lower bound
    //In case two nodes of the tree have the same lower-bound value, use the index of the city to break
    //the tie, lower indices should be given higher priority. (it is still possible to have a draw, but very
    //unlikely, hence you may safely ignore that situation)
    //tady pouzivat pouze datovy typ double (ne float)

    bool operator>(const Node& other) const {
        if(lower_bound == other.lower_bound){
            // lower indices have higher priority
            return currentCity < other.currentCity;
        }
        return lower_bound > other.lower_bound;
    }
};

// Calculate the lower-bound for the tour
double computeLowerBound(vector<vector<int>>& distances, int city) {
    double lb = 0;
    for(int i = 0; i < distances.size(); i++){

        // Initialize the lowest and second-lowest values
        double lowest1 = numeric_limits<double>::max();
        double lowest2 = numeric_limits<double>::max();

        // Loop through the row and update the lowest and second-lowest values
        for (int currentValue : distances[city]) {
            if (currentValue != 0) {
                if (currentValue < lowest1) {
                    lowest2 = lowest1;
                    lowest1 = currentValue;
                } else if (currentValue < lowest2) {
                    lowest2 = currentValue;
                }
            }
        }
        lb = (lowest1+lowest2) / 2;
    }
    return lb;
}

int tsp(vector<vector<int>>& distances, double bestTourCost){
    // Compute initial lower bound for root node
    double initialLB = computeLowerBound(distances, 0);

    // Queue init
    PriorityQueue<Node> queue;

    // First node
    Node first({0},0.0,initialLB,1,0);
    queue.push(first);

    // Compute lower bound for all children and insert them in queue
    for(int i = 1; i < distances.size(); i++){
        double value = computeLowerBound(distances, i);
        cout << value<< endl;

        // insert nodes
        //TODO: how to set path for children nodes to contain previously visited cities
        //TODO: how to access different nodes in queue
        //
        Node node({0, i},distances[0][i], value, 2, 1);
        queue.push(node);
    }

    /*
    while(!queue.empty()){
        if(first.lower_bound >= bestTourCost){
            // print tour
            cout << "FINISHED" << endl;
            return 0;
            //return make_tuple(first.path, bestTourCost);
        }
        if(first.nodes == distances.size()){
            if(first.cost + distances[0] < bestTourCost){

            }
        }
    }
    */

    return 0;
}

vector<vector<int>> parse_inputs(const string& filename) {
    ifstream inputFile(filename);
    int num_cities, roads;
    inputFile >> num_cities >> roads;

    vector<vector<int>> distances(num_cities, vector<int>(num_cities, 0));

    // fill the array with values from the input file
    int row, col, val;
    while (inputFile >> row >> col >> val) {
        distances[row][col] = val;
    }
    // mirror the array
    for (int i = 0; i < num_cities; i++) {
        for (int j = i+1; j < num_cities; j++) {
            distances[j][i] = distances[i][j];
        }
    }
    // print the resulting array
    for (int i = 0; i < num_cities; i++) {
        for (int j = 0; j < num_cities; j++) {
            cout << distances[i][j] << " ";
        }
        cout << endl;
    }
    inputFile.close();
    return distances;
}