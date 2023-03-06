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

    // Compare nodes based on their lower bound
    bool operator>(const Node& other) const {
        if(lower_bound == other.lower_bound){
            // lower indices have higher priority
            return currentCity < other.currentCity;
        }
        return lower_bound > other.lower_bound;
    }
};
pair<pair<int, double>,pair<int, double>> getLowestCosts(vector<vector<int>>& distances, int city){
    pair<int, double> lowest1 = {0, numeric_limits<double>::max()};
    pair<int, double> lowest2 = {0, numeric_limits<double>::max()};

    for(int i = 0; i < distances[city].size(); i++){
        double currentValue = distances[city][i];

        if (currentValue != 0) {
            if (currentValue < lowest1.second) {
                // pass values
                lowest2.second = lowest1.second;
                // pass city indexes
                lowest2.first = lowest1.first;

                lowest1.second = currentValue;
                lowest1.first = i;
            } else if (currentValue < lowest2.second) {
                lowest2.second = currentValue;
                lowest2.first = i;
            }
        }
    }
    return {lowest1, lowest2};
}


tuple<double, pair<int,double>,pair<int,double>> computeInitLowerBound(vector<vector<int>>& distances, int city) {
    double lb = 0.0;

    pair<pair<int, double>,pair<int, double>> result = getLowestCosts(distances, city);
    pair<int, double> lowest1 = result.first;
    pair<int, double> lowest2 = result.second;
    lb = (lowest1.second + lowest2.second) / 2;

    return {lb, lowest1, lowest2};
}


double computeLowerBound(vector<vector<int>>& distances,
                         vector<pair<pair<int, double>, pair<int,double>>>& lowestPairs,
                         int city1,
                         int city2,
                         double lb) {
    double cf = 0;
    double ct = 0;
    double cost = distances[city1][city2];
    pair<int, double> lowest1 = {0, 0};
    pair<int, double> lowest2 = {0, 0};

    ////////////////////
    // CF

    // Check if lowest1 and lowest2 are computed
    auto it = find_if(lowestPairs.begin(), lowestPairs.end(),
                           [city1, city2](const auto& p) {
                               return p.first.first == city1 || p.second.first == city2;
                           });
    // Lowest1 and lowest2 were found
    if (it != lowestPairs.end()) {
        // Retrieve the found pair
        pair<pair<int, double>, pair<int,double>> found_pair = *it;
        lowest1 = found_pair.first;
        lowest2 = found_pair.second;
    }
    // Lowest1 and lowest2 were not found
    else{
        // Compute lowest
        pair<pair<int, double>,pair<int, double>> result = getLowestCosts(distances, city1);
        lowest1 = result.first;
        lowest2 = result.second;
    }

    // Main condition
    if(cost >= lowest2.second){
        cf = lowest2.second;
    }
    else{
        cf = lowest1.second;
    }

    ////////////////////
    // CT
    // Check if lowest1 and lowest2 are computed
    it = find_if(lowestPairs.begin(), lowestPairs.end(),
                      [city1, city2](const auto& p) {
                          return p.first.first == city2 || p.second.first == city1;
                      });
    // Lowest1 and lowest2 were found
    if (it != lowestPairs.end()) {
        // Retrieve the found pair
        pair<pair<int, double>, pair<int,double>> found_pair = *it;
        lowest1 = found_pair.first;
        lowest2 = found_pair.second;
    }
    // Lowest1 and lowest2 were not found
    else{
        // Compute lowest
        pair<pair<int, double>,pair<int, double>> result = getLowestCosts(distances, city2);
        lowest1 = result.first;
        lowest2 = result.second;
    }

    // Main condition
    if(cost >= lowest2.second){
        ct = lowest2.second;
    }
    else{
        ct = lowest1.second;
    }

    return lb + cost - (cf+ct)/2;
}

int tsp(vector<vector<int>>& distances, double bestTourCost){
    //Saving lowest cost of travel
    // vector with tuples of lowest costs containing tuples with index of the city and cost of travel
    vector<pair<pair<int, double>, pair<int, double>>> lowestPairs(distances.size());

    // Queue init
    PriorityQueue<Node> queue;
    double rootLB = 0.0;

    // Compute lower bound for all children and insert them in queue
    for(int i = 0; i < distances.size(); i++){
        tuple<double, pair<int,double>,pair<int,double>> result = computeInitLowerBound(distances, i);
        double value = get<0>(result);
        pair<int,double> lowest1 = get<1>(result);
        pair<int,double> lowest2 = get<2>(result);

        cout << value << endl;
        rootLB += value;

        // Saving computed lowest costs
        lowestPairs.push_back({lowest1, lowest2});

        if(i == 0){
            continue;
        }
        cout << "Node " << i << ": " << distances[0][i] << " | " << value << endl;
        Node node({0, i},distances[0][i], value, 2, i);
        queue.push(node);
    }
    cout << rootLB << endl;

    // First node
    Node first({0},0.0,rootLB,1,0);
    queue.push(first);

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
    /*
    for (int i = 0; i < num_cities; i++) {
        for (int j = 0; j < num_cities; j++) {
            cout << distances[i][j] << " ";
        }
        cout << endl;
    }*/
    inputFile.close();
    return distances;
}