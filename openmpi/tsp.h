#ifndef TRAVELING_SALESPERSON_TSP_H
#define TRAVELING_SALESPERSON_TSP_H

#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <limits>

using namespace std;

class Node {
public:
    vector<int> tour = {};
    double cost = 0.0;
    double lower_bound = 0.0;
    int nodes = 0;
    int currentCity = 0;

    Node() : cost(0), lower_bound(0), nodes(0), currentCity(0) {}

    Node(vector<int> tour, double cost, double lower_bound, int nodes, int currentCity) :
            tour(move(tour)),
            cost(cost),
            lower_bound(lower_bound),
            nodes(nodes),
            currentCity(currentCity) {}

    // Compare nodes based on their lower bound
    bool operator>(const Node& other) const {
        if(lower_bound == other.lower_bound){
            // lower indices have higher priority
            return currentCity < other.currentCity;
        }
        return lower_bound > other.lower_bound;
    }
};

pair<double,double> getLowestCosts(vector<vector<double>>& distances, int city);

double computeInitLowerBound(vector<pair<double, double>>& lowestCosts);

double computeLowerBound(vector<vector<double>>& distances,
                         vector<pair<double,double>>& lowestPairs,
                         int city1,
                         int city2,
                         double lb);

pair<vector<int>,double> tsp(pair<vector<vector<double>>,vector<pair<double,double>>>& inputs, double maxTourCost);

pair<vector<vector<double>>,vector<pair<double,double>>>  parse_inputs(const string& filename);

#endif //TRAVELING_SALESPERSON_TSP_H
