#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

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

pair<double,double> getLowestCosts(vector<vector<double>>& distances, int city){
    double lowest1 = numeric_limits<double>::max();
    double lowest2 = numeric_limits<double>::max();

    for(int i = 0; i < distances[city].size(); i++){
        double currentValue = distances[city][i];

        if (currentValue != 0) {
            if (currentValue < lowest1) {
                lowest2 = lowest1;
                lowest1 = currentValue;

            } else if (currentValue < lowest2) {
                lowest2 = currentValue;
            }
        }
    }
    return {lowest1, lowest2};
}


double computeInitLowerBound(vector<pair<double, double>>& lowestCosts) {
    double lowest1 = 0;
    double lowest2 = 0;
    double sum = 0;
    for(int i = 0; i < lowestCosts.size(); i++){
        lowest1 = lowestCosts[i].first;
        lowest2 = lowestCosts[i].second;
        sum += lowest1 + lowest2;
    }
    return sum / 2;
}


double computeLowerBound(vector<vector<double>>& distances,
                         vector<pair<double,double>>& lowestPairs,
                         int city1,
                         int city2,
                         double lb) {
    double cf = 0;
    double ct = 0;
    double cost = distances[city1][city2];

    // CF
    double lowest1 = lowestPairs[city1].first;
    double lowest2 = lowestPairs[city1].second;

    if(cost >= lowest2){
        cf = lowest2;
    }
    else{
        cf = lowest1;
    }

    // CT
    lowest1 = lowestPairs[city2].first;
    lowest2 = lowestPairs[city2].second;

    if(cost >= lowest2){
        ct = lowest2;
    }
    else{
        ct = lowest1;
    }
    return lb + cost - (cf+ct)/2;
}

pair<vector<int>,double> tsp(pair<vector<vector<double>>,vector<pair<double,double>>>& inputs, double maxTourCost){
    vector<vector<double>> distances = get<0>(inputs);
    vector<pair<double,double>> lowestCosts = get<1>(inputs);
    vector<int> bestTour;
    double bestTourCost = maxTourCost;

    // Queue init
    PriorityQueue<Node> queue;

    // Initial lower bound
    double rootLB = computeInitLowerBound(lowestCosts);

    // Root node
    Node root({0},0,rootLB,1,0);
    queue.push(root);

    while(!queue.empty()){
        Node currentNode = queue.pop();

        if(currentNode.lower_bound >= bestTourCost){
            // tour is not complete -> no solution
            if(bestTour.size() < distances[0].size()){
                cout << "NO SOLUTION" << endl;
                return {{},0};
            }
            else{
                bestTour.push_back(0);
                return {bestTour, bestTourCost};
            }
        }

        if(currentNode.nodes == distances[0].size()){
            if(currentNode.cost + distances[currentNode.currentCity][0] < bestTourCost){
                bestTour = currentNode.tour;
                bestTourCost = currentNode.cost + distances[currentNode.currentCity][0];
            }
        }
        else{
            for(int v = 0; v < distances[currentNode.currentCity].size(); v++){
                // Check if city was visited
                if(find(currentNode.tour.begin(), currentNode.tour.end(), v) == currentNode.tour.end())
                {
                    // Skip same cities
                    if(v == currentNode.currentCity || distances[currentNode.currentCity][v] == INFINITY){
                        continue;
                    }
                    double newBound = computeLowerBound(distances,lowestCosts, currentNode.currentCity,
                                                        v, currentNode.lower_bound);

                    // Is higher than best so far
                    if(newBound > bestTourCost){
                        continue;
                    }
                    vector<int> newTour = currentNode.tour;
                    newTour.push_back(v);
                    double newCost = currentNode.cost + distances[currentNode.currentCity][v];

                    Node newNode = Node(newTour, newCost, newBound, currentNode.nodes + 1, v);
                    queue.push(newNode);
                }
            }
        }
    }
    bestTour.push_back(0);
    return {bestTour,bestTourCost};
}

pair<vector<vector<double>>,vector<pair<double,double>>>  parse_inputs(const string& filename) {
    ifstream inputFile(filename);
    int num_cities, roads, row, col;
    double val;
    inputFile >> num_cities >> roads;
    vector<vector<double>> distances(num_cities, vector<double>(num_cities, INFINITY));
    vector<pair<double,double>> lowestCosts(num_cities);

    // fill the array with values from the input file
    while (inputFile >> row >> col >> val) {
        distances[row][col] = val;
    }
    // mirror the array
    for (int i = 0; i < num_cities; i++) {
        for (int j = i+1; j < num_cities; j++) {
            distances[j][i] = distances[i][j];
        }
        pair<double,double> result = getLowestCosts(distances, i);
        lowestCosts[i] = result;
    }
    inputFile.close();
    return {distances, lowestCosts};
}