#include <string>
#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <utility>
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


double computeInitLowerBound(vector<pair<pair<int, double>,pair<int, double>>> lowestCosts, int city) {
    pair<int, double> lowest1 = lowestCosts[city].first;
    pair<int, double> lowest2 = lowestCosts[city].second;
    return (lowest1.second + lowest2.second) / 2;
}


double computeLowerBound(vector<vector<int>>& distances,
                         vector<pair<pair<int, double>, pair<int,double>>>& lowestPairs,
                         int city1,
                         int city2,
                         double lb) {
    double cf = 0;
    double ct = 0;
    double cost = distances[city1][city2];

    ////////////////////
    // CF
    pair<int, double> lowest1 = lowestPairs[city1].first;
    pair<int, double> lowest2 = lowestPairs[city1].second;

    // Main condition
    if(cost >= lowest2.second){
        cf = lowest2.second;
    }
    else{
        cf = lowest1.second;
    }

    ////////////////////
    // CT
    lowest1 = lowestPairs[city2].first;
    lowest2 = lowestPairs[city2].second;

    // Main condition
    if(cost >= lowest2.second){
        ct = lowest2.second;
    }
    else{
        ct = lowest1.second;
    }
    return lb + cost - (cf+ct)/2;
}

pair<vector<int>,double> tsp(pair<vector<vector<int>>,vector<pair<pair<int, double>,pair<int, double>>>>& inputs, double maxTourCost){
    vector<vector<int>> distances = get<0>(inputs);
    vector<pair<pair<int, double>,pair<int, double>>> lowestCosts = get<1>(inputs);
    vector<int> bestTour = {};
    double rootLB = 0.0;
    double bestTourCost = maxTourCost;

    // Queue init
    PriorityQueue<Node> queue;

    // Compute lower bound for all children and insert them in queue
    for(int i = 0; i < distances.size(); i++){
        double newBound = computeInitLowerBound(lowestCosts, i);

        cout << newBound << endl;
        rootLB += newBound;

        if(i == 0){
            continue;
        }
        cout << "Node " << i << ": " << distances[0][i] << " | " << newBound << endl;
        Node node({0, i},distances[0][i], newBound, 2, i);
        queue.push(node);
    }

    // Root node
    Node root({0},0,rootLB,1,0);
    queue.push(root);

    // Algorithm
    while(!queue.empty()){
        Node currentNode = queue.pop();

        if(currentNode.lower_bound >= bestTourCost){
            cout << "FINISHED" << endl;
            currentNode.tour.push_back(0);
            return {currentNode.tour, bestTourCost};
        }
        //cout << "Nodes: " << currentNode.nodes << endl;
        if(currentNode.nodes == distances[0].size()){
            if(currentNode.cost + distances[currentNode.currentCity][0] < bestTourCost){
                //bestTour.push_back(currentNode.currentCity);
                bestTour = currentNode.tour;
                cout << "Current city: " << currentNode.currentCity << " Lower Bound: "<< currentNode.lower_bound << " Cost: " << currentNode.cost << endl;
                bestTourCost = currentNode.cost + distances[currentNode.currentCity][0];
            }
        }
        else{
            for(int v = 0; v < distances[currentNode.currentCity].size(); v++){
                // Child nodes already created
                if(currentNode.currentCity == 0){
                    continue;
                }
                //cout << currentNode.cost << endl;
                // Check if city was visited
                if(find(currentNode.tour.begin(), currentNode.tour.end(), v) == currentNode.tour.end())
                {
                    // Skip same cities
                    if(v == currentNode.currentCity){
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

                    cout << "v: " << v << " Cities: "<< currentNode.nodes << " Tour: ";
                    for(int i = 0; i < newTour.size(); i++)
                        cout << newTour[i] << ' ';
                    cout << ""<< endl;
                    double newCost = currentNode.cost + distances[currentNode.currentCity][v];

                    Node newNode = Node(newTour, newCost, newBound, currentNode.nodes + 1, v);
                    queue.push(newNode);
                }
            }
        }
    }
    if(bestTourCost > maxTourCost || bestTour.size() < distances[0].size()){
        cout << "NO SOLUTION" << endl;
        //TODO: fix this
    }
    bestTour.push_back(0);
    return {bestTour,bestTourCost};
}

pair<vector<vector<int>>,vector<pair<pair<int, double>,pair<int, double>>>>  parse_inputs(const string& filename) {
    ifstream inputFile(filename);
    int num_cities, roads;
    inputFile >> num_cities >> roads;

    vector<vector<int>> distances(num_cities, vector<int>(num_cities, 0));
    vector<pair<pair<int, double>,pair<int, double>>> lowestCosts(num_cities);

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
        //TODO: compute all smallest distances here
        pair<pair<int, double>,pair<int, double>> result = getLowestCosts(distances, i);
        lowestCosts[i] = result;
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
    return {distances, lowestCosts};
}