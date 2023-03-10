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

pair<pair<int, double>,pair<int, double>> getLowestCosts(vector<vector<double>>& distances, int city){
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


double computeInitLowerBound(vector<pair<pair<int, double>,pair<int, double>>>& lowestCosts, int city) {
    double lowest1 = 0;
    double lowest2 = 0;
    double sum = 0;
    for(int i = 0; i < lowestCosts.size(); i++){
        lowest1 = lowestCosts[i].first.second;
        lowest2 = lowestCosts[i].second.second;
        sum += lowest1 + lowest2;
    }
    return sum / 2;
}


double computeLowerBound(vector<vector<double>>& distances,
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

pair<vector<int>,double> tsp(pair<vector<vector<double>>,vector<pair<pair<int, double>,pair<int, double>>>>& inputs, double maxTourCost){
    vector<vector<double>> distances = get<0>(inputs);
    vector<pair<pair<int, double>,pair<int, double>>> lowestCosts = get<1>(inputs);
    vector<int> bestTour;
    double bestTourCost = maxTourCost;

    // Queue init
    PriorityQueue<Node> queue;

    // Initial lower bound
    double rootLB = computeInitLowerBound(lowestCosts, 0);

    // Root node
    Node root({0},0,rootLB,1,0);
    queue.push(root);

    // Algorithm
    while(!queue.empty()){
        Node currentNode = queue.pop();

        if(currentNode.lower_bound >= bestTourCost){
            // tour is not complete -> no solution
            if(bestTour.size() < distances[0].size()){
                cout << "NO SOLUTION" << endl;
                return {{},0};
            }
            else{
                cout << "FINISHED" << endl;
                bestTour.push_back(0);
                return {bestTour, bestTourCost};
            }
        }

        if(currentNode.nodes == distances[0].size()){
            if(currentNode.cost + distances[currentNode.currentCity][0] < bestTourCost){
                bestTour = currentNode.tour;

                cout << "BestTour: ";
                for(int i = 0; i < bestTour.size(); i++)
                    cout << bestTour[i] << ' ';
                cout << " BestTour Cost: " << bestTourCost << endl;
                bestTourCost = currentNode.cost + distances[currentNode.currentCity][0];
            }
        }
        else{
            for(int v = 0; v < distances[currentNode.currentCity].size(); v++){
                // Child nodes already created
                /*
                if(currentNode.currentCity == 0){
                    continue;
                }*/
                //cout << currentNode.cost << endl;
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
                        //cout << "Skipped" << endl;
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

pair<vector<vector<double>>,vector<pair<pair<int, double>,pair<int, double>>>>  parse_inputs(const string& filename) {
    ifstream inputFile(filename);
    int num_cities, roads;
    inputFile >> num_cities >> roads;
    vector<vector<double>> distances(num_cities, vector<double>(num_cities, INFINITY));
    vector<pair<pair<int, double>,pair<int, double>>> lowestCosts(num_cities);

    // fill the array with values from the input file
    int row, col;
    double val;
    while (inputFile >> row >> col >> val) {
        distances[row][col] = val;
    }
    // mirror the array
    for (int i = 0; i < num_cities; i++) {
        for (int j = i+1; j < num_cities; j++) {
            distances[j][i] = distances[i][j];
        }
        pair<pair<int, double>,pair<int, double>> result = getLowestCosts(distances, i);
        lowestCosts[i] = result;
    }

    // print the resulting array

    for (int i = 0; i < num_cities; i++) {
        for (int j = 0; j < num_cities; j++) {
            cout << distances[i][j] << " ";
        }
        cout << endl;
    }
    inputFile.close();
    return {distances, lowestCosts};
}