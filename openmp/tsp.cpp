#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>

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
    vector<int> bestTour = {};
    double bestTourCost = maxTourCost;
    int nthreads = omp_get_max_threads();

    // Queues init
    vector<PriorityQueue<Node>> queues(nthreads); // Create one queue per thread

    // Compute starting nodes - root and first layer
    vector<Node> startNodes;
    double lb;

    // Initial lower bound
    double rootLB = computeInitLowerBound(lowestCosts);

    for(int i = 1; i < 4; i++){
        lb = computeLowerBound(distances, lowestCosts,0, i, rootLB);
        double cost = distances[0][i]; // current is 0
        Node newNode = Node({0, i}, cost, lb, 2, i);
        startNodes.push_back(newNode);
    }

    // Insert starting nodes into the thread queues
    for(int i = 0; i < queues.size(); i++) {
        int node = i % startNodes.size();
        queues[i].push(startNodes[node]);
    }

    vector<pair<vector<int>,double>> results; // {bestTour, bestTourCost}
    vector<int> threadTour;
    Node currentNode;
    double newCost;
    vector<int> newTour;
    double newBound;
    int threadId;
    double threadCost;

#pragma omp parallel shared(results, queues, bestTourCost, bestTour) private(newCost, newTour, newBound, currentNode, threadId, threadCost)
    {
        threadId = omp_get_thread_num();
        threadCost = maxTourCost;
        while(true){
            if (queues[threadId].empty()) {
                bool workStolen = false;
                for (int i = 0; i < nthreads; ++i) {
                    int stealQueue = (threadId + i) % nthreads;
                    if (!queues[stealQueue].empty()) {
                        currentNode = queues[stealQueue].pop();
                        workStolen = true;
                        break;
                    }
                }
                if (!workStolen) {
                    break;
                }
            } else {
                currentNode = queues[threadId].pop();
            }
            if(currentNode.lower_bound >= threadCost){
                // tour is not complete -> no solution
                if(threadTour.size() < distances[0].size()){
                    results.push_back({{},0});
                    break;
                }
                else{
                    results.push_back({threadTour, threadCost});
                    break;
                }
            }

            if(currentNode.nodes == distances[0].size()){
                if(currentNode.cost + distances[currentNode.currentCity][0] < threadCost){
                    threadTour = currentNode.tour;
#pragma omp atomic write
                    threadCost = currentNode.cost + distances[currentNode.currentCity][0];
                }
            }
            else{
#pragma omp parallel for
                for(int v = 0; v < distances[currentNode.currentCity].size(); v++){
                    // Check if city was visited
                    if(find(currentNode.tour.begin(), currentNode.tour.end(), v) == currentNode.tour.end())
                    {
                        // Skip same cities
                        if(v == currentNode.currentCity || distances[currentNode.currentCity][v] == INFINITY){
                            continue;
                        }
                        newBound = computeLowerBound(distances,lowestCosts, currentNode.currentCity,
                                                     v, currentNode.lower_bound);

                        // Is higher than best so far
                        if(newBound > threadCost){
                            continue;
                        }
                        newTour = currentNode.tour;
                        newTour.push_back(v);
#pragma omp atomic write
                        newCost = currentNode.cost + distances[currentNode.currentCity][v];

                        Node newNode = Node(newTour, newCost, newBound, currentNode.nodes + 1, v);
#pragma omp critical
                        {
                            queues[threadId].push(newNode);
                        };
                    }
                }
            }
        }
        // Compare results
#pragma omp master
#pragma omp critical
        {
            int best = 0;
            double cost = results[0].second;
            for (int i = 0; i < results.size(); i++) {
                if (results[i].second < cost) {
                    best = i;
                    cost = results[i].second;
                }
            }
            bestTour = results[best].first;
            bestTourCost = cost;
        }
    }

    if(bestTour.size() < distances[0].size()){
        cout << "NO SOLUTION" << endl;
    }
    else{
        bestTour.push_back(0);
    }
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