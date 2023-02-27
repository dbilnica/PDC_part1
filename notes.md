Class city {
Int index;
double lb; // two mins, select lower one
Path paths[] = {Path2, Path3…}
}

Class Path {
City first;
City second;
double cost;
}

tour = {0,2,3,1} - indexes of visited cities, with tour we also save current tours cost

For path in paths:
	Calculate lower bound for each path of the city
Then move to next city

* Use double for everything
* Every city can be visited only once
* Tour, ukladam si indexy kde jsem byl (budu vedet, ze tama uz nemam jit)
* Prunning should be in the queue


#include <string>
#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>

using namespace std;

class Path{
public:
int firstCity;
int secondCity;
double cost;
};

class City{
public:
explicit City(int index) : index(index) {}
int index;
double lowerBound = 0;
vector<Path> paths;

    // Add a method to compute the lower-bound for this city
    double computeLowerBound() const {
        if (paths.size() < 2) {
            // City has less than two paths, so lower-bound is zero
            return 0;
        }
        double lb = paths[0].cost + paths[1].cost;
        double min1 = paths[0].cost;
        double min2 = paths[1].cost;
        for (int i = 2; i < paths.size(); i++) {
            double cost = paths[i].cost;
            if (cost < min1) {
                min2 = min1;
                min1 = cost;
            } else if (cost < min2) {
                min2 = cost;
            }
            lb += cost;
        }
        // Divide by 2 since each edge is counted twice
        lb /= 2;
        return lb - (min1 + min2) / 2;
    }
};

void printCitiesAndRoads(const vector<City>& cities) {
for (const auto& city : cities) {
std::cout << "City " << city.index << " | Lower Bound: " << city.lowerBound << ":\n";
for (const auto& road : city.paths) {
std::cout << "  Road to city " << (road.firstCity == city.index ? road.secondCity : road.firstCity) << " with cost " << road.cost << "\n";
}
}
}

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

// Calculate the initial lower-bound for the tour
double computeInitialLowerBound(const City& city) {
double lb = 0;

    double min1 = INFINITY;
    double min2 = INFINITY;
    for (const Path& road : city.paths) {
        if (road.cost < min1) {
            min2 = min1;
            min1 = road.cost;
        } else if (road.cost < min2) {
            min2 = road.cost;
        }
    }
    cout << min1 << min2 << endl;
    lb += min1 + min2;

    // Divide by 2 since each edge is counted twice
    lb /= 2;
    return lb;
}
int tsp(vector<vector<int>>& distances, double maxValue){
double bestTourCost = maxValue;
double initialLB = 0;
// Compute initial lower bound for children



    /*
    for (City& city : cities) {
        city.lowerBound = computeInitialLowerBound(city);
    }
    printCitiesAndRoads(cities);*/

    // Queue init
    PriorityQueue<Node> queue;

    // First node
    Node first({0},0.0,initialLB,1,0);
    queue.push(first);

    while(!queue.empty()){
        if(first.lower_bound >= bestTourCost){
            // print tour
            cout << "FINISHED" << endl;
            return 0;
            //return make_tuple(first.path, bestTourCost);
        }
        if(first.nodes == cities.size()){
            //if(first.cost)
        }

    }

    return 0;
}

//make function lower bound
// napr. pro nulu mame lowest bound 1 a 5
// najdeme si cislo v kolecku a najdeme dve nejnizsi cisla (hrany)
// F origin F = 0



vector<vector<int>> parse_inputs(const string& filename) {/*
ifstream inputFile(filename);
if (!inputFile.is_open()) {
throw runtime_error("Failed to open input file: " + filename);
}

    int numCities, numRoads;
    inputFile >> numCities >> numRoads;
    vector<City> cities;

    for (int i = 0; i < numCities; i++) {
        cities.emplace_back(i);
    }

    for (int i = 0; i < numRoads; i++) {
        int firstCity, secondCity;
        double cost;
        inputFile >> firstCity >> secondCity >> cost;
        Path path = {firstCity, secondCity, cost};
        cities[firstCity].paths.push_back(path);
        cities[secondCity].paths.push_back(path);
    }*/
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