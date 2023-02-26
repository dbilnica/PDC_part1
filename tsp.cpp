#include <string>
#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
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
};

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

string tsp(const vector<City>& roads, double maxValue){
    double bestTourCost = maxValue;

    // Compute initial lower bound
    
    // Queue init
    PriorityQueue<Node> queue;

    // First node
    Node first({0},0.0,1.23,1,0);
    queue.push(first);

    return "o";
}

//make function lower bound
// napr. pro nulu mame lowest bound 1 a 5
// najdeme si cislo v kolecku a najdeme dve nejnizsi cisla (hrany)
// F origin F = 0

void printCitiesAndRoads(const std::vector<City>& cities) {
    for (const auto& city : cities) {
        std::cout << "City " << city.index << ":\n";
        for (const auto& road : city.paths) {
            std::cout << "  Road to city " << (road.firstCity == city.index ? road.secondCity : road.firstCity) << " with cost " << road.cost << "\n";
        }
    }
}

vector<City> parse_inputs(const string& filename) {
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
    }

    printCitiesAndRoads(cities);
    inputFile.close();
    return cities;
}