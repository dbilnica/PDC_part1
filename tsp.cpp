#include <string>
#include "queue.hpp"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

using namespace std;

class Road{
    public:
        int firstCity;
        int secondCity;
        double cost;
};

class Node {
public:
    int path;
    double cost;
    double lower_bound;
    int nodes;
    int currentCity;

    Node(int path, double cost, double lower_bound, int nodes, int currentCity) {
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

    bool operator>(const Node& other) const {
        if(lower_bound == other.lower_bound){
            // lower indices have higher priority
            return currentCity < other.currentCity;
        }
        return lower_bound > other.lower_bound;
    }
};

string tsp(const vector<Road>& roads, double maxValue){

    // Queue init
    PriorityQueue<Node> queue;

    // First node
    Node first(0,0.0,1.23,1,0);
    queue.push(first);

    return "o";
}

vector<Road> parse_inputs(const string& filename) {
    vector<Road> roads;
    ifstream input_file(filename);
    if (!input_file.is_open()) {
        throw runtime_error("Failed to open input file: " + filename);
    }

    int firstCity, secondCity;
    double cost;
    while (input_file >> firstCity >> secondCity >> cost) {
        Road road;
        road.firstCity = firstCity;
        road.secondCity = secondCity;
        road.cost = cost;
        roads.push_back(road);
    }
    input_file.close();
    return roads;
}