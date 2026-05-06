#ifndef SOLUTION_H
#define SOLUTION_H

#include <string>
#include "Graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

class Solution{
    public:
        void read_benchmark(Graph &graph, string benchmark_name);
        void sa_partition(Graph &graph, set<int> &X, set<int> &Y);
        void output_partition(Graph &graph, set<int> &X, string output_file);
    private:
        int calculate_cut_size(Graph &graph, map<int, int> &partition_map);
        int calculate_delta(Graph &graph, map<int, int> &partition_map, int node_id);
        void initialize_partition(Graph &graph, map<int, int> &partition_map, int &size_X, int &size_Y);
        bool is_balanced(int size_X, int size_Y, int total_size);
};

#endif
