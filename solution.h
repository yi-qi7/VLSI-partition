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
#include <queue>
#include <cmath>
#include <algorithm>
#include <random>
#include <climits>

using namespace std;

class Solution{
    public:
        void read_benchmark(Graph &graph, string benchmark_name);
        void fm_partition(Graph &graph, set<int> &X, set<int> &Y);
        void output_partition(Graph &graph, set<int> &X, string output_file);
    private:
        int calculate_gain(Node *node, map<int, int> &partition_map);
        int calculate_cut_size(Graph &graph, map<int, int> &partition_map);
        bool is_balanced(int size_X, int size_Y, int total_size, int moving_from);
        void initialize_partition(Graph &graph, map<int, int> &partition_map, 
                                 int &size_X, int &size_Y, unsigned int seed);
        void initialize_partition_bfs(Graph &graph, map<int, int> &partition_map,
                                      int &size_X, int &size_Y, unsigned int seed);
        void update_gains_incremental(Node *moved_node, map<int, int> &partition_map, 
                                     map<int, int> &gains, vector<bool> &locked);
        void fm_pass(Graph &graph, map<int, int> &partition_map, int &size_X, int &size_Y);
};

#endif
