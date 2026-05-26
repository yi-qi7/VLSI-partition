// ./main ISPD_benchmark/ibm01.hgr

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <filesystem>
#include "Net.h"
#include "Node.h"
#include "Graph.h"
#include "evaluate.h"
#include "solution.h"

using namespace std;

string get_output_filename(const string& benchmark_path) {
    filesystem::path p(benchmark_path);
    string stem = p.stem().string();
    return stem + "_partition.txt";
}

int main(int argc, char **argv) {

    Solution solution;

    if(argc != 2) {
        cout << "Usage: ./main benchmark_file" << endl;
        exit(-1);
    }
    string benchmark_name = argv[1];
    Graph graph;

    solution.read_benchmark(graph, benchmark_name);    

    cout << "Num nodes: " << graph.get_node_num() << endl;
    cout << "Num nets: " << graph.get_net_num() << endl;

    set<int> X, Y;
    
    cout << "Running FM partition algorithm..." << endl;
    solution.fm_partition(graph, X, Y);
    
    cout << "Partition completed." << endl;
    cout << "Size of partition X: " << X.size() << endl;
    cout << "Size of partition Y: " << Y.size() << endl;
    
    double ratio_X = (double)X.size() / graph.get_node_num();
    double ratio_Y = (double)Y.size() / graph.get_node_num();
    cout << "Ratio X: " << ratio_X << endl;
    cout << "Ratio Y: " << ratio_Y << endl;
    
    string output_file = get_output_filename(benchmark_name);
    solution.output_partition(graph, X, output_file);
    cout << "Partition result saved to: " << output_file << endl;
    
    int cut = calculate_cut(graph, X, Y);
    cout << "Cut size: " << cut << endl;

    return 0;
}
