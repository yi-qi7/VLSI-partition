#include "solution.h"

#define EPSILON 0.02
#define INITIAL_TEMP 10000.0
#define FINAL_TEMP 0.001
#define COOLING_RATE 0.995
#define ITERATIONS_PER_TEMP 100
#define RANDOM_SEED 42

void Solution::read_benchmark(Graph &graph, string benchmark_name) {
    ifstream file(benchmark_name);

    if(!file.is_open()) {
        cerr << "Failed to open the file!" << endl;
        exit(-1);
    }

    int edge_num, node_num;
    string line;
    getline(file >> ws, line);
    istringstream iss(line);
    iss >> edge_num;
    iss >> node_num;

    for(int i = 0; i < edge_num; i++) {
        getline(file, line);
        istringstream iss(line);
        int node_id;
        
        Net *net = graph.add_net(i);

        while(iss >> node_id) {
            Node *node = graph.get_or_create_node(node_id);
            node->add_net(net);
            net->add_node(node);
        }
    }

    file.close();
}

void Solution::initialize_partition(Graph &graph, map<int, int> &partition_map, int &size_X, int &size_Y) {
    vector<Node*> nodes = graph.get_nodes();
    int total_size = nodes.size();
    int target_X = total_size / 2;
    
    mt19937 g(RANDOM_SEED);
    shuffle(nodes.begin(), nodes.end(), g);
    
    size_X = 0;
    size_Y = 0;
    
    for(int i = 0; i < total_size; i++) {
        if(i < target_X) {
            partition_map[nodes[i]->get_index()] = 0;
            size_X++;
        } else {
            partition_map[nodes[i]->get_index()] = 1;
            size_Y++;
        }
    }
}

int Solution::calculate_cut_size(Graph &graph, map<int, int> &partition_map) {
    int cut = 0;
    for(Net *net : graph.get_nets()) {
        bool has_X = false, has_Y = false;
        for(Node *node : net->get_nodes()) {
            if(partition_map[node->get_index()] == 0) has_X = true;
            else has_Y = true;
            if(has_X && has_Y) {
                cut++;
                break;
            }
        }
    }
    return cut;
}

int Solution::calculate_delta(Graph &graph, map<int, int> &partition_map, int node_id) {
    int delta = 0;
    int part = partition_map[node_id];
    
    Node *node = nullptr;
    for(Node *n : graph.get_nodes()) {
        if(n->get_index() == node_id) {
            node = n;
            break;
        }
    }
    
    for(Net *net : node->get_nets()) {
        int count_X = 0, count_Y = 0;
        for(Node *n : net->get_nodes()) {
            if(partition_map[n->get_index()] == 0) count_X++;
            else count_Y++;
        }
        
        if(part == 0) {
            if(count_X == 1 && count_Y > 0) delta--;
            else if(count_Y == 0) delta++;
        } else {
            if(count_Y == 1 && count_X > 0) delta--;
            else if(count_X == 0) delta++;
        }
    }
    
    return delta;
}

bool Solution::is_balanced(int size_X, int size_Y, int total_size) {
    double ratio_X = (double)size_X / total_size;
    double ratio_Y = (double)size_Y / total_size;
    return (ratio_X >= 0.5 - EPSILON && ratio_X <= 0.5 + EPSILON &&
            ratio_Y >= 0.5 - EPSILON && ratio_Y <= 0.5 + EPSILON);
}

void Solution::sa_partition(Graph &graph, set<int> &X, set<int> &Y) {
    vector<Node*> nodes = graph.get_nodes();
    int total_size = nodes.size();
    
    mt19937 gen(RANDOM_SEED);
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    map<int, int> partition_map;
    int size_X, size_Y;
    initialize_partition(graph, partition_map, size_X, size_Y);
    
    int current_cut = calculate_cut_size(graph, partition_map);
    
    map<int, int> best_partition = partition_map;
    int best_cut = current_cut;
    
    double temperature = INITIAL_TEMP;
    
    vector<int> node_ids;
    for(Node *n : nodes) {
        node_ids.push_back(n->get_index());
    }
    
    while(temperature > FINAL_TEMP) {
        for(int iter = 0; iter < ITERATIONS_PER_TEMP; iter++) {
            uniform_int_distribution<int> node_dist(0, total_size - 1);
            
            if(dis(gen) < 0.5) {
                int idx = node_dist(gen);
                int node_id = node_ids[idx];
                
                int new_size_X = size_X;
                int new_size_Y = size_Y;
                
                if(partition_map[node_id] == 0) {
                    new_size_X--;
                    new_size_Y++;
                } else {
                    new_size_X++;
                    new_size_Y--;
                }
                
                if(!is_balanced(new_size_X, new_size_Y, total_size)) {
                    continue;
                }
                
                int delta = calculate_delta(graph, partition_map, node_id);
                partition_map[node_id] = 1 - partition_map[node_id];
                
                if(delta < 0 || dis(gen) < exp(-delta / temperature)) {
                    size_X = new_size_X;
                    size_Y = new_size_Y;
                    current_cut += delta;
                    
                    if(current_cut < best_cut) {
                        best_cut = current_cut;
                        best_partition = partition_map;
                    }
                } else {
                    partition_map[node_id] = 1 - partition_map[node_id];
                }
            } else {
                int idx1 = node_dist(gen);
                int idx2 = node_dist(gen);
                
                int node1 = node_ids[idx1];
                int node2 = node_ids[idx2];
                
                if(partition_map[node1] == partition_map[node2]) {
                    continue;
                }
                
                int new_size_X = size_X;
                int new_size_Y = size_Y;
                
                if(partition_map[node1] == 0) {
                    new_size_X--;
                    new_size_Y++;
                } else {
                    new_size_X++;
                    new_size_Y--;
                }
                
                if(partition_map[node2] == 0) {
                    new_size_X--;
                    new_size_Y++;
                } else {
                    new_size_X++;
                    new_size_Y--;
                }
                
                if(!is_balanced(new_size_X, new_size_Y, total_size)) {
                    continue;
                }
                
                int delta1 = calculate_delta(graph, partition_map, node1);
                partition_map[node1] = 1 - partition_map[node1];
                int delta2 = calculate_delta(graph, partition_map, node2);
                partition_map[node2] = 1 - partition_map[node2];
                
                int total_delta = delta1 + delta2;
                
                if(total_delta < 0 || dis(gen) < exp(-total_delta / temperature)) {
                    size_X = new_size_X;
                    size_Y = new_size_Y;
                    current_cut += total_delta;
                    
                    if(current_cut < best_cut) {
                        best_cut = current_cut;
                        best_partition = partition_map;
                    }
                } else {
                    partition_map[node1] = 1 - partition_map[node1];
                    partition_map[node2] = 1 - partition_map[node2];
                }
            }
        }
        
        temperature *= COOLING_RATE;
    }
    
    partition_map = best_partition;
    
    X.clear();
    Y.clear();
    for(auto &p : partition_map) {
        if(p.second == 0) {
            X.insert(p.first);
        } else {
            Y.insert(p.first);
        }
    }
}

void Solution::output_partition(Graph &graph, set<int> &X, string output_file) {
    ofstream out(output_file);
    if(!out.is_open()) {
        cerr << "Failed to open output file!" << endl;
        return;
    }
    
    vector<Node*> nodes = graph.get_nodes();
    sort(nodes.begin(), nodes.end(), [](Node* a, Node* b) {
        return a->get_index() < b->get_index();
    });
    
    for(Node *node : nodes) {
        if(X.find(node->get_index()) != X.end()) {
            out << 0 << endl;
        } else {
            out << 1 << endl;
        }
    }
    
    out.close();
}
