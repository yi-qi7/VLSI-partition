//  ./main ISPD_benchmark/ibm01.hgr

#include "solution.h"
#include <omp.h>

#define EPSILON 0.02
#define MAX_PASSES 50
#define MULTI_START 20

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

void Solution::initialize_partition(Graph &graph, map<int, int> &partition_map, 
                                    int &size_X, int &size_Y, unsigned int seed) {
    vector<Node*> nodes = graph.get_nodes();
    int total_size = nodes.size();
    int target_X = total_size / 2;
    
    mt19937 g(seed);
    shuffle(nodes.begin(), nodes.end(), g);
    
    size_X = 0;
    size_Y = 0;
    partition_map.clear();
    
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

void Solution::initialize_partition_bfs(Graph &graph, map<int, int> &partition_map,
                                        int &size_X, int &size_Y, unsigned int seed) {
    vector<Node*> nodes = graph.get_nodes();
    int total_size = nodes.size();
    int target_X = total_size / 2;
    
    mt19937 g(seed);
    uniform_int_distribution<int> dist(0, total_size - 1);
    int start_node_idx = dist(g);
    
    partition_map.clear();
    vector<bool> visited(total_size + 1, false);
    queue<Node*> q;
    
    Node* start_node = nodes[start_node_idx];
    q.push(start_node);
    visited[start_node->get_index()] = true;
    
    size_X = 0;
    size_Y = 0;
    
    while(!q.empty() && size_X < target_X) {
        Node* current = q.front();
        q.pop();
        
        if(size_X < target_X) {
            partition_map[current->get_index()] = 0;
            size_X++;
        }
        
        vector<Node*> neighbors;
        for(Net* net : current->get_nets()) {
            for(Node* neighbor : net->get_nodes()) {
                if(!visited[neighbor->get_index()]) {
                    neighbors.push_back(neighbor);
                }
            }
        }
        
        shuffle(neighbors.begin(), neighbors.end(), g);
        for(Node* neighbor : neighbors) {
            if(!visited[neighbor->get_index()]) {
                visited[neighbor->get_index()] = true;
                q.push(neighbor);
            }
        }
    }
    
    for(Node* node : nodes) {
        if(partition_map.find(node->get_index()) == partition_map.end()) {
            partition_map[node->get_index()] = 1;
            size_Y++;
        }
    }
}

int Solution::calculate_gain(Node *node, map<int, int> &partition_map) {
    int gain = 0;
    int node_part = partition_map[node->get_index()];
    
    for(Net *net : node->get_nets()) {
        int count_in_X = 0, count_in_Y = 0;
        
        for(Node *n : net->get_nodes()) {
            if(partition_map[n->get_index()] == 0) {
                count_in_X++;
            } else {
                count_in_Y++;
            }
        }
        
        if(node_part == 0) {
            if(count_in_X == 1 && count_in_Y > 0) gain++;
            if(count_in_Y == 0) gain--;
        } else {
            if(count_in_Y == 1 && count_in_X > 0) gain++;
            if(count_in_X == 0) gain--;
        }
    }
    
    return gain;
}

int Solution::calculate_cut_size(Graph &graph, map<int, int> &partition_map) {
    int cut = 0;
    vector<Net*> nets = graph.get_nets();
    
    for(Net* net : nets) {
        bool has_in_X = false, has_in_Y = false;
        for(Node* node : net->get_nodes()) {
            if(partition_map[node->get_index()] == 0) has_in_X = true;
            else has_in_Y = true;
            if(has_in_X && has_in_Y) {
                cut++;
                break;
            }
        }
    }
    
    return cut;
}

bool Solution::is_balanced(int size_X, int size_Y, int total_size, int moving_from) {
    int new_size_X = size_X;
    int new_size_Y = size_Y;
    
    if(moving_from == 0) {
        new_size_X--;
        new_size_Y++;
    } else {
        new_size_X++;
        new_size_Y--;
    }
    
    double ratio_X = (double)new_size_X / total_size;
    
    return (ratio_X >= 0.5 - EPSILON && ratio_X <= 0.5 + EPSILON);
}

void Solution::update_gains_incremental(Node *moved_node, map<int, int> &partition_map, 
                                       map<int, int> &gains, vector<bool> &locked) {
    int moved_part = partition_map[moved_node->get_index()];
    int from_part = 1 - moved_part;
    
    for(Net *net : moved_node->get_nets()) {
        int count_in_from = 0, count_in_to = 0;
        Node* critical_node_from = nullptr;
        Node* critical_node_to = nullptr;
        
        for(Node *n : net->get_nodes()) {
            if(partition_map[n->get_index()] == from_part) {
                count_in_from++;
                if(count_in_from == 1) critical_node_from = n;
            } else {
                count_in_to++;
                if(count_in_to == 1) critical_node_to = n;
            }
        }
        
        if(count_in_from == 0) {
            for(Node *n : net->get_nodes()) {
                if(!locked[n->get_index()] && partition_map[n->get_index()] == moved_part) {
                    gains[n->get_index()]--;
                }
            }
        } else if(count_in_from == 1) {
            if(critical_node_from && !locked[critical_node_from->get_index()]) {
                gains[critical_node_from->get_index()]++;
            }
        }
        
        if(count_in_to == 0) {
            for(Node *n : net->get_nodes()) {
                if(!locked[n->get_index()] && partition_map[n->get_index()] == from_part) {
                    gains[n->get_index()]++;
                }
            }
        } else if(count_in_to == 1) {
            if(critical_node_to && !locked[critical_node_to->get_index()]) {
                gains[critical_node_to->get_index()]--;
            }
        }
    }
}

void Solution::fm_pass(Graph &graph, map<int, int> &partition_map, int &size_X, int &size_Y) {
    vector<Node*> nodes = graph.get_nodes();
    int total_size = nodes.size();
    
    map<int, int> gains;
    vector<bool> locked(total_size + 1, false);
    int locked_count = 0;
    
    for(Node *node : nodes) {
        gains[node->get_index()] = calculate_gain(node, partition_map);
    }
    
    int current_size_X = size_X;
    int current_size_Y = size_Y;
    map<int, int> best_partition = partition_map;
    int cumulative_gain = 0;
    int best_cumulative = 0;
    
    while(locked_count < total_size) {
        int best_gain = INT_MIN;
        int best_node_idx = -1;
        
        for(Node *node : nodes) {
            int idx = node->get_index();
            if(!locked[idx]) {
                int node_part = partition_map[idx];
                if(is_balanced(current_size_X, current_size_Y, total_size, node_part)) {
                    if(gains[idx] > best_gain) {
                        best_gain = gains[idx];
                        best_node_idx = idx;
                    }
                }
            }
        }
        
        if(best_node_idx == -1) break;
        
        int old_part = partition_map[best_node_idx];
        partition_map[best_node_idx] = 1 - old_part;
        
        if(old_part == 0) {
            current_size_X--;
            current_size_Y++;
        } else {
            current_size_X++;
            current_size_Y--;
        }
        
        locked[best_node_idx] = true;
        locked_count++;
        cumulative_gain += best_gain;
        
        if(cumulative_gain > best_cumulative) {
            best_cumulative = cumulative_gain;
            best_partition = partition_map;
        }
        
        Node *best_node = nullptr;
        for(Node *n : nodes) {
            if(n->get_index() == best_node_idx) {
                best_node = n;
                break;
            }
        }
        
        update_gains_incremental(best_node, partition_map, gains, locked);
    }
    
    partition_map = best_partition;
    size_X = 0;
    size_Y = 0;
    for(auto &p : partition_map) {
        if(p.second == 0) size_X++;
        else size_Y++;
    }
}

void Solution::fm_partition(Graph &graph, set<int> &X, set<int> &Y) {
    vector<Node*> nodes = graph.get_nodes();
    
    map<int, int> best_partition;
    int best_cut = INT_MAX;
    
    int completed_runs = 0;
    int total_runs = MULTI_START;
    
    cout << "\nProgress: [";
    int bar_width = 50;
    for(int i = 0; i < bar_width; i++) cout << " ";
    cout << "] 0/" << total_runs << " runs" << flush;
    
    #pragma omp parallel
    {
        map<int, int> local_partition;
        int local_size_X, local_size_Y;
        int local_best_cut = INT_MAX;
        map<int, int> local_best_partition;
        
        #pragma omp for schedule(dynamic)
        for(int run = 0; run < MULTI_START; run++) {
            unsigned int seed = 42 + run * 1000;
            
            if(run % 2 == 0) {
                initialize_partition(graph, local_partition, local_size_X, local_size_Y, seed);
            } else {
                initialize_partition_bfs(graph, local_partition, local_size_X, local_size_Y, seed);
            }
            
            for(int pass = 0; pass < MAX_PASSES; pass++) {
                int prev_cut = calculate_cut_size(graph, local_partition);
                fm_pass(graph, local_partition, local_size_X, local_size_Y);
                int new_cut = calculate_cut_size(graph, local_partition);
                
                if(new_cut < local_best_cut) {
                    local_best_cut = new_cut;
                    local_best_partition = local_partition;
                }
                
                if(new_cut >= prev_cut && pass > 5) {
                    break;
                }
            }
            
            #pragma omp critical
            {
                completed_runs++;
                int pos = bar_width * completed_runs / total_runs;
                
                cout << "\rProgress: [";
                for(int i = 0; i < bar_width; i++) {
                    if(i < pos) cout << "=";
                    else if(i == pos) cout << ">";
                    else cout << " ";
                }
                cout << "] " << completed_runs << "/" << total_runs << " runs" << flush;
            }
        }
        
        #pragma omp critical
        {
            if(local_best_cut < best_cut) {
                best_cut = local_best_cut;
                best_partition = local_best_partition;
            }
        }
    }
    
    cout << endl << endl;
    
    X.clear();
    Y.clear();
    for(auto &p : best_partition) {
        if(p.second == 0) {
            X.insert(p.first);
        } else {
            Y.insert(p.first);
        }
    }
    
    cout << "Best cut size found: " << best_cut << endl;
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
