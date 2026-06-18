/**
 * @file    topopart_initial.cpp
 * @brief   TopoPart Module 4: Initial Partition (paper=grow dual-mode)
 */

#include "topopart_initial.h"
#include "topopart_utils.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <unordered_set>
#include <queue>
using namespace std;

// ================================================================
//  EvaluateCutDelta
// ================================================================
long long EvaluateCutDelta(const CoarsenLevel& level,
                            const FPGAGraph& fg,
                            const PartitionResult& part,
                            int vi, int target_f) {
    long long delta = 0;
    int old_f = part.node2fpga[vi];
    for (int nb : level.super_graph.adj[vi]) {
        int nb_f = part.node2fpga[nb];
        if (nb_f < 0) continue;
        bool was_cut = (old_f >= 0) && (old_f != nb_f);
        bool new_cut = (target_f != nb_f);
        if (!was_cut && new_cut) delta += 1;
        else if (was_cut && !new_cut) delta -= 1;
    }
    return delta;
}

// ================================================================
//  CountTopoViolations
// ================================================================
static int CountTopoViolations(int vi, int target_f,
                                const CoarsenLevel& level,
                                const FPGAGraph& fg,
                                const PartitionResult& part) {
    int violations = 0;
    for (int vk : level.super_graph.adj[vi]) {
        int v_hat_k = part.node2fpga[vk];
        if (v_hat_k < 0) continue;
        if (fg.dist_all[target_f][v_hat_k] > 1) ++violations;
    }
    return violations;
}

// Forward declarations
static PartitionResult InitialPartitionPaper(CoarsenLevel&, const FPGAGraph&, int);
static PartitionResult InitialPartitionGrow(CoarsenLevel&, const FPGAGraph&);

// ================================================================
//  InitialPartition dispatcher
// ================================================================
PartitionResult InitialPartition(CoarsenLevel& coarsest_level,
                                  const FPGAGraph& fg,
                                  int max_traceback,
                                  const string& mode) {
    if (mode == "grow") {
        return InitialPartitionGrow(coarsest_level, fg);
    }
    return InitialPartitionPaper(coarsest_level, fg, max_traceback);
}

// ================================================================
//  Paper Mode: ICCAD 2021 Algorithm 2 (priority queue + traceback)
// ================================================================
static PartitionResult InitialPartitionPaper(CoarsenLevel& coarsest_level,
                                              const FPGAGraph& fg,
                                              int max_traceback) {
    int N = coarsest_level.super_graph.node_num;
    int K = fg.fpga_num;
    auto& cddt = coarsest_level.super_cddt;
    auto& graph = coarsest_level.super_graph;
    auto& weight = coarsest_level.super_weight;

    cout << "\n========== Module 4: InitialPartition (Paper: Algorithm 2 + Traceback) ==========" << endl;
    cout << "  Super-nodes: " << N << ", FPGAs: " << K << endl;

    PartitionResult part;
    part.init(N, K);

    auto wload_of = [&](int f) -> long long {
        long long wl = 0;
        for (int nid : part.fpga2nodes[f]) wl += (nid < N) ? weight[nid] : 1;
        return wl;
    };

    // Step 1: Assign fixed nodes
    for (int vi : graph.fixed_nodes) {
        part.assign(vi, graph.node2fpga.at(vi));
    }

    // Step 2: Build priority queue Q (smallest candidate set first)
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> Q;
    for (int vi = 0; vi < N; ++vi) {
        if (part.node2fpga[vi] >= 0) continue;
        int sz = cddt[vi].size();
        if (sz == 0) throw runtime_error("Node " + to_string(vi) + " has empty Cddt!");
        Q.push({sz, vi});
    }

    // Step 3: Main assignment loop with traceback
    int assigned = 0, tracebacks = 0;
    while (!Q.empty()) {
        auto [sz, vi] = Q.top(); Q.pop();
        if (part.node2fpga[vi] >= 0) continue;

        int cur_sz = cddt[vi].size();
        if (cur_sz != sz) {
            if (cur_sz == 0) goto do_traceback;
            Q.push({cur_sz, vi}); continue;
        }

        // Find best FPGA: min cut delta, topology-valid
        {
            int best_f = -1;
            long long best_wl = LLONG_MAX;
            int best_viol = INT_MAX;

            uint64_t bits = cddt[vi].bits;
            while (bits) {
                int f; _BitScanForward64((unsigned long*)&f, bits);
                bits &= bits - 1;

                if (wload_of(f) + 1 > fg.resource_cap[f]) continue;
                int viol = CountTopoViolations(vi, f, coarsest_level, fg, part);
                long long wl = wload_of(f);
                if (viol < best_viol || (viol == best_viol && wl < best_wl)) {
                    best_viol = viol; best_wl = wl; best_f = f;
                }
            }

            if (best_f < 0) {
                // No valid FPGA: traceback
                goto do_traceback;
            }

            if (best_viol > 0) {
                // Has violations: try removing worst FPGAs and retry
                cddt[vi].remove_fpga(best_f, K);
                if (cddt[vi].empty()) goto do_traceback;
                Q.push({cddt[vi].size(), vi});
                continue;
            }

            // Valid assignment found
            part.assign(vi, best_f);
            ++assigned;

            // Propagate constraints: backup ALL neighbors first
            vector<int> affected;
            vector<uint64_t> backups;
            bool conflict = false;
            for (int vk : graph.adj[vi]) {
                if (part.node2fpga[vk] >= 0) continue;
                uint64_t adj_mask = 0;
                for (int af : fg.adj[best_f]) adj_mask |= (1ULL << af);
                adj_mask |= (1ULL << best_f);
                affected.push_back(vk);
                backups.push_back(cddt[vk].bits);
                cddt[vk].intersect_with_bits(adj_mask);
                if (cddt[vk].empty()) conflict = true;
            }

            if (conflict) {
                // Undo: restore ALL neighbors
                part.unassign(vi, best_f);
                --assigned;
                for (size_t i = 0; i < affected.size(); ++i)
                    cddt[affected[i]].bits = backups[i];
                cddt[vi].remove_fpga(best_f, K);
                goto do_traceback;
            }
            continue;
        }

    do_traceback:
        ++tracebacks;
        if (tracebacks > max_traceback) {
            throw runtime_error("Traceback limit " + to_string(max_traceback) + " exceeded!");
        }
        if (!cddt[vi].empty()) {
            Q.push({cddt[vi].size(), vi});
        } else {
            throw runtime_error("Node " + to_string(vi) + " all FPGAs failed!");
        }
    }

    int fv = 0;
    for (int vi = 0; vi < N; ++vi)
        fv += CountTopoViolations(vi, part.node2fpga[vi], coarsest_level, fg, part);
    cout << "  Assigned: " << assigned << " nodes, violations: " << (fv/2)
         << ", tracebacks: " << tracebacks << endl;
    part.violation_cnt = fv / 2;

    for (int i = 0; i < N; ++i)
        if (part.node2fpga[i] < 0)
            throw runtime_error("Unassigned node " + to_string(i));

    cout << "========== Module 4 Done ==========\n" << endl;
    return part;
}

// ================================================================
//  Grow Mode: Region-growing + topology-preserving spread
// ================================================================
static PartitionResult InitialPartitionGrow(CoarsenLevel& coarsest_level,
                                             const FPGAGraph& fg) {
    int N = coarsest_level.super_graph.node_num;
    int K = fg.fpga_num;
    auto& graph = coarsest_level.super_graph;
    auto& weight = coarsest_level.super_weight;

    cout << "\n========== Module 4: InitialPartition (Grow: Region-Growing + Spread) ==========" << endl;
    cout << "  Super-nodes: " << N << ", FPGAs: " << K << endl;

    PartitionResult part;
    part.init(N, K);

    auto wload_of = [&](int f) -> long long {
        long long wl = 0;
        for (int nid : part.fpga2nodes[f]) wl += (nid < N) ? weight[nid] : 1;
        return wl;
    };

    // Step 1: BFS region-growing (topology-first)
    cout << "  Region-growing..." << endl;
    vector<bool> assigned(N, false);
    queue<int> bfs_q;

    for (int vi : graph.fixed_nodes) {
        part.assign(vi, graph.node2fpga.at(vi));
        assigned[vi] = true;
    }

    vector<int> seed_order(N);
    for (int i = 0; i < N; ++i) seed_order[i] = i;
    for (int i = N-1; i > 0; --i) { int j = rand() % (i+1); swap(seed_order[i], seed_order[j]); }

    for (int si = 0; si < N; ++si) {
        int seed = seed_order[si];
        if (assigned[seed]) continue;
        int seed_f = 0; long long min_wl = wload_of(0);
        for (int f = 1; f < K; ++f) { long long wl = wload_of(f); if (wl < min_wl) { min_wl = wl; seed_f = f; } }
        part.assign(seed, seed_f); assigned[seed] = true;
        bfs_q.push(seed);
        while (!bfs_q.empty()) {
            int u = bfs_q.front(); bfs_q.pop();
            for (int v : graph.adj[u]) {
                if (assigned[v]) continue;
                part.assign(v, part.node2fpga[u]); assigned[v] = true;
                bfs_q.push(v);
            }
        }
    }
    for (int vi = 0; vi < N; ++vi) if (!assigned[vi]) {
        int bf = 0; long long bw = wload_of(0);
        for (int f = 1; f < K; ++f) { long long w = wload_of(f); if (w < bw) { bw = w; bf = f; } }
        part.assign(vi, bf);
    }

    // Step 2: Topology-preserving spread (adaptive convergence)
    cout << "  Topology-preserving spread..." << endl;
    int spread_total = 0, stall = 0, pass;
    for (pass = 0; ; ++pass) {
        int pm = 0;
        for (int vi = 0; vi < N; ++vi) {
            int old_f = part.node2fpga[vi];
            if (old_f < 0 || graph.fixed_nodes.count(vi)) continue;
            if (wload_of(old_f) <= fg.resource_cap[old_f]) continue;
            if (CountTopoViolations(vi, old_f, coarsest_level, fg, part) > 0) continue;
            int best_f = -1; long long best_wl = LLONG_MAX;
            for (int tf = 0; tf < K; ++tf) {
                if (tf == old_f) continue;
                if (wload_of(tf) + 1 > fg.resource_cap[tf]) continue;
                if (wload_of(tf) >= wload_of(old_f)) continue;
                if (CountTopoViolations(vi, tf, coarsest_level, fg, part) > 0) continue;
                if (wload_of(tf) < best_wl) { best_wl = wload_of(tf); best_f = tf; }
            }
            if (best_f >= 0) { part.unassign(vi, old_f); part.assign(vi, best_f); ++spread_total; ++pm; }
        }
        if (pm == 0) { ++stall; if (stall >= 5) break; } else stall = 0;
    }
    cout << "  Spread: " << spread_total << " moves (converged after " << (pass+1) << " passes)" << endl;

    int fv = 0, fo = 0;
    for (int vi = 0; vi < N; ++vi) fv += CountTopoViolations(vi, part.node2fpga[vi], coarsest_level, fg, part);
    for (int f = 0; f < K; ++f) if (wload_of(f) > fg.resource_cap[f]) ++fo;
    cout << "  Final: viol=" << (fv/2) << " over=" << fo << endl;
    part.violation_cnt = fv / 2;

    for (int i = 0; i < N; ++i) if (part.node2fpga[i] < 0) throw runtime_error("Unassigned " + to_string(i));
    cout << "========== Module 4 Done ==========\n" << endl;
    return part;
}

