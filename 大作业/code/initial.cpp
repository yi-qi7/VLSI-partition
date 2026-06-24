/**
 * @file    initial.cpp
 * @brief   模块4：初始划分（论文 Algorithm 2：T_vec + Rj + 回溯）
 */

#include "initial.h"
#include "utils.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <queue>
using namespace std;

// ================================================================
long long EvaluateCutDelta(const CoarsenLevel& level,
                            const PartitionResult& part,
                            int vi, int target_f) {
    long long delta = 0;
    int old_f = part.node2fpga[vi];
    for (int nb : level.super_graph.adj[vi]) {
        int nb_f = part.node2fpga[nb];
        if (nb_f < 0) continue;
        bool was_cut = (old_f >= 0) && (old_f != nb_f);
        bool new_cut = (target_f != nb_f);
        if (!was_cut && new_cut)      delta += 1;
        else if (was_cut && !new_cut) delta -= 1;
    }
    return delta;
}

// ================================================================
static int CountTopoViolations(int vi, int target_f,
                                const CoarsenLevel& level,
                                const FPGAGraph& fg,
                                const PartitionResult& part) {
    int violations = 0;
    for (int vk : level.super_graph.adj[vi]) {
        int vk_f = part.node2fpga[vk];
        if (vk_f < 0) continue;
        if (fg.dist_all[target_f][vk_f] > 1) ++violations;
    }
    return violations;
}

static PartitionResult InitialPartitionPaper(CoarsenLevel&, const FPGAGraph&, int);
static PartitionResult InitialPartitionGrow(CoarsenLevel&, const FPGAGraph&);

PartitionResult InitialPartition(CoarsenLevel& coarsest_level,
                                  const FPGAGraph& fg,
                                  int max_traceback,
                                  const string& mode) {
    if (mode == "grow")
        return InitialPartitionGrow(coarsest_level, fg);
    return InitialPartitionPaper(coarsest_level, fg, max_traceback);
}

// ================================================================
//  Paper 模式：论文 Algorithm 2 完整实现
//   Q (候选集大小) + Rj (割边增量) + T_vec (引用计数) + Traceback
// ================================================================
static PartitionResult InitialPartitionPaper(CoarsenLevel& coarsest_level,
                                              const FPGAGraph& fg,
                                              int /*max_traceback*/) {
    int N = coarsest_level.super_graph.node_num;
    int K = fg.fpga_num;
    auto& cddt   = coarsest_level.super_cddt;
    auto& graph  = coarsest_level.super_graph;
    auto& weight = coarsest_level.super_weight;

    cout << "\n========== Module 4: InitialPartition (Algorithm 2) ==========" << endl;
    cout << "  Super-nodes: " << N << ", FPGAs: " << K << endl;

    PartitionResult part;
    part.init(N, K);

    auto wload_of = [&](int f) -> long long {
        long long wl = 0;
        for (int nid : part.fpga2nodes[f])
            wl += (nid < N) ? weight[nid] : 1;
        return wl;
    };

    // ---- Line 1-2: 分配固定节点 ----
    for (int vi : graph.fixed_nodes)
        part.assign(vi, graph.node2fpga.at(vi));

    // ---- Line 5: 构建 Q (候选集大小升序) ----
    typedef pair<int, int> QEntry;
    priority_queue<QEntry, vector<QEntry>, greater<QEntry>> Q;
    for (int vi = 0; vi < N; ++vi) {
        if (part.node2fpga[vi] >= 0) continue;
        if (cddt[vi].empty())
            throw runtime_error("Node " + to_string(vi) + " has empty Cddt!");
        Q.push({cddt[vi].size(), vi});
    }

    int assigned = 0, tracebacks = 0;

    // ---- Line 9: 主循环 ----
    while (!Q.empty()) {
        // Line 10: vj ← Q.pop()
        auto [sz, vj] = Q.top(); Q.pop();
        if (part.node2fpga[vj] >= 0) continue;

        int cur_sz = cddt[vj].size();
        if (cur_sz != (int)sz) {
            if (cur_sz > 0) Q.push({cur_sz, vj});
            continue;
        }
        if (cur_sz == 0) continue;

        // ---- Line 6-8, 11: 构建 Rj (割边增量升序) ----
        vector<pair<long long, int>> Rj;
        {
            uint64_t bits = cddt[vj].bits;
            while (bits) {
                unsigned long f_ul;
                _BitScanForward64(&f_ul, bits);
                int f = (int)f_ul;
                bits &= bits - 1;
                // 动态松弛：10% 弹性为 Traceback 提供腾挪空间
                if (wload_of(f) + weight[vj] > fg.resource_cap[f] * 1.3) continue;
                long long delta = EvaluateCutDelta(coarsest_level, part, vj, f);
                Rj.push_back({delta, f});
            }
        }
        if (Rj.empty()) continue;
        sort(Rj.begin(), Rj.end());

        bool assigned_ok = false;

        for (auto& [cut_delta, v_hat_j] : Rj) {
            // Line 12: part(vj) ← v̂j
            part.assign(vj, v_hat_j);
            ++assigned;

            // ---- Line 14-20: 约束传播 (T_vec 引用计数) ----
            bool conflict = false;
            vector<int> affected;
            // Bug fix: to_remove = 所有不在 keep_mask 中的 K 个 FPGA
            //   (而非仅当前 bits 中的) — 确保多重约束下 T_vec 正确叠加
            uint64_t all_k_mask = (K >= 64) ? ~0ULL : ((1ULL << K) - 1);

            for (int vk : graph.adj[vj]) {
                if (part.node2fpga[vk] >= 0) continue;
                affected.push_back(vk);

                // 保留掩码: {v̂j} ∪ adj(v̂j)
                uint64_t keep_mask = (1ULL << v_hat_j);
                for (int af : fg.adj[v_hat_j])
                    keep_mask |= (1ULL << af);

                // Bug fix: to_remove = 所有 K 个 FPGA 中不在 keep_mask 的
                uint64_t to_remove = all_k_mask & ~keep_mask;

                // Line 15: Tk[î]-- for EACH removed FPGA
                uint64_t tr = to_remove;
                while (tr) {
                    unsigned long f_ul;
                    _BitScanForward64(&f_ul, tr);
                    int f = (int)f_ul;
                    tr &= tr - 1;
                    if (f < K) cddt[vk].T_vec[f]--;
                }

                // 更新 bits: 只保留 T_vec[f] > 0 的 FPGA
                uint64_t new_bits = 0;
                for (int f = 0; f < K; ++f)
                    if (cddt[vk].T_vec[f] > 0)
                        new_bits |= (1ULL << f);
                cddt[vk].bits = new_bits;

                // Line 16-18: |Cddt(vk)| == 0 → break (论文 Line 18)
                if (cddt[vk].bits == 0) {
                    conflict = true;
                    break;  // 检测到冲突立即中断，不再处理后续邻居
                }
            }

            // ---- Line 21-27: Traceback ----
            if (conflict) {
                // Line 22: 从 Cddt(vj) 中永久移除 v̂j
                cddt[vj].remove_fpga(v_hat_j);
                // Line 23: part(vj) ← -1
                part.unassign(vj, v_hat_j);
                --assigned;

                // Line 25-26: 恢复 T_vec (Tk[î]++)
                // to_restore = 所有不在 keep_mask 中的 FPGA（与传播时对称）
                uint64_t keep_mask_r = (1ULL << v_hat_j);
                for (int af : fg.adj[v_hat_j])
                    keep_mask_r |= (1ULL << af);
                uint64_t to_restore = all_k_mask & ~keep_mask_r;

                for (int vk : affected) {
                    uint64_t tr2 = to_restore;
                    while (tr2) {
                        unsigned long f_ul;
                        _BitScanForward64(&f_ul, tr2);
                        int f = (int)f_ul;
                        tr2 &= tr2 - 1;
                        if (f < K) cddt[vk].T_vec[f]++;
                    }

                    // 从 T_vec 重建 bits 并重新入队 Q
                    uint64_t restored = 0;
                    for (int f = 0; f < K; ++f)
                        if (cddt[vk].T_vec[f] > 0)
                            restored |= (1ULL << f);
                    cddt[vk].bits = restored;
                    // Line 27: Add to Q
                    Q.push({cddt[vk].size(), vk});
                }

                ++tracebacks;
            } else {
                // Line 20: 邻居重新入队 Q
                for (int vk : affected)
                    Q.push({cddt[vk].size(), vk});

                // 分配成功，跳出 FPGA 循环
                assigned_ok = true;
                break;
            }
        }

        // Line 24: 若未成功分配，vj 重新入队
        if (!assigned_ok && !cddt[vj].empty())
            Q.push({cddt[vj].size(), vj});
    }

    // ---- 兜底：贪心分配未分配节点 ----
    int greedy = 0;
    for (int vi = 0; vi < N; ++vi) {
        if (part.node2fpga[vi] >= 0) continue;
        int bf = 0;
        long long bw = LLONG_MAX;
        for (int f = 0; f < K; ++f) {
            long long wl = wload_of(f);
            if (wl < bw) { bw = wl; bf = f; }
        }
        part.assign(vi, bf);
        ++assigned;
        ++greedy;
    }

    // ---- 统计 ----
    int fv = 0;
    for (int vi = 0; vi < N; ++vi)
        fv += CountTopoViolations(vi, part.node2fpga[vi], coarsest_level, fg, part);
    cout << "  Assigned: " << assigned
         << ", viol: " << (fv / 2)
         << ", tracebacks: " << tracebacks
         << ", greedy: " << greedy << endl;
    part.violation_cnt = fv / 2;

    for (int i = 0; i < N; ++i)
        if (part.node2fpga[i] < 0)
            throw runtime_error("Unassigned node " + to_string(i));

    cout << "========== Module 4 Done ==========\n" << endl;
    return part;
}

// ================================================================
//  Grow 模式：BFS 区域生长 + 违规感知负载扩散
// ================================================================
static PartitionResult InitialPartitionGrow(CoarsenLevel& coarsest_level,
                                             const FPGAGraph& fg) {
    int N = coarsest_level.super_graph.node_num;
    int K = fg.fpga_num;
    auto& graph  = coarsest_level.super_graph;
    auto& weight = coarsest_level.super_weight;

    cout << "\n========== Module 4: InitialPartition (Grow) ==========" << endl;
    cout << "  Super-nodes: " << N << ", FPGAs: " << K << endl;

    PartitionResult part;
    part.init(N, K);

    auto wload_of = [&](int f) -> long long {
        long long wl = 0;
        for (int nid : part.fpga2nodes[f])
            wl += (nid < N) ? weight[nid] : 1;
        return wl;
    };

    cout << "  Region-growing..." << endl;
    vector<bool> assigned(N, false);
    queue<int> bfs_q;

    for (int vi : graph.fixed_nodes) {
        part.assign(vi, graph.node2fpga.at(vi));
        assigned[vi] = true;
    }
    for (int vi : graph.fixed_nodes) {
        bfs_q.push(vi);
        while (!bfs_q.empty()) {
            int u = bfs_q.front(); bfs_q.pop();
            for (int v : graph.adj[u]) {
                if (assigned[v]) continue;
                part.assign(v, part.node2fpga[u]);
                assigned[v] = true;
                bfs_q.push(v);
            }
        }
    }

    vector<int> seeds;
    for (int i = 0; i < N; ++i)
        if (!assigned[i]) seeds.push_back(i);
    for (int i = (int)seeds.size() - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        swap(seeds[i], seeds[j]);
    }

    for (int seed : seeds) {
        if (assigned[seed]) continue;
        int sf = 0;
        long long mw = wload_of(0);
        for (int f = 1; f < K; ++f) {
            long long w = wload_of(f);
            if (w < mw) { mw = w; sf = f; }
        }
        part.assign(seed, sf);
        assigned[seed] = true;
        bfs_q.push(seed);
        while (!bfs_q.empty()) {
            int u = bfs_q.front(); bfs_q.pop();
            for (int v : graph.adj[u]) {
                if (assigned[v]) continue;
                part.assign(v, part.node2fpga[u]);
                assigned[v] = true;
                bfs_q.push(v);
            }
        }
    }

    for (int vi = 0; vi < N; ++vi) {
        if (assigned[vi]) continue;
        int bf = 0;
        long long bw = wload_of(0);
        for (int f = 1; f < K; ++f) {
            long long w = wload_of(f);
            if (w < bw) { bw = w; bf = f; }
        }
        part.assign(vi, bf);
    }

    cout << "  Violation-aware spread..." << endl;
    int spread_total = 0, stall = 0, pass;
    for (pass = 0; ; ++pass) {
        int pm = 0;
        for (int vi = 0; vi < N; ++vi) {
            int old_f = part.node2fpga[vi];
            if (old_f < 0 || graph.fixed_nodes.count(vi)) continue;
            if (wload_of(old_f) <= fg.resource_cap[old_f]) continue;

            int old_v = CountTopoViolations(vi, old_f, coarsest_level, fg, part);
            int best_f = -1;
            long long best_score = LLONG_MAX;

            for (int tf = 0; tf < K; ++tf) {
                if (tf == old_f) continue;
                if (wload_of(tf) + 1 > fg.resource_cap[tf]) continue;
                if (wload_of(tf) >= wload_of(old_f)) continue;
                int nv = CountTopoViolations(vi, tf, coarsest_level, fg, part);
                if (nv > old_v) continue;
                long long sc = (long long)nv * 1000000LL + wload_of(tf);
                if (sc < best_score) { best_score = sc; best_f = tf; }
            }
            if (best_f >= 0) {
                part.unassign(vi, old_f);
                part.assign(vi, best_f);
                ++spread_total;
                ++pm;
            }
        }
        if (pm == 0) { ++stall; if (stall >= 5) break; }
        else         stall = 0;
    }
    cout << "  Spread: " << spread_total << " (" << (pass + 1) << " passes)" << endl;

    int fv = 0, fo = 0;
    for (int vi = 0; vi < N; ++vi)
        fv += CountTopoViolations(vi, part.node2fpga[vi], coarsest_level, fg, part);
    for (int f = 0; f < K; ++f)
        if (wload_of(f) > fg.resource_cap[f]) ++fo;
    cout << "  Final: viol=" << (fv / 2) << " over=" << fo << endl;
    part.violation_cnt = fv / 2;

    for (int i = 0; i < N; ++i)
        if (part.node2fpga[i] < 0)
            throw runtime_error("Unassigned " + to_string(i));

    cout << "========== Module 4 Done ==========\n" << endl;
    return part;
}
