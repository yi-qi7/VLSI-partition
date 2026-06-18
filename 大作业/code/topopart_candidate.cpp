/**
 * @file    topopart_candidate.cpp
 * @brief   TopoPart 模块2：候选 FPGA 传播实现（Algorithm 1）
 *
 * 论文 Algorithm 1: Candidate FPGA Set Propagation
 *
 * 定理 III.1 判别依据：
 *   vi → v̂i, vj → v̂j 无拓扑违规 ⇔ dist_circuit(vi, vj) ≥ dist_fpga(v̂i, v̂j)
 *
 * 传播收缩逻辑：
 *   已知 vi 固定绑定 v̂i，对于可移动节点 vj（电路距离 k = dist(vi, vj)）：
 *     vj 只能分配在 Ŝ(v̂i, k) 中（距离 v̂i ≤ k 的 FPGA），
 *     否则：v̂j ∉ Ŝ(v̂i, k) → dist_fpga(v̂i, v̂j) > k = dist_circuit(vi, vj)
 *     → 违反定理 III.1 的 x ≥ y 条件 → 拓扑违规。
 */

#include "topopart_candidate.h"
#include "topopart_utils.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>

vector<CandidateSet> CandidateFPGAPropagation(CircuitGraph& g, FPGAGraph& fg) {
    int N = g.node_num;
    int K = fg.fpga_num;

    cout << "\n========== Module 2: CandidateFPGAPropagation ==========" << endl;
    cout << "  Circuit nodes: " << N << ", FPGAs: " << K << endl;
    cout << "  Fixed nodes: " << g.fixed_nodes.size() << endl;

    // ----------------------------------------------------------
    // 步骤1：前置条件检查（调用者应已完成）
    //   - fg.dist_all 已计算
    //   - fg.max_dist 已计算
    //   - fg.S_cache 已构建
    //   - g.adj 已构建（按需 BFS，无需 g.dist_all）
    // ----------------------------------------------------------
    if (fg.S_cache.empty()) {
        throw runtime_error("CandidateFPGAPropagation: FPGA S_cache not initialized!"
                            " Please call BuildFPGADist(fg) first.");
    }
    // g.dist_all 不再使用（已改用按需 BFS 节省内存），无需检查

    // ----------------------------------------------------------
    // 步骤2：初始化候选集
    //   - 固定节点：仅保留绑定的那个 FPGA
    //   - 可移动节点：候选集 = 全部 FPGA
    // ----------------------------------------------------------
    g.init_move_nodes();

    vector<CandidateSet> cddt(N, CandidateSet(K));
    for (int i = 0; i < N; ++i) {
        if (g.fixed_nodes.count(i)) {
            // 固定节点：候选集仅包含绑定 FPGA
            int bound_fpga = g.node2fpga.at(i);
            cddt[i].bits |= (1ULL << bound_fpga);
            cddt[i].T_vec[bound_fpga] = 1;
        } else {
            // 可移动节点：候选集 = 全部 FPGA
            cddt[i].reset_all(K);
        }
    }

    // ----------------------------------------------------------
    // 步骤3：初始化队列 q
    //   队列元素：(固定节点 id, 绑定 FPGA id)
    // ----------------------------------------------------------
    queue<pair<int, int>> q;
    for (int vi : g.fixed_nodes) {
        int v_hat_i = g.node2fpga.at(vi);
        q.push({vi, v_hat_i});
    }

    // 记录哪些可移动节点已转为固定节点（避免重复入队）
    // fixed_nodes 集合已包含所有固定节点（包括传播中新固定的）
    // 我们通过检查 cddt[vj].is_singleton() 来判断

    // ----------------------------------------------------------
    // 步骤4：主传播循环（BFS 风格）
    // ----------------------------------------------------------
    int propagation_round = 0;
    int new_fixed_count = 0;

    while (!q.empty()) {
        ++propagation_round;
        auto [vi, v_hat_i] = q.front();
        q.pop();

        // a. d = maxDist(v̂i)：FPGA v̂i 到所有 FPGA 的最大最短距离
        int d = fg.max_dist[v_hat_i];

        // ★ 内存优化：不再使用预计算的 N×N 距离矩阵（N=300K 时需 ~360GB）
        // 改用按需 BFS，以 vi 为源点，深度 < d
        // 定理 III.1：对于 BFS 发现的每个节点 vj（距离 k < d），
        //   收缩候选集 Cddt[vj] = Cddt[vj] ∩ Ŝ(v̂i, k)
        BFSFromSource(g.adj, vi, d - 1,
            [&](int vj, int k) {
                // 跳过源点自身
                if (vj == vi) return;
                // 跳过固定节点
                if (g.fixed_nodes.count(vj)) return;
                // 跳过已经是单候选的节点（已在队列中或已处理）
                if (cddt[vj].is_singleton()) return;

                // b. 取出 Ŝ(v̂i, k)：FPGA 距离 ≤ k 的 FPGA 集合
                const unordered_set<int>& S_hat = fg.get_S(v_hat_i, k);

                // c. 收缩候选集：Cddt[vj] = Cddt[vj] ∩ Ŝ(v̂i, k)
                bool changed = cddt[vj].intersect_with_check(S_hat, K);

                // d. 分支判定
                if (cddt[vj].empty()) {
                    // 候选集为空 → 无可行拓扑划分解
                    throw runtime_error(
                        "CandidateFPGAPropagation: node " + to_string(vj) +
                        " candidate set is empty! No feasible topology partition.\n"
                        "  Propagation source: fixed node " + to_string(vi) +
                        " -> FPGA " + to_string(v_hat_i) +
                        ", circuit distance k = " + to_string(k));
                }

                if (cddt[vj].is_singleton() && changed) {
                    // 候选集收缩为单元素 → vj 转为新固定节点
                    int only_fpga = cddt[vj].get_only_fpga();
                    g.fixed_nodes.insert(vj);
                    g.node2fpga[vj] = only_fpga;
                    g.move_nodes.erase(vj);
                    q.push({vj, only_fpga});
                    ++new_fixed_count;
                }
            });
    }

    // ----------------------------------------------------------
    // 步骤5：输出统计
    // ----------------------------------------------------------
    cout << "  Propagation rounds: " << propagation_round << endl;
    cout << "  Newly fixed nodes: " << new_fixed_count << endl;
    cout << "  Total fixed nodes: " << g.fixed_nodes.size() << endl;

    // 统计候选集大小分布
    vector<int> cddt_size_dist(K + 1, 0);
    for (int i = 0; i < N; ++i) {
        int sz = cddt[i].size();
        if (sz <= K) cddt_size_dist[sz]++;
    }
    cout << "  Candidate set size distribution:" << endl;
    for (int sz = 1; sz <= min(K, 10); ++sz) {
        if (cddt_size_dist[sz] > 0)
            cout << "    |Cddt|=" << sz << ": " << cddt_size_dist[sz] << " nodes" << endl;
    }
    if (cddt_size_dist[0] > 0)
        cout << "    WARNING |Cddt|=0 (empty): " << cddt_size_dist[0] << " nodes!!!" << endl;

    cout << "========== Module 2 Done ==========\n" << endl;
    return cddt;
}
