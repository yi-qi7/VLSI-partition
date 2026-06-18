/**
 * @file    topopart_refine.cpp
 * @brief   TopoPart 模块5：逐层解粗化与约束感知细化 UncoarsenRefine 实现
 *
 * 论文 §IV-D: Uncoarsening and Constraint-Aware Refinement
 *
 * 核心思想：
 *   从最粗层向原始层逐层还原，每层执行边界节点迁移优化（类 FM 算法），
 *   全程受候选集 Cddt 约束 + 资源约束 + 拓扑约束三重保护。
 */

#include "topopart_refine.h"
#include "topopart_utils.h"
#include "topopart_initial.h"  // 复用 EvaluateCutDelta
#include <iostream>
#include <algorithm>
#include <set>
#include <climits>
#include <cstdlib>

using namespace std;

// ================================================================
//  单层细化迭代
// ================================================================

bool RefineSingleLevel(CoarsenLevel& level, const FPGAGraph& fg, int max_iter) {
    int N = level.super_graph.node_num;
    int K = fg.fpga_num;
    auto& part = level.level_part;
    auto& cddt = level.super_cddt;
    auto& graph = level.super_graph;

    bool global_improved = false;

    // 电路距离：对于相邻节点，距离恒为 1（无需全点对矩阵）
    // 拓扑约束检查直接使用 dist=1

    for (int iter = 0; iter < max_iter; ++iter) {
        // ----------------------------------------------------------
        // 步骤1：标记边界节点
        //   边界节点定义：存在至少一个邻居分配在不同 FPGA
        // ----------------------------------------------------------
        vector<bool> is_boundary(N, false);
        int boundary_cnt = 0;

        for (int vi = 0; vi < N; ++vi) {
            int vi_f = part.node2fpga[vi];
            if (vi_f < 0) continue;  // 未分配跳过

            for (int vj : graph.adj[vi]) {
                int vj_f = part.node2fpga[vj];
                if (vj_f >= 0 && vj_f != vi_f) {
                    is_boundary[vi] = true;
                    ++boundary_cnt;
                    break;
                }
            }
        }

        if (boundary_cnt == 0) {
            // 无边界节点 → 已是最优
            break;
        }

        // ----------------------------------------------------------
        // 步骤2：对每个边界节点尝试迁移优化
        // ----------------------------------------------------------
        bool round_improved = false;

        for (int vi = 0; vi < N; ++vi) {
            if (!is_boundary[vi]) continue;

            int old_f = part.node2fpga[vi];
            if (old_f < 0) continue;

            // 固定节点不允许移动
            if (graph.fixed_nodes.count(vi)) continue;

            // ----------------------------------------------------------
            // 步骤2a：遍历候选 FPGA — 拓扑约束硬守卫 + 割边/均衡优化
            //
            // 硬约束：迁移不能增加拓扑违规数
            // 软目标：在满足硬约束的前提下，最小化 (cut_delta - balance_gain*2)
            // ----------------------------------------------------------
            long long best_gain = 1;  // 至少需要 gain < 0 才迁移
            int best_f = -1;

            // 当前拓扑违规数
            int old_violations = 0;
            for (int vk : graph.adj[vi]) {
                int v_hat_k = part.node2fpga[vk];
                if (v_hat_k < 0) continue;
                if (fg.dist_all[old_f][v_hat_k] > 1) ++old_violations;
            }

            int cur_load_old = (int)part.fpga2nodes[old_f].size();

            uint64_t bits = cddt[vi].bits;
            while (bits) {
                int target_f; _BitScanForward64((unsigned long*)&target_f, bits);
                bits &= bits - 1;
                if (target_f == old_f) continue;

                // 资源约束检查
                if ((int)part.fpga2nodes[target_f].size() >= fg.resource_cap[target_f]) {
                    continue;
                }

                // ★ 硬约束：不能增加拓扑违规
                int new_violations = 0;
                for (int vk : graph.adj[vi]) {
                    int v_hat_k = part.node2fpga[vk];
                    if (v_hat_k < 0) continue;
                    if (fg.dist_all[target_f][v_hat_k] > 1) ++new_violations;
                }
                if (new_violations > old_violations) continue;  // 硬拒绝

                // 计算割边变化
                long long cut_delta = EvaluateCutDelta(level, fg, part, vi, target_f);

                // 负载均衡改善（正值 = 负载更均衡）
                int cur_load_target = (int)part.fpga2nodes[target_f].size();
                long long balance_gain = (cur_load_old - 1) - cur_load_target;

                // 拓扑改善（正值 = 违规减少）
                long long topo_improve = (old_violations - new_violations);

                // 复合增益（负值 = 总体改善）
                // 优先减少违规，其次优化割边，最后均衡负载
                long long gain = cut_delta - topo_improve * 100 - balance_gain * 2;

                if (gain < best_gain) {
                    best_gain = gain;
                    best_f = target_f;
                }
            }

            // ----------------------------------------------------------
            // 步骤2b：执行最优迁移（gain < 0 才迁移）
            // ----------------------------------------------------------
            if (best_f >= 0 && best_gain < 0) {
                part.unassign(vi, old_f);
                part.assign(vi, best_f);

                part.cut_size += best_gain;  // 复合 gain（近似）
                round_improved = true;
            }
        }

        if (!round_improved) {
            // 本轮无改进 → 终止该层细化
            break;
        }
        global_improved = true;
    }

    return global_improved;
}

// ================================================================
//  逐层解粗化与细化主函数
// ================================================================

PartitionResult UncoarsenRefine(vector<CoarsenLevel>& levels, const FPGAGraph& fg) {
    int L = (int)levels.size();
    int K = fg.fpga_num;

    cout << "\n========== Module 5: UncoarsenRefine ==========" << endl;
    cout << "  Levels: " << L << endl;

    if (L < 2) {
        // 只有原始层，直接对该层做细化
        cout << "  Single level only, refining directly..." << endl;
        RefineSingleLevel(levels[0], fg, 20);
        return levels[0].level_part;
    }

    // ----------------------------------------------------------
    // 步骤1：在最粗层执行初始划分细化
    // ----------------------------------------------------------
    cout << "  Coarsest level (L=" << L - 1 << ", |V|="
         << levels.back().super_graph.node_num << ") refining..." << endl;
    RefineSingleLevel(levels.back(), fg, 20);

    // ----------------------------------------------------------
    // 步骤2：逐层解粗化（从最粗层向原始层）
    // ----------------------------------------------------------
    for (int lv = L - 2; lv >= 0; --lv) {
        auto& coarse = levels[lv + 1];  // 较粗层
        auto& fine   = levels[lv];      // 较细层（当前目标层）

        int fine_N = fine.super_graph.node_num;

        cout << "  Uncoarsening L=" << lv << " (|V|=" << fine_N
             << "), inheriting from L=" << (lv + 1) << "..." << endl;

        // ----------------------------------------------------------
        // 步骤2a：继承上层划分 → 拆分超节点
        //   超节点中所有子节点继承超节点的 FPGA 分配
        // ----------------------------------------------------------
        fine.level_part.init(fine_N, K);

        for (int sid = 0; sid < coarse.super_graph.node_num; ++sid) {
            int super_fpga = coarse.level_part.node2fpga[sid];
            if (super_fpga < 0) {
                // 超节点未分配（不应出现，但做防御）
                cerr << "  WARNING: super-node " << sid << " unassigned!" << endl;
                continue;
            }

            // 拆分超节点：所有子节点继承 FPGA
            auto it = coarse.super2origin.find(sid);
            if (it == coarse.super2origin.end()) {
                cerr << "  WARNING: super-node " << sid << " has no child mapping!" << endl;
                continue;
            }

            for (int child : it->second) {
                fine.level_part.assign(child, super_fpga);
                // 应用固定节点约束（覆盖继承值）
                if (fine.super_graph.fixed_nodes.count(child)) {
                    int fixed_f = fine.super_graph.node2fpga.at(child);
                    if (fixed_f != super_fpga) {
                        fine.level_part.unassign(child, super_fpga);
                        fine.level_part.assign(child, fixed_f);
                    }
                }
            }
        }

        // Refine + KL-FM at this level
        cout << "    Refining L=" << lv << "..." << endl;
        RefineSingleLevel(fine, fg, 10);
        // Quick KL-FM at each level for cut optimization
        if (lv <= 3) {  // Only at finer levels where it matters
            int lvl_fm = 0;
            for (int fp = 0; fp < 5; ++fp) {
                int pfm = 0;
                for (int vi = 0; vi < fine_N; ++vi) {
                    int old_f = fine.level_part.node2fpga[vi];
                    if (old_f < 0 || fine.super_graph.fixed_nodes.count(vi)) continue;
                    int old_v = 0;
                    for (int vk : fine.super_graph.adj[vi]) {
                        int vk_f = fine.level_part.node2fpga[vk];
                        if (vk_f >= 0 && fg.dist_all[old_f][vk_f] > 1) ++old_v;
                    }
                    if (old_v > 0) continue;
                    long long best_cd = 0; int best_tf = -1;
                    for (int tf = 0; tf < K; ++tf) {
                        if (tf == old_f) continue;
                        if ((int)fine.level_part.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                        int nv = 0;
                        for (int vk : fine.super_graph.adj[vi]) { int vk_f = fine.level_part.node2fpga[vk]; if (vk_f >= 0 && fg.dist_all[tf][vk_f] > 1) ++nv; }
                        if (nv > 0) continue;
                        long long cd = EvaluateCutDelta(fine, fg, fine.level_part, vi, tf);
                        if (cd < best_cd) { best_cd = cd; best_tf = tf; }
                    }
                    if (best_tf >= 0) {
                        fine.level_part.unassign(vi, old_f);
                        fine.level_part.assign(vi, best_tf);
                        ++lvl_fm; ++pfm;
                    }
                }
                if (pfm == 0) break;
            }
            if (lvl_fm > 0) cout << "      KL-FM moves: " << lvl_fm << endl;
        }
    }

    // ============================================================
    // L0: Alternating balance + KL-FM cycles (3 rounds)
    // ============================================================
    {
        auto& lv0 = levels[0];
        int N0 = lv0.super_graph.node_num;
        int K0 = fg.fpga_num;
        auto& part0 = lv0.level_part;
        auto& graph0 = lv0.super_graph;
        long long avg_load = N0 / K0;

        for (int round = 0; round < 3; ++round) {
            cout << "  L0 Round " << round << "..." << endl;

            // Phase A: Topology-preserving balance
            int bal_total = 0, bal_stall = 0;
            for (int pass = 0; ; ++pass) {
                int pass_bal = 0;
                for (int vi = 0; vi < N0; ++vi) {
                    int old_f = part0.node2fpga[vi];
                    if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                    if ((int)part0.fpga2nodes[old_f].size() <= fg.resource_cap[old_f]) continue;
                    int old_v = 0;
                    for (int vk : graph0.adj[vi]) { int vk_f = part0.node2fpga[vk]; if (vk_f >= 0 && fg.dist_all[old_f][vk_f] > 1) ++old_v; }
                    if (old_v > 0) continue;
                    int best_f = -1; long long best_load = LLONG_MAX;
                    for (int tf = 0; tf < K0; ++tf) {
                        if (tf == old_f) continue;
                        long long tl = (int)part0.fpga2nodes[tf].size();
                        if (tl + 1 > fg.resource_cap[tf] || tl >= (int)part0.fpga2nodes[old_f].size()) continue;
                        int nv = 0;
                        for (int vk : graph0.adj[vi]) { int vk_f = part0.node2fpga[vk]; if (vk_f >= 0 && fg.dist_all[tf][vk_f] > 1) ++nv; }
                        if (nv > 0) continue;
                        if (tl < best_load) { best_load = tl; best_f = tf; }
                    }
                    if (best_f >= 0) { part0.unassign(vi, old_f); part0.assign(vi, best_f); ++bal_total; ++pass_bal; }
                }
                if (pass_bal == 0) { ++bal_stall; if (bal_stall >= 3) break; } else bal_stall = 0;
            }

            // Phase B: KL-FM cut optimization
            int fm_total = 0;
        for (int pass = 0; ; ++pass) {  // adaptive: breaks when kept==0
            // Compute gains for all boundary nodes
            struct Move { long long gain; int node; int target; };
            vector<Move> candidates;
            for (int vi = 0; vi < N0; ++vi) {
                int old_f = part0.node2fpga[vi];
                if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                int old_v = 0;
                for (int vk : graph0.adj[vi]) {
                    int vk_f = part0.node2fpga[vk];
                    if (vk_f >= 0 && fg.dist_all[old_f][vk_f] > 1) ++old_v;
                }
                if (old_v > 0) continue;
                for (int tf = 0; tf < K0; ++tf) {
                    if (tf == old_f) continue;
                    if ((int)part0.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                    int nv = 0;
                    for (int vk : graph0.adj[vi]) {
                        int vk_f = part0.node2fpga[vk];
                        if (vk_f >= 0 && fg.dist_all[tf][vk_f] > 1) ++nv;
                    }
                    if (nv > 0) continue;
                    long long g = -EvaluateCutDelta(lv0, fg, part0, vi, tf);
                    if (g > 0) candidates.push_back({g, vi, tf});
                }
            }

            if (candidates.empty()) break;

            // Sort by gain descending
            sort(candidates.begin(), candidates.end(),
                 [](const Move& a, const Move& b) { return a.gain > b.gain; });

            // KL pass: make moves, track cumulative gain
            vector<bool> locked(N0, false);
            vector<int> moves_vi, moves_tf;
            vector<long long> gains;
            long long cumulative = 0;
            long long best_cum = 0;
            int best_idx = -1;

            for (auto& mv : candidates) {
                int vi = mv.node, tf = mv.target;
                if (locked[vi]) continue;
                int old_f = part0.node2fpga[vi];
                if (old_f < 0) continue;

                // Re-verify (neighbors might have been moved)
                int old_v = 0;
                for (int vk : graph0.adj[vi]) {
                    int vk_f = part0.node2fpga[vk];
                    if (vk_f >= 0 && fg.dist_all[old_f][vk_f] > 1) ++old_v;
                }
                if (old_v > 0) continue;
                if ((int)part0.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                int nv = 0;
                for (int vk : graph0.adj[vi]) {
                    int vk_f = part0.node2fpga[vk];
                    if (vk_f >= 0 && fg.dist_all[tf][vk_f] > 1) ++nv;
                }
                if (nv > 0) continue;

                long long g = -EvaluateCutDelta(lv0, fg, part0, vi, tf);
                // Accept move even if gain < 0 (KL allows temporary worsening)
                part0.unassign(vi, old_f);
                part0.assign(vi, tf);
                locked[vi] = true;
                cumulative += g;
                moves_vi.push_back(vi);
                moves_tf.push_back(old_f);
                gains.push_back(g);

                if (cumulative > best_cum) {
                    best_cum = cumulative;
                    best_idx = (int)moves_vi.size();
                }
            }

            // Rollback moves after best prefix
            int kept = (best_idx > 0) ? best_idx : 0;
            for (int i = (int)moves_vi.size() - 1; i >= kept; --i) {
                part0.unassign(moves_vi[i], part0.node2fpga[moves_vi[i]]);
                part0.assign(moves_vi[i], moves_tf[i]);
            }

            fm_total += kept;
            cout << "    KL Pass " << pass << ": tried=" << moves_vi.size()
                 << " kept=" << kept << " best_cum=" << best_cum << endl;
            if (kept == 0) break;
        }
        cout << "  L0 KL-FM total moves: " << fm_total << endl;

        // Smart re-balance: move from overloaded to underloaded (topology-preserving)
        cout << "  Smart re-balance..." << endl;
        int rebal = 0;
        for (int pass = 0; pass < 20; ++pass) {
            int pr = 0;
            for (int vi = 0; vi < N0; ++vi) {
                int old_f = part0.node2fpga[vi];
                if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                int src_sz = (int)part0.fpga2nodes[old_f].size();
                if (src_sz <= avg_load) continue;
                int old_v = 0;
                for (int vk : graph0.adj[vi]) {
                    int vk_f = part0.node2fpga[vk];
                    if (vk_f >= 0 && fg.dist_all[old_f][vk_f] > 1) ++old_v;
                }
                if (old_v > 0) continue;
                int best_f = -1; int best_sz = INT_MAX;
                for (int tf = 0; tf < K0; ++tf) {
                    if (tf == old_f) continue;
                    int tgt_sz = (int)part0.fpga2nodes[tf].size();
                    if (tgt_sz >= src_sz) continue;
                    if (tgt_sz + 1 > fg.resource_cap[tf]) continue;
                    int nv = 0;
                    for (int vk : graph0.adj[vi]) {
                        int vk_f = part0.node2fpga[vk];
                        if (vk_f >= 0 && fg.dist_all[tf][vk_f] > 1) ++nv;
                    }
                    if (nv > 0) continue;
                    long long cd = EvaluateCutDelta(lv0, fg, part0, vi, tf);
                    // Accept small cut increase for better balance
                    if (cd <= 2 && tgt_sz < best_sz) { best_sz = tgt_sz; best_f = tf; }
                }
                if (best_f >= 0) {
                    part0.unassign(vi, old_f);
                    part0.assign(vi, best_f);
                    ++rebal; ++pr;
                }
            }
            if (pr == 0) break;
        }
        cout << "  Re-balance moves: " << rebal << endl;

        // Final quick KL-FM pass to clean up
        for (int pass = 0; pass < 5; ++pass) {
            int pf = 0;
            for (int vi = 0; vi < N0; ++vi) {
                int old_f = part0.node2fpga[vi];
                if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                int old_v = 0;
                for (int vk : graph0.adj[vi]) { int vk_f = part0.node2fpga[vk]; if (vk_f >= 0 && fg.dist_all[old_f][vk_f] > 1) ++old_v; }
                if (old_v > 0) continue;
                long long best_cd = 0; int best_tf = -1;
                for (int tf = 0; tf < K0; ++tf) {
                    if (tf == old_f) continue;
                    if ((int)part0.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                    int nv = 0;
                    for (int vk : graph0.adj[vi]) { int vk_f = part0.node2fpga[vk]; if (vk_f >= 0 && fg.dist_all[tf][vk_f] > 1) ++nv; }
                    if (nv > 0) continue;
                    long long cd = EvaluateCutDelta(lv0, fg, part0, vi, tf);
                    if (cd < best_cd) { best_cd = cd; best_tf = tf; }
                }
                if (best_tf >= 0) { part0.unassign(vi, old_f); part0.assign(vi, best_tf); ++pf; }
            }
            if (pf == 0) break;
        }
        }  // end round loop
    }

    // ----------------------------------------------------------
    // 步骤4：返回原始层划分结果
    // ----------------------------------------------------------
    auto& result = levels[0].level_part;

    // 计算最终割边数
    result.cut_size = 0;
    for (int vi = 0; vi < levels[0].super_graph.node_num; ++vi) {
        int vi_f = result.node2fpga[vi];
        for (int vj : levels[0].super_graph.adj[vi]) {
            int vj_f = result.node2fpga[vj];
            if (vi_f >= 0 && vj_f >= 0 && vi_f != vj_f) {
                ++result.cut_size;
            }
        }
    }
    result.cut_size /= 2;  // 每条割边被计算两次

    cout << "  Final cut size (original level): " << result.cut_size << endl;
    cout << "========== Module 5 Done ==========\n" << endl;

    return result;
}
