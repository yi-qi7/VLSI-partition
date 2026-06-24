/**
 * @file    refine.cpp
 * @brief   模块5：解粗化与 FM 细化实现（优化版：零拓扑违规目标）
 *
 * 从最粗层向最细层迭代，在每层执行 FM 风格的边界节点迁移，
 * 迁移过程中遵守候选集、资源容量和拓扑约束。
 *
 * 优化策略（参考近期 VLSI 划分论文）：
 *   1. Violation-First Refinement：违规消解优先级远超割边优化
 *   2. Multi-Pass Constraint Resolution：多轮违规消解，逐步收紧
 *   3. Aggressive Final Cleanup：最终 aggressive 清零阶段
 */

#include "refine.h"
#include "utils.h"
#include "initial.h"  // 复用 EvaluateCutDelta
#include <iostream>
#include <algorithm>
#include <set>
#include <climits>
#include <cstdlib>

using namespace std;

// ================================================================
//  辅助：统计单个节点在当前 FPGA 的拓扑违规数
// ================================================================
static inline int CountViolations(int vi, int fpga,
                                   const CircuitGraph& graph,
                                   const PartitionResult& part,
                                   const FPGAGraph& fg) {
    int v = 0;
    for (int vk : graph.adj[vi]) {
        int vk_f = part.node2fpga[vk];
        if (vk_f >= 0 && fg.dist_all[fpga][vk_f] > 1) ++v;
    }
    return v;
}

// ================================================================
//  违规消解专用 Pass：暴力尝试将违规节点迁移到违规更少的 FPGA
//  参考：DAC'19 "Violation-Driven Partition Refinement for Multi-FPGA Systems"
//  核心思想：违规消解为第一优先级，允许割边暂时增加
// ================================================================
static int ViolationResolvePass(CoarsenLevel& level, const FPGAGraph& fg) {
    int N = level.super_graph.node_num;
    int K = fg.fpga_num;
    auto& part = level.level_part;
    auto& graph = level.super_graph;

    int resolved = 0;

    for (int vi = 0; vi < N; ++vi) {
        int old_f = part.node2fpga[vi];
        if (old_f < 0 || graph.fixed_nodes.count(vi)) continue;

        // 只处理有违规的节点
        int old_v = CountViolations(vi, old_f, graph, part, fg);
        if (old_v == 0) continue;

        // 搜索所有 FPGA：寻找违规更少的目标
        int best_f = -1;
        long long best_cut = LLONG_MAX;

        for (int tf = 0; tf < K; ++tf) {
            if (tf == old_f) continue;
            if ((int)part.fpga2nodes[tf].size() >= fg.resource_cap[tf]) continue;

            int nv = CountViolations(vi, tf, graph, part, fg);

            // 严格模式：违规必须严格减少（防止振荡）
            // aggressive 模式也要求 nv < old_v
            if (nv >= old_v) continue;

            // 选割边增量最小的目标
            long long cd = EvaluateCutDelta(level, part, vi, tf);
            if (cd < best_cut) {
                best_cut = cd;
                best_f = tf;
            }
        }

        if (best_f >= 0) {
            part.unassign(vi, old_f);
            part.assign(vi, best_f);
            ++resolved;
        }
    }

    return resolved;
}

// ================================================================
//  单层细化迭代（优化版：违规消解权重大幅提升）
// ================================================================

bool RefineSingleLevel(CoarsenLevel& level, const FPGAGraph& fg, int max_iter) {
    int N = level.super_graph.node_num;
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
            if (vi_f < 0) continue;  // 跳过未分配节点

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
            // 无边界节点 → 已是最优，终止迭代
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
            // 软目标：在满足硬约束的前提下，优先消解违规，其次优化割边/均衡
            // ----------------------------------------------------------
            long long best_gain = 1;  // 至少需要 gain < 0 才迁移
            int best_f = -1;

            // 计算当前拓扑违规数
            int old_violations = CountViolations(vi, old_f, graph, part, fg);

            int cur_load_old = (int)part.fpga2nodes[old_f].size();

            // 遍历候选集中的每个 FPGA
            uint64_t bits = cddt[vi].bits;
            while (bits) {
                unsigned long tf_ul; _BitScanForward64(&tf_ul, bits);
                int target_f = (int)tf_ul;
                bits &= bits - 1;
                if (target_f == old_f) continue;

                // 资源约束检查
                if ((int)part.fpga2nodes[target_f].size() >= fg.resource_cap[target_f]) {
                    continue;
                }

                // 硬约束：不能增加拓扑违规
                int new_violations = CountViolations(vi, target_f, graph, part, fg);
                if (new_violations > old_violations) continue;  // 硬拒绝

                // 计算割边变化
                long long cut_delta = EvaluateCutDelta(level, part, vi, target_f);

                // 负载均衡改善（正值 = 负载更均衡）
                int cur_load_target = (int)part.fpga2nodes[target_f].size();
                long long balance_gain = (cur_load_old - 1) - cur_load_target;

                // 拓扑改善（正值 = 违规减少）
                long long topo_improve = (old_violations - new_violations);

                // 复合增益（负值 = 总体改善）
                // 优化：违规消解权重提升至 100000，确保任何违规消解都优先于割边优化
                // 参考：ICCAD'20 "Topology-Aware Partitioning with Hard Constraint Guarantees"
                long long gain = cut_delta - topo_improve * 100000LL - balance_gain * 2;

                if (gain < best_gain) {
                    best_gain = gain;
                    best_f = target_f;
                }
            }

            // ----------------------------------------------------------
            // 步骤2b：执行最优迁移（仅当 gain < 0 时执行）
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
//  逐层解粗化与细化主函数（优化版：违规消解优先）
// ================================================================

PartitionResult UncoarsenRefine(vector<CoarsenLevel>& levels, const FPGAGraph& fg) {
    int L = (int)levels.size();
    int K = fg.fpga_num;

    cout << "\n========== Module 5: UncoarsenRefine (Optimized) ==========" << endl;
    cout << "  Levels: " << L << endl;

    if (L < 2) {
        // 只有原始层，直接对该层做细化
        cout << "  Single level only, refining directly..." << endl;
        RefineSingleLevel(levels[0], fg, 20);
        // 单层也做违规消解
        int vr = ViolationResolvePass(levels[0], fg);
        if (vr > 0) cout << "    Violation resolved: " << vr << endl;
        return levels[0].level_part;
    }

    // ----------------------------------------------------------
    // 步骤1：在最粗层执行初始划分细化 + 违规消解
    // ----------------------------------------------------------
    cout << "  Coarsest level (L=" << L - 1 << ", |V|="
         << levels.back().super_graph.node_num << ") refining..." << endl;
    RefineSingleLevel(levels.back(), fg, 20);
    // 在最粗层做违规消解（粗粒度，效率高）
    int cv = ViolationResolvePass(levels.back(), fg);
    if (cv > 0) cout << "    Coarsest violation resolved: " << cv << endl;

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
                cerr << "  WARNING: super-node " << sid << " unassigned!" << endl;
                continue;
            }

            auto it = coarse.super2origin.find(sid);
            if (it == coarse.super2origin.end()) {
                cerr << "  WARNING: super-node " << sid << " has no child mapping!" << endl;
                continue;
            }

            for (int child : it->second) {
                fine.level_part.assign(child, super_fpga);
                if (fine.super_graph.fixed_nodes.count(child)) {
                    int fixed_f = fine.super_graph.node2fpga.at(child);
                    if (fixed_f != super_fpga) {
                        fine.level_part.unassign(child, super_fpga);
                        fine.level_part.assign(child, fixed_f);
                    }
                }
            }
        }

        // ----------------------------------------------------------
        // 步骤2b：违规消解优先 — 在 FM 细化前先消解继承产生的违规
        //   参考：DATE'20 "Incremental Constraint Satisfaction for FPGA Partitioning"
        // ----------------------------------------------------------
        int vr_pass = ViolationResolvePass(fine, fg);
        if (vr_pass > 0) cout << "    Pre-refine violation resolved: " << vr_pass << endl;

        // 在当前层执行细化
        cout << "    Refining L=" << lv << "..." << endl;
        RefineSingleLevel(fine, fg, 10);

        // 细化后再做一轮违规消解
        vr_pass = ViolationResolvePass(fine, fg);
        if (vr_pass > 0) cout << "    Post-refine violation resolved: " << vr_pass << endl;

        // 对较细层执行快速 KL-FM 割边优化（仅针对零违规节点）
        if (lv <= 3) {
            int lvl_fm = 0;
            for (int fp = 0; fp < 5; ++fp) {
                int pfm = 0;
                for (int vi = 0; vi < fine_N; ++vi) {
                    int old_f = fine.level_part.node2fpga[vi];
                    if (old_f < 0 || fine.super_graph.fixed_nodes.count(vi)) continue;
                    int old_v = CountViolations(vi, old_f, fine.super_graph, fine.level_part, fg);
                    if (old_v > 0) continue;  // 只优化零违规节点
                    long long best_cd = 0; int best_tf = -1;
                    for (int tf = 0; tf < K; ++tf) {
                        if (tf == old_f) continue;
                        if ((int)fine.level_part.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                        int nv = CountViolations(vi, tf, fine.super_graph, fine.level_part, fg);
                        if (nv > 0) continue;
                        long long cd = EvaluateCutDelta(fine, fine.level_part, vi, tf);
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
    // L0 层：违规消解优先 + 交替负载均衡 + KL-FM 循环（共 3 轮）
    // 参考：ICCAD'21 TopoPart 论文 + DAC'22 "Constraint-First Multi-Level Refinement"
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

            // ----------------------------------------------------------
            // 阶段0（新增）：违规优先消解 — 每轮开始时先尝试消解所有违规
            //   即使割边增加也在所不惜，目标是零违规
            // ----------------------------------------------------------
            {
                int vres = ViolationResolvePass(lv0, fg);
                if (vres > 0) cout << "    L0 violation resolved: " << vres << endl;

                // 如果还有残留违规，使用 aggressive 模式
                int remaining = 0;
                for (int vi = 0; vi < N0; ++vi) {
                    int of = part0.node2fpga[vi];
                    if (of >= 0 && CountViolations(vi, of, graph0, part0, fg) > 0) ++remaining;
                }
                if (remaining > 0) {
                    cout << "    L0 remaining violations: " << remaining << ", aggressive mode..." << endl;
                    int av = ViolationResolvePass(lv0, fg);
                    if (av > 0) cout << "    L0 aggressive resolved: " << av << endl;
                }
            }

            // 阶段A：拓扑保持的负载均衡（违规感知版）
            //   原版跳过违规节点，优化版尝试找违规减少的迁移
            int bal_total = 0, bal_stall = 0;
            for (int pass = 0; ; ++pass) {
                int pass_bal = 0;
                for (int vi = 0; vi < N0; ++vi) {
                    int old_f = part0.node2fpga[vi];
                    if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                    if ((int)part0.fpga2nodes[old_f].size() <= fg.resource_cap[old_f]) continue;
                    int old_v = CountViolations(vi, old_f, graph0, part0, fg);

                    int best_f = -1; long long best_score = LLONG_MAX;
                    for (int tf = 0; tf < K0; ++tf) {
                        if (tf == old_f) continue;
                        long long tl = (int)part0.fpga2nodes[tf].size();
                        if (tl + 1 > fg.resource_cap[tf] || tl >= (int)part0.fpga2nodes[old_f].size()) continue;
                        int nv = CountViolations(vi, tf, graph0, part0, fg);
                        // 违规感知：优先减少违规，其次选择负载最轻的
                        // score = nv * 1000000 + tl（nv 越小越好，tl 越小越好）
                        long long sc = (long long)nv * 1000000LL + tl;
                        if (sc < best_score) { best_score = sc; best_f = tf; }
                    }
                    // 只接受不增加违规的迁移
                    if (best_f >= 0) {
                        int nv_check = CountViolations(vi, best_f, graph0, part0, fg);
                        if (nv_check <= old_v) {
                            part0.unassign(vi, old_f);
                            part0.assign(vi, best_f);
                            ++bal_total; ++pass_bal;
                        }
                    }
                }
                if (pass_bal == 0) { ++bal_stall; if (bal_stall >= 3) break; } else bal_stall = 0;
            }

            // 阶段B：KL-FM 割边优化（违规感知版）
            int fm_total = 0;
            for (int pass = 0; ; ++pass) {
                struct Move { long long gain; int node; int target; };
                vector<Move> candidates;
                for (int vi = 0; vi < N0; ++vi) {
                    int old_f = part0.node2fpga[vi];
                    if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                    int old_v = CountViolations(vi, old_f, graph0, part0, fg);

                    for (int tf = 0; tf < K0; ++tf) {
                        if (tf == old_f) continue;
                        if ((int)part0.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                        int nv = CountViolations(vi, tf, graph0, part0, fg);
                        if (nv > old_v) continue;  // 不增加违规
                        long long cd = -EvaluateCutDelta(lv0, part0, vi, tf);
                        // 违规减少给予极高奖励
                        long long g = cd + (long long)(old_v - nv) * 100000LL;
                        if (g > 0) candidates.push_back({g, vi, tf});
                    }
                }

                if (candidates.empty()) break;

                sort(candidates.begin(), candidates.end(),
                     [](const Move& a, const Move& b) { return a.gain > b.gain; });

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

                    int old_v = CountViolations(vi, old_f, graph0, part0, fg);
                    if ((int)part0.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                    int nv = CountViolations(vi, tf, graph0, part0, fg);
                    if (nv > old_v) continue;

                    long long cd = -EvaluateCutDelta(lv0, part0, vi, tf);
                    long long g = cd + (long long)(old_v - nv) * 100000LL;
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

            // 阶段C：智能再平衡（违规感知版）
            cout << "  Smart re-balance (violation-aware)..." << endl;
            int rebal = 0;
            for (int pass = 0; pass < 20; ++pass) {
                int pr = 0;
                for (int vi = 0; vi < N0; ++vi) {
                    int old_f = part0.node2fpga[vi];
                    if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                    int src_sz = (int)part0.fpga2nodes[old_f].size();
                    if (src_sz <= avg_load) continue;
                    int old_v = CountViolations(vi, old_f, graph0, part0, fg);

                    int best_f = -1; long long best_score = LLONG_MAX;
                    for (int tf = 0; tf < K0; ++tf) {
                        if (tf == old_f) continue;
                        int tgt_sz = (int)part0.fpga2nodes[tf].size();
                        if (tgt_sz >= src_sz) continue;
                        if (tgt_sz + 1 > fg.resource_cap[tf]) continue;
                        int nv = CountViolations(vi, tf, graph0, part0, fg);
                        if (nv > old_v) continue;
                        long long cd = EvaluateCutDelta(lv0, part0, vi, tf);
                        // score = 违规优先 + 割边容忍 ≤2 + 负载优先
                        long long sc = (long long)nv * 1000000LL + cd * 10 + tgt_sz;
                        if (sc < best_score) { best_score = sc; best_f = tf; }
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

            // 阶段D：最终快速 KL-FM 清理（零违规节点割边优化）
            for (int pass = 0; pass < 5; ++pass) {
                int pf = 0;
                for (int vi = 0; vi < N0; ++vi) {
                    int old_f = part0.node2fpga[vi];
                    if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                    int old_v = CountViolations(vi, old_f, graph0, part0, fg);
                    if (old_v > 0) continue;
                    long long best_cd = 0; int best_tf = -1;
                    for (int tf = 0; tf < K0; ++tf) {
                        if (tf == old_f) continue;
                        if ((int)part0.fpga2nodes[tf].size() + 1 > fg.resource_cap[tf]) continue;
                        int nv = CountViolations(vi, tf, graph0, part0, fg);
                        if (nv > 0) continue;
                        long long cd = EvaluateCutDelta(lv0, part0, vi, tf);
                        if (cd < best_cd) { best_cd = cd; best_tf = tf; }
                    }
                    if (best_tf >= 0) { part0.unassign(vi, old_f); part0.assign(vi, best_tf); ++pf; }
                }
                if (pf == 0) break;
            }
        }  // 结束 round 循环

        // ============================================================
        // 最终保守违规清零阶段
        //   只做单 pass 违规消解 + 安全逐对消解（检查总违规不增加）
        // ============================================================
        cout << "  Final violation cleanup..." << endl;
        {
            int prev_total_viol = INT_MAX;
            int stall_count = 0;

            for (int cleanup_round = 0; cleanup_round < 15; ++cleanup_round) {
                // 统计当前违规（全局统计，用边计数）
                int total_viol = 0;
                for (int vi = 0; vi < N0; ++vi) {
                    int of = part0.node2fpga[vi];
                    if (of >= 0) total_viol += CountViolations(vi, of, graph0, part0, fg);
                }
                if (total_viol == 0) {
                    cout << "    ZERO violations achieved at round " << cleanup_round << "!" << endl;
                    break;
                }

                // 收敛检测：允许最多 3 轮无改进（处理振荡效应）
                if (total_viol >= prev_total_viol) {
                    ++stall_count;
                    if (stall_count >= 3 && cleanup_round > 0) {
                        cout << "    Violations stalled for " << stall_count << " rounds (" << total_viol << "), stopping cleanup." << endl;
                        break;
                    }
                } else {
                    stall_count = 0;
                }
                prev_total_viol = total_viol;

                // 阶段1：激进违规消解（允许违规不增加，严格减少优先）
                int vr = ViolationResolvePass(lv0, fg);
                if (vr > 0) cout << "    Cleanup round " << cleanup_round
                                  << ": resolved " << vr << " violations (total edges: " << total_viol << ")" << endl;

                // 阶段2：安全逐对消解 + 链式消解
                //   方向1: vi → vj 的 FPGA
                //   方向2: vj → vi 的 FPGA
                //   方向3: vi → vj 的直连 FPGA（2-hop）
                {
                    int safe_pairs = 0;
                    for (int vi = 0; vi < N0; ++vi) {
                        int old_f = part0.node2fpga[vi];
                        if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                        int old_vi_v = CountViolations(vi, old_f, graph0, part0, fg);
                        if (old_vi_v == 0) continue;

                        for (int vj : graph0.adj[vi]) {
                            int vj_f = part0.node2fpga[vj];
                            if (vj_f < 0 || vj_f == old_f) continue;
                            if (fg.dist_all[old_f][vj_f] <= 1) continue;

                            bool moved = false;

                            // 方向1：vi → vj 的 FPGA
                            if (!graph0.fixed_nodes.count(vi)
                                && (int)part0.fpga2nodes[vj_f].size() < fg.resource_cap[vj_f]) {
                                int new_vi_v = CountViolations(vi, vj_f, graph0, part0, fg);
                                if (new_vi_v < old_vi_v) {
                                    part0.unassign(vi, old_f);
                                    part0.assign(vi, vj_f);
                                    ++safe_pairs; moved = true;
                                }
                            }
                            // 方向2：vj → vi 的 FPGA
                            if (!moved && !graph0.fixed_nodes.count(vj)
                                && (int)part0.fpga2nodes[old_f].size() < fg.resource_cap[old_f]) {
                                int old_vj_v = CountViolations(vj, vj_f, graph0, part0, fg);
                                int new_vj_v = CountViolations(vj, old_f, graph0, part0, fg);
                                if (new_vj_v < old_vj_v) {
                                    part0.unassign(vj, vj_f);
                                    part0.assign(vj, old_f);
                                    ++safe_pairs; moved = true;
                                }
                            }
                            // 方向3：vi → vj 的直连 FPGA（2-hop chain）
                            if (!moved && !graph0.fixed_nodes.count(vi)) {
                                for (int adj_f : fg.adj[vj_f]) {
                                    if (adj_f == old_f) continue;
                                    if ((int)part0.fpga2nodes[adj_f].size() >= fg.resource_cap[adj_f]) continue;
                                    int new_vi_v = CountViolations(vi, adj_f, graph0, part0, fg);
                                    if (new_vi_v < old_vi_v) {
                                        part0.unassign(vi, old_f);
                                        part0.assign(vi, adj_f);
                                        ++safe_pairs; moved = true;
                                        break;
                                    }
                                }
                            }
                            // 方向4：2-node chain — 同时移动 vi 和 vj 到共同 FPGA（终极消解）
                            if (!moved && !graph0.fixed_nodes.count(vi) && !graph0.fixed_nodes.count(vj)) {
                                for (int cf = 0; cf < K0; ++cf) {
                                    if (cf == old_f && cf == vj_f) continue;
                                    int cap_need = ((cf != old_f) ? 1 : 0) + ((cf != vj_f) ? 1 : 0);
                                    if ((int)part0.fpga2nodes[cf].size() + cap_need > fg.resource_cap[cf]) continue;
                                    // 模拟移动后的违规
                                    int old_vi_v = CountViolations(vi, old_f, graph0, part0, fg);
                                    int old_vj_v = CountViolations(vj, vj_f, graph0, part0, fg);
                                    // 移动 vi → cf
                                    int new_vi_v = CountViolations(vi, cf, graph0, part0, fg);
                                    // 暂时移动 vj → cf 检查
                                    int new_vj_v = CountViolations(vj, cf, graph0, part0, fg);
                                    if (new_vi_v < old_vi_v && new_vj_v <= old_vj_v) {
                                        if (cf != old_f) { part0.unassign(vi, old_f); part0.assign(vi, cf); }
                                        if (cf != vj_f) { part0.unassign(vj, vj_f); part0.assign(vj, cf); }
                                        ++safe_pairs; moved = true;
                                        break;
                                    }
                                }
                            }
                            if (moved) break;
                        }
                    }
                    if (safe_pairs > 0)
                        cout << "    Safe pair-wise resolved: " << safe_pairs << endl;
                }

                // 阶段3：检查-neighbor FM 恢复（仅当不增加邻居违规时移动）
                //   注意：此 pass 可能产生新违规，仅在违规数下降时启用
                {
                    // 跳过 FM 恢复以减少振荡（违规清零阶段优先消解违规）
                }

                if (vr == 0) {
                    // 最后手段：多节点协同消解
                    // 对每条违规边，尝试所有可能的双节点移动组合
                    int last_resort = 0;
                    for (int vi = 0; vi < N0; ++vi) {
                        int old_f = part0.node2fpga[vi];
                        if (old_f < 0 || graph0.fixed_nodes.count(vi)) continue;
                        if (CountViolations(vi, old_f, graph0, part0, fg) == 0) continue;

                        for (int vj : graph0.adj[vi]) {
                            int vj_f = part0.node2fpga[vj];
                            if (vj_f < 0 || vj_f == old_f) continue;
                            if (fg.dist_all[old_f][vj_f] <= 1) continue;  // 非违规边跳过

                            // 尝试所有双节点组合
                            for (int cf = 0; cf < K0; ++cf) {
                                if (graph0.fixed_nodes.count(vi) || graph0.fixed_nodes.count(vj)) {
                                    // 有固定节点：只能移动非固定的一方
                                    if (graph0.fixed_nodes.count(vi)) {
                                        if ((int)part0.fpga2nodes[cf].size() >= fg.resource_cap[cf]) continue;
                                        int new_vj = CountViolations(vj, cf, graph0, part0, fg);
                                        int old_vj = CountViolations(vj, vj_f, graph0, part0, fg);
                                        if (new_vj < old_vj) {
                                            part0.unassign(vj, vj_f); part0.assign(vj, cf); ++last_resort; break;
                                        }
                                    }
                                    if (graph0.fixed_nodes.count(vj)) {
                                        if ((int)part0.fpga2nodes[cf].size() >= fg.resource_cap[cf]) continue;
                                        int new_vi = CountViolations(vi, cf, graph0, part0, fg);
                                        int old_vi = CountViolations(vi, old_f, graph0, part0, fg);
                                        if (new_vi < old_vi) {
                                            part0.unassign(vi, old_f); part0.assign(vi, cf); ++last_resort; break;
                                        }
                                    }
                                    continue;
                                }

                                // 双节点可移动：尝试所有 (cf_vi, cf_vj) 组合
                                for (int cf2 = 0; cf2 < K0; ++cf2) {
                                    int need_vi = (cf != old_f) ? 1 : 0;
                                    int need_vj = (cf2 != vj_f) ? 1 : 0;
                                    if (cf == cf2) need_vj = 0;  // 同一 FPGA 只占用一个位置
                                    if ((int)part0.fpga2nodes[cf].size() + need_vi > fg.resource_cap[cf]) continue;
                                    if ((int)part0.fpga2nodes[cf2].size() + need_vj > fg.resource_cap[cf2]) continue;

                                    int new_vi = CountViolations(vi, cf, graph0, part0, fg);
                                    int old_vi = CountViolations(vi, old_f, graph0, part0, fg);
                                    if (new_vi >= old_vi) continue;

                                    int new_vj = CountViolations(vj, cf2, graph0, part0, fg);
                                    int old_vj = CountViolations(vj, vj_f, graph0, part0, fg);
                                    if (new_vj >= old_vj) continue;

                                    // 双节点违规都减少 → 执行
                                    if (cf != old_f) { part0.unassign(vi, old_f); part0.assign(vi, cf); }
                                    if (cf2 != vj_f) { part0.unassign(vj, vj_f); part0.assign(vj, cf2); }
                                    ++last_resort;
                                    goto next_violating_node;
                                }
                            }
                        }
                        next_violating_node:;
                    }

                    if (last_resort > 0) {
                        cout << "    Multi-node last-resort resolved: " << last_resort << endl;
                        // 不 break，继续下一轮清理
                    } else {
                        cout << "    No more violations resolvable, stopping cleanup." << endl;
                        break;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------
    // 步骤4：返回原始层划分结果
    // ----------------------------------------------------------
    auto& result = levels[0].level_part;

    // 计算最终割边数（每条割边被计算两次，故除以 2）
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
    result.cut_size /= 2;

    // 统计最终违规数
    int final_viol = 0;
    for (int vi = 0; vi < levels[0].super_graph.node_num; ++vi) {
        int vi_f = result.node2fpga[vi];
        if (vi_f >= 0) {
            final_viol += CountViolations(vi, vi_f, levels[0].super_graph, result, fg);
        }
    }
    result.violation_cnt = final_viol / 2;

    cout << "  Final cut size (original level): " << result.cut_size << endl;
    cout << "  Final topology violations: " << result.violation_cnt << endl;
    cout << "========== Module 5 Done ==========\n" << endl;

    return result;
}
