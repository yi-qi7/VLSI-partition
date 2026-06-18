/**
 * @file    topopart_coarsen.cpp
 * @brief   TopoPart 模块3：候选集感知多层粗化 Coarsening 实现
 *
 * 合并规则（双条件同时满足才允许合并）：
 *   条件1：两节点电路互连权重高（有直接边相连）
 *   条件2：两节点候选集交集大小 ≥ β
 *
 * 超节点候选集：super_cddt = cddt_a ∩ cddt_b
 *
 * 迭代粗化流程：
 *   每轮扫描所有节点对，贪心匹配满足双条件的节点对，
 *   合并为超节点，直到图规模无法再压缩。
 */

#include "topopart_coarsen.h"
#include "topopart_utils.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <set>

// ================================================================
//  辅助：判断两个节点是否可以合并
// ================================================================

/**
 * @brief 检查节点 u 和 v 是否满足合并双条件
 * @param u, v      节点 id
 * @param adj       当前层邻接表
 * @param cddt      当前层候选集数组
 * @param beta      合并阈值
 * @return true 表示允许合并
 */
static bool CanMerge(int u, int v,
                     const vector<vector<int>>& adj,
                     const vector<CandidateSet>& cddt,
                     int beta) {
    // 条件1：两节点必须有边相连（互连权重高）
    // 简化：检查是否存在直接邻接关系
    bool has_edge = false;
    for (int nb : adj[u]) {
        if (nb == v) { has_edge = true; break; }
    }
    if (!has_edge) return false;

    // 条件2：候选集交集大小 ≥ β (bitset 按位运算)
    uint64_t inter = cddt[u].bits & cddt[v].bits;
    return __popcnt64(inter) >= beta;
}

/**
 * @brief 计算两候选集交集（返回新候选集）
 */
static CandidateSet ComputeSuperCddt(const CandidateSet& a, const CandidateSet& b, int K) {
    CandidateSet result(K);
    result.bits = a.bits & b.bits;
    // 复制 T_vec
    result.T_vec = a.T_vec;
    for (int f = 0; f < K; ++f) {
        if (!result.contains(f)) result.T_vec[f] = 0;
    }
    return result;
}

// ================================================================
//  单层粗化：从当前图构建超图
// ================================================================

/**
 * @brief 执行单层粗化
 * @param curr_graph     当前层电路图
 * @param curr_cddt      当前层候选集
 * @param beta           合并阈值
 * @param next_level     输出：下一粗化层级
 * @return true 表示成功进行了合并（图规模缩小），false 表示无法再合并
 */
static bool DoOneCoarsenLevel(const CircuitGraph& curr_graph,
                               const vector<CandidateSet>& curr_cddt,
                               const vector<int>& prev_weight,
                               int beta,
                               CoarsenLevel& next_level) {
    int n = curr_graph.node_num;
    int K = (int)curr_cddt[0].T_vec.size();

    // ----------------------------------------------------------
    // 步骤1：贪心匹配
    //   使用 visited 标记已被合并的节点
    //   使用 match 记录匹配对 match[u] = v, match[v] = u
    // ----------------------------------------------------------
    vector<bool> visited(n, false);
    vector<int> match(n, -1);

    // 按节点度降序排序，优先匹配高度节点（更多连接 → 合并收益更大）
    vector<pair<int, int>> node_degree;  // (度, id)
    for (int i = 0; i < n; ++i) {
        node_degree.push_back({(int)curr_graph.adj[i].size(), i});
    }
    sort(node_degree.begin(), node_degree.end(),
         [](const auto& a, const auto& b) { return a.first > b.first; });

    for (auto& [deg, u] : node_degree) {
        if (visited[u]) continue;
        if (curr_cddt[u].is_singleton()) {
            // 单候选节点暂不合并（避免过早固定）
            // 标记为已访问但单独成为一个超节点
            visited[u] = true;
            continue;
        }

        // 寻找最佳合并候选：候选集交集最大的邻居
        int best_v = -1;
        int best_inter = beta - 1;  // 至少需要 beta

        for (int v : curr_graph.adj[u]) {
            if (visited[v]) continue;
            if (curr_cddt[v].is_singleton()) continue;

            if (!CanMerge(u, v, curr_graph.adj, curr_cddt, beta)) continue;

            // 计算交集大小 (bitset 按位运算)
            int inter = (int)__popcnt64(curr_cddt[u].bits & curr_cddt[v].bits);
            if (inter > best_inter) {
                best_inter = inter;
                best_v = v;
            }
        }

        if (best_v >= 0) {
            // 匹配成功
            match[u] = best_v;
            match[best_v] = u;
            visited[u] = true;
            visited[best_v] = true;
        } else {
            // 无合适匹配，单独成超节点
            visited[u] = true;
        }
    }

    // ----------------------------------------------------------
    // 步骤2：构建超节点映射
    //   每个超节点由 1~2 个原始节点合并而成
    // ----------------------------------------------------------
    int super_cnt = 0;
    vector<int> node2super(n, -1);       // 原始节点 → 超节点 id
    unordered_map<int, vector<int>> super2origin;  // 超节点 id → 子节点列表

    for (int i = 0; i < n; ++i) {
        if (node2super[i] >= 0) continue;

        if (match[i] >= 0 && node2super[match[i]] < 0) {
            // 合并对
            int sid = super_cnt++;
            node2super[i] = sid;
            node2super[match[i]] = sid;
            super2origin[sid] = {i, match[i]};
        } else if (match[i] < 0) {
            // 单独成超节点
            int sid = super_cnt++;
            node2super[i] = sid;
            super2origin[sid] = {i};
        }
    }

    // 检查是否有效压缩
    if (super_cnt >= n) {
        // 没有合并发生
        return false;
    }

    // ----------------------------------------------------------
    // 步骤3：构建超节点邻接表
    //   若超节点 A 包含原始节点 a1,a2，超节点 B 包含原始节点 b1,b2，
    //   且 ai 与 bj 在原图有边 → A 与 B 有边
    // ----------------------------------------------------------
    next_level.super_graph.node_num = super_cnt;
    next_level.super_graph.adj.assign(super_cnt, vector<int>());
    next_level.super_graph.fixed_nodes.clear();
    next_level.super_graph.node2fpga.clear();

    set<pair<int, int>> super_edges;  // 去重用

    for (int u = 0; u < n; ++u) {
        int su = node2super[u];
        for (int v : curr_graph.adj[u]) {
            int sv = node2super[v];
            if (su != sv) {
                // 无向边，统一 (min, max) 去重
                auto edge = make_pair(min(su, sv), max(su, sv));
                if (super_edges.insert(edge).second) {
                    next_level.super_graph.adj[su].push_back(sv);
                    next_level.super_graph.adj[sv].push_back(su);
                }
            }
        }
    }

    // ----------------------------------------------------------
    // 步骤4：构建超节点候选集
    //   super_cddt = cddt_a ∩ cddt_b（若为合并对）
    //   单节点超节点 = 原始候选集
    // ----------------------------------------------------------
    next_level.super_cddt.resize(super_cnt, CandidateSet(K));

    for (int sid = 0; sid < super_cnt; ++sid) {
        const auto& children = super2origin[sid];

        if (children.size() == 1) {
            // 单节点超节点：直接继承候选集
            int orig = children[0];
            next_level.super_cddt[sid].bits = curr_cddt[orig].bits;
            next_level.super_cddt[sid].T_vec = curr_cddt[orig].T_vec;
        } else {
            // 合并超节点：候选集 = 交集
            int a = children[0], b = children[1];
            next_level.super_cddt[sid] = ComputeSuperCddt(curr_cddt[a], curr_cddt[b], K);

            // 超节点候选集为空 → 此合并无效
            if (next_level.super_cddt[sid].empty()) {
                cerr << "  WARNING super-node " << sid << " has empty candidate set!"
                     << " child " << a << " (|Cddt|=" << curr_cddt[a].size()
                     << ") + " << b << " (|Cddt|=" << curr_cddt[b].size() << ")"
                     << endl;
                return false;
            }
        }

        // 检查是否超节点继承固定属性
        if (next_level.super_cddt[sid].is_singleton()) {
            int only_fpga = next_level.super_cddt[sid].get_only_fpga();
            next_level.super_graph.fixed_nodes.insert(sid);
            next_level.super_graph.node2fpga[sid] = only_fpga;
        }
    }

    // ----------------------------------------------------------
    // 步骤5：初始化超节点可移动集合
    // ----------------------------------------------------------
    next_level.super_graph.init_move_nodes();

    // 保留原始节点映射
    next_level.super2origin = std::move(super2origin);

    // 计算超节点权重 = 子节点权重之和
    next_level.super_weight.resize(super_cnt, 0);
    for (int i = 0; i < n; ++i) {
        int sid = node2super[i];
        int w = (i < (int)prev_weight.size()) ? prev_weight[i] : 1;
        next_level.super_weight[sid] += w;
    }

    // 初始化该层划分结果
    next_level.level_part.init(super_cnt, K);

    // 应用固定约束
    for (int sid : next_level.super_graph.fixed_nodes) {
        int f = next_level.super_graph.node2fpga[sid];
        next_level.level_part.assign(sid, f);
    }

    cout << "  [Coarsen] " << n << " nodes -> " << super_cnt
         << " super-nodes (ratio " << (double)n / super_cnt << ":1)"
         << ", β=" << beta << endl;

    return true;
}

// ================================================================
//  主粗化函数
// ================================================================

vector<CoarsenLevel> Coarsening(const CircuitGraph& origin_g,
                                 const vector<CandidateSet>& node_cddt,
                                 const FPGAGraph& fg,
                                 int beta_init) {
    int K = fg.fpga_num;
    cout << "\n========== Module 3: Cddt-Aware Coarsening ==========" << endl;
    cout << "  Initial beta = " << beta_init << endl;

    vector<CoarsenLevel> levels;

    // Level 0：原始层
    {
        CoarsenLevel lv0;
        lv0.super_graph = origin_g;
        lv0.super_cddt = node_cddt;
        lv0.level_part.init(origin_g.node_num, K);

        // 构建原始节点映射（恒等映射）+ 权重（每个原始节点权重 = 1）
        lv0.super_weight.resize(origin_g.node_num, 1);
        for (int i = 0; i < origin_g.node_num; ++i) {
            lv0.super2origin[i] = {i};
        }
        // 应用固定约束
        for (int vi : origin_g.fixed_nodes) {
            int f = origin_g.node2fpga.at(vi);
            lv0.level_part.assign(vi, f);
        }

        levels.push_back(std::move(lv0));
    }

    // 迭代粗化循环
    int beta = beta_init;
    int min_super_size = max(5, K * 2);  // 最少超节点数（保证足够自由度）
    int stall_count = 0;                  // 压缩停滞计数器
    int prev_node_cnt = levels.back().super_graph.node_num;

    for (int round = 1; ; ++round) {
        const auto& prev = levels.back();
        CoarsenLevel next;

        bool success = DoOneCoarsenLevel(prev.super_graph, prev.super_cddt,
                                          prev.super_weight, beta, next);
        if (!success) {
            cout << "  [Coarsen] Round " << round << ": no merges, coarsening stops." << endl;
            break;
        }

        // 检查候选集是否有空
        bool has_empty = false;
        for (int i = 0; i < next.super_graph.node_num; ++i) {
            if (next.super_cddt[i].empty()) {
                cerr << "  WARNING super-node " << i << " has empty candidate set!" << endl;
                has_empty = true;
            }
        }
        if (has_empty) {
            cout << "  [Coarsen] Empty candidate set detected, rolling back this round." << endl;
            break;
        }

        levels.push_back(std::move(next));

        // 终止条件：图规模足够小
        if (levels.back().super_graph.node_num <= min_super_size) {
            cout << "  [Coarsen] Reached minimum super-node count "
                 << levels.back().super_graph.node_num << ", coarsening stops." << endl;
            break;
        }

        // 早期终止：压缩率持续过低（连续 3 轮压缩率 < 1.02）
        int curr_node_cnt = levels.back().super_graph.node_num;
        double ratio = (double)prev_node_cnt / curr_node_cnt;
        if (ratio < 1.02) {
            ++stall_count;
            if (stall_count >= 3) {
                cout << "  [Coarsen] Compression stalled (ratio " << ratio
                     << " < 1.02 for " << stall_count << " rounds), stopping early." << endl;
                break;
            }
        } else {
            stall_count = 0;
        }
        prev_node_cnt = curr_node_cnt;

        // 安全上限
        if (round >= 50) {
            cout << "  [Coarsen] Reached max coarsening rounds (50), stopping." << endl;
            break;
        }
    }

    cout << "  Coarsening levels: " << levels.size()
         << " (coarsest level nodes: " << levels.back().super_graph.node_num << ")" << endl;
    cout << "========== Module 3 Done ==========\n" << endl;

    return levels;
}

// ================================================================
//  带重试的粗化
// ================================================================

vector<CoarsenLevel> CoarseningWithRetry(const CircuitGraph& origin_g,
                                          const vector<CandidateSet>& node_cddt,
                                          const FPGAGraph& fg,
                                          const TopoPartConfig& cfg) {
    int beta = cfg.get_beta(fg.fpga_num);

    for (int retry = 0; retry < cfg.max_coarsen_retry; ++retry) {
        if (cfg.verbose) {
            cout << "\n[CoarseningWithRetry] Attempt #" << (retry + 1)
                 << ", beta = " << beta << endl;
        }

        auto levels = Coarsening(origin_g, node_cddt, fg, beta);

        // 检查最粗层是否有空候选集
        bool valid = true;
        for (auto& cddt : levels.back().super_cddt) {
            if (cddt.empty()) { valid = false; break; }
        }

        if (valid) {
            return levels;
        }

        // 增大 β 重试
        beta = max(beta + 1, (int)(beta * cfg.beta_grow_factor));
        cout << "  [CoarseningWithRetry] Empty candidate set in coarsening, increasing beta -> " << beta << endl;
    }

    throw runtime_error("CoarseningWithRetry: exceeded max retries, "
                        "cannot find valid coarsening plan!");
}
