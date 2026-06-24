/**
 * @file    utils.cpp
 * @brief   工具函数实现
 */

#include "utils.h"
#include <queue>
#include <iostream>

// ================================================================
//  1. Dijkstra 全点对最短距离计算
// ================================================================

vector<vector<int>> CalcAllPairShortestPath(const vector<vector<int>>& adj, int node_cnt) {
    // 初始化距离矩阵，默认值为无穷大
    const int INF = INT_MAX / 2;
    vector<vector<int>> dist(node_cnt, vector<int>(node_cnt, INF));

    // 对每个源点运行 Dijkstra 算法
    for (int src = 0; src < node_cnt; ++src) {
        dist[src][src] = 0;

        // 小顶堆优先队列：(距离, 节点id)
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
        pq.push({0, src});

        while (!pq.empty()) {
            auto [d, u] = pq.top();
            pq.pop();

            // 懒惰删除：若已有更短路径则跳过
            if (d > dist[src][u]) continue;

            for (int v : adj[u]) {
                int nd = d + 1;  // 无权图，边权 = 1
                if (nd < dist[src][v]) {
                    dist[src][v] = nd;
                    pq.push({nd, v});
                }
            }
        }
    }

    return dist;
}

// ================================================================
//  2. 构建距离 ≤ x 的节点集合
// ================================================================

unordered_set<int> BuildS(const vector<vector<int>>& dist_all, int center, int max_x) {
    unordered_set<int> result;
    int n = (int)dist_all.size();
    for (int v = 0; v < n; ++v) {
        if (dist_all[center][v] <= max_x) {
            result.insert(v);
        }
    }
    return result;
}

// ================================================================
//  3. FPGA 图距离预计算
// ================================================================

void BuildFPGADist(FPGAGraph& fg) {
    int K = fg.fpga_num;

    // 步骤1：计算全点对最短距离
    fg.dist_all = CalcAllPairShortestPath(fg.adj, K);

    // 步骤2：计算每个 FPGA 的 max_dist
    fg.max_dist.resize(K, 0);
    for (int i = 0; i < K; ++i) {
        int max_d = 0;
        for (int j = 0; j < K; ++j) {
            if (fg.dist_all[i][j] > max_d && fg.dist_all[i][j] < INT_MAX / 2) {
                max_d = fg.dist_all[i][j];
            }
        }
        fg.max_dist[i] = max_d;
    }

    // 步骤3：预缓存所有 Ŝ(v̂i, x)
    fg.build_S_cache();

    cout << "[BuildFPGADist] FPGA all-pair distance computed: "
         << K << " FPGAs" << endl;
}

// ================================================================
//  4. 电路图距离预计算
// ================================================================

void BuildCircuitDist(CircuitGraph& g) {
    // 内存优化：不再预计算 N×N 全点对距离矩阵（N=300K 时需约 360GB 内存）
    // 电路距离在需要时按以下规则获取：
    //   - 相邻节点（遍历 adj[]）：距离 = 1
    //   - 候选传播阶段：使用 BFSFromSource 按需计算
    // 仅初始化 dist_all 为空（下游代码通过 dist_all.empty() 判断是否已计算）
    g.dist_all.clear();
    cout << "[BuildCircuitDist] Circuit distance: using on-demand BFS (O(1) memory)"
         << " for " << g.node_num << " nodes" << endl;
}
