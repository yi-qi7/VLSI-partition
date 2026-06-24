#ifndef UTILS_H
#define UTILS_H

/**
 * @file    utils.h
 * @brief   通用工具函数（最短路径、候选集构建）
 *
 * 包含：
 *   1. CalcAllPairShortestPath  —— Dijkstra 全点对最短路径
 *   2. BuildS                   —— 构建距离中心点 ≤ x 的节点集合
 *   3. BuildFPGADist            —— 构建 FPGA 图全点对距离 + max_dist + S_cache
 *   4. BuildCircuitDist         —— 电路图距离预计算（按需 BFS，O(1) 内存）
 */

#include "types.h"
#include <queue>
#include <algorithm>
#include <climits>

// ================================================================
//  1. Dijkstra 全点对最短距离
// ================================================================
/**
 * @brief 计算无向图全节点对最短距离（Dijkstra 算法）
 * @param adj       邻接表，adj[u] 存储节点 u 的所有邻居
 * @param node_cnt  节点总数
 * @return 全节点对最短距离矩阵 dist[u][v] = 节点 u 到 v 的最短距离
 *
 * 复杂度：O(node_cnt * (E + V log V))
 * 用于：FPGA 图 max_dist 计算、电路节点距离计算
 */
vector<vector<int>> CalcAllPairShortestPath(const vector<vector<int>>& adj, int node_cnt);

// ================================================================
//  2. 构建距离 ≤ x 的节点集合
// ================================================================
/**
 * @brief 构建 S(center, max_x)：距离中心点 ≤ max_x 的节点集合
 * @param dist_all  全节点对最短距离矩阵
 * @param center    中心节点 id
 * @param max_x     最大距离阈值
 * @return 距离 center ≤ max_x 的所有节点 id 集合
 *
 * 用于：FPGA 图 Ŝ(v̂i, x) 构建、电路图 S(vi, d) 构建
 */
unordered_set<int> BuildS(const vector<vector<int>>& dist_all, int center, int max_x);

// ================================================================
//  2b. 按需 BFS（替代全点对距离矩阵，避免 O(N²) 内存爆炸）
// ================================================================
/**
 * @brief 从源点运行 BFS，对每个访问到的节点调用回调
 * @param adj         邻接表
 * @param src         源节点 id
 * @param max_depth   BFS 最大深度（≤0 表示无限制）
 * @param callback    回调函数：callback(node_id, distance)
 *
 * 内存：O(V)（仅 BFS 队列 + visited 数组）
 * 用于：替代 N×N 全点对距离矩阵
 */
template<typename Func>
void BFSFromSource(const vector<vector<int>>& adj, int src, int max_depth, Func&& callback) {
    int n = (int)adj.size();
    vector<int> dist(n, -1);
    queue<int> q;

    dist[src] = 0;
    q.push(src);
    callback(src, 0);

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        int nd = dist[u] + 1;

        if (max_depth > 0 && nd > max_depth) continue;

        for (int v : adj[u]) {
            if (dist[v] >= 0) continue;  // 已访问过，跳过
            dist[v] = nd;
            q.push(v);
            callback(v, nd);
        }
    }
}

// ================================================================
//  3. FPGA 图距离预计算
// ================================================================
/**
 * @brief 为 FPGA 图计算全点对距离、max_dist、预缓存所有 Ŝ(v̂i, x)
 * @param fg  FPGA 图（传入引用，修改内部 dist_all / max_dist / S_cache）
 *
 * 步骤：
 *   1. 调用 CalcAllPairShortestPath 计算 dist_all
 *   2. 计算每个 FPGA 的 max_dist[v̂i] = max_j dist_all[v̂i][j]
 *   3. 对每个 FPGA v̂i 和每个 x ∈ [0, max_dist[v̂i]]，构建 Ŝ(v̂i, x) 并缓存
 */
void BuildFPGADist(FPGAGraph& fg);

// ================================================================
//  4. 电路图距离预计算（可选，按需调用）
// ================================================================
/**
 * @brief 为电路图计算全点对最短距离
 * @param g  电路图（传入引用，修改内部 dist_all）
 *
 * 当前实现采用按需 BFS 策略，不再预计算 N×N 矩阵，
 * 以节省大规模电路（N>100K）时的内存
 */
void BuildCircuitDist(CircuitGraph& g);

#endif // UTILS_H
