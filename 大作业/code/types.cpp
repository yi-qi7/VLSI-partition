/**
 * @file    types.cpp
 * @brief   核心数据结构成员函数实现
 */

#include "types.h"
#include "utils.h"
#include <stdexcept>
#include <algorithm>

// ================================================================
//  FPGAGraph 成员函数实现
// ================================================================

void FPGAGraph::build_S_cache() {
    // 预分配 S_cache 空间
    S_cache.clear();
    S_cache.resize(fpga_num);

    for (int center = 0; center < fpga_num; ++center) {
        int max_d = max_dist[center];
        S_cache[center].resize(max_d + 1);

        // 按距离分桶：bucket[d] 为距 center 恰好为 d 的所有 FPGA
        vector<vector<int>> bucket(max_d + 1);
        for (int other = 0; other < fpga_num; ++other) {
            int d = dist_all[center][other];
            if (d <= max_d) {
                bucket[d].push_back(other);
            }
        }

        // 前缀累积构建 Ŝ 缓存：Ŝ(v̂i, x) = ∪_{d=0}^{x} bucket[d]
        unordered_set<int> accum;
        for (int x = 0; x <= max_d; ++x) {
            for (int f : bucket[x]) {
                accum.insert(f);
            }
            S_cache[center][x] = accum;
        }
    }
}

const unordered_set<int>& FPGAGraph::get_S(int center_fpga, int x) const {
    // 边界保护：将 x 裁剪到 [0, max_dist] 范围内
    if (x > max_dist[center_fpga]) {
        x = max_dist[center_fpga];
    }
    if (x < 0) x = 0;
    return S_cache[center_fpga][x];
}

// ================================================================
//  CircuitGraph 成员函数实现
// ================================================================

void CircuitGraph::init_move_nodes() {
    // 可移动节点 = 全部节点 - 固定节点
    move_nodes.clear();
    for (int i = 0; i < node_num; ++i) {
        if (fixed_nodes.find(i) == fixed_nodes.end()) {
            move_nodes.insert(i);
        }
    }
}

void CircuitGraph::build_dist() {
    // 委托给工具函数计算全点对最短路径
    dist_all = CalcAllPairShortestPath(adj, node_num);
}

// ================================================================
//  CandidateSet 成员函数实现
// ================================================================

CandidateSet::CandidateSet(int fpga_num) : bits(0) {
    T_vec.assign(fpga_num, 0);
}

void CandidateSet::reset_all(int fpga_num) {
    // 将所有 FPGA 置入候选集（bits 全 1）
    bits = (fpga_num >= 64) ? ~0ULL : ((1ULL << fpga_num) - 1);
    T_vec.assign(fpga_num, 0);
    for (int f = 0; f < fpga_num; ++f) T_vec[f] = 1;
}

void CandidateSet::intersect_with(const unordered_set<int>& other) {
    // 位图交集运算，同步更新 T_vec 计数向量
    uint64_t other_bits = 0;
    for (int f : other) other_bits |= (1ULL << f);
    uint64_t old = bits;
    bits &= other_bits;
    uint64_t removed = old ^ bits;
    while (removed) {
        unsigned long f_ul; _BitScanForward64(&f_ul, removed);
        int f = (int)f_ul;
        if (f < (int)T_vec.size()) T_vec[f] = 0;
        removed &= removed - 1;
    }
}

bool CandidateSet::intersect_with_check(const unordered_set<int>& other) {
    // 执行交集并返回是否发生了变化
    uint64_t old = bits;
    intersect_with(other);
    return bits != old;
}

void CandidateSet::intersect_with_bits(uint64_t other_bits) {
    // 直接使用位图进行交集（更快，跳过 T_vec 同步）
    bits &= other_bits;
}

void CandidateSet::remove_fpga(int f) {
    // 从候选集位图中永久移除指定 FPGA（T_vec 归零防止回溯时被意外恢复）
    bits &= ~(1ULL << f);
    if (f < (int)T_vec.size()) T_vec[f] = -1000000;  // 永久排除
}

// get_only_fpga() 在 types.h 中以内联方式定义

// ================================================================
//  PartitionResult 成员函数实现
// ================================================================

void PartitionResult::init(int node_cnt, int fpga_cnt) {
    // 初始化指定规模的划分结果：所有节点未分配（-1）
    node2fpga.assign(node_cnt, -1);
    fpga2nodes.assign(fpga_cnt, vector<int>());
    cut_size = 0;
    violation_cnt = 0;
}

void PartitionResult::assign(int vi, int f) {
    // 将节点 vi 分配到 FPGA f，含越界检查
    if (vi < 0 || vi >= (int)node2fpga.size()) {
        throw out_of_range("PartitionResult::assign: node id out of range");
    }
    if (f < 0 || f >= (int)fpga2nodes.size()) {
        throw out_of_range("PartitionResult::assign: FPGA id out of range");
    }
    node2fpga[vi] = f;
    fpga2nodes[f].push_back(vi);
}

void PartitionResult::unassign(int vi, int f) {
    // 撤销节点 vi 的分配，从 FPGA f 的节点列表中移除
    if (vi < 0 || vi >= (int)node2fpga.size()) return;
    node2fpga[vi] = -1;
    auto& vec = fpga2nodes[f];
    vec.erase(remove(vec.begin(), vec.end(), vi), vec.end());
}
