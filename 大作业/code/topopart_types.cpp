/**
 * @file    topopart_types.cpp
 * @brief   TopoPart 数据结构成员函数实现
 */

#include "topopart_types.h"
#include "topopart_utils.h"
#include <stdexcept>
#include <algorithm>

// ================================================================
//  FPGAGraph 成员函数
// ================================================================

void FPGAGraph::build_S_cache() {
    // S_cache[center_fpga][x]：距离 center_fpga ≤ x 的所有 FPGA 集合
    // 预分配：S_cache.resize(fpga_num)
    S_cache.clear();
    S_cache.resize(fpga_num);

    for (int center = 0; center < fpga_num; ++center) {
        int max_d = max_dist[center];
        S_cache[center].resize(max_d + 1);

        // 按距离桶分组：bucket[d] = 距离 center 恰好为 d 的 FPGA
        vector<vector<int>> bucket(max_d + 1);
        for (int other = 0; other < fpga_num; ++other) {
            int d = dist_all[center][other];
            if (d <= max_d) {
                bucket[d].push_back(other);
            }
        }

        // 前缀累积：Ŝ(v̂i, x) = ∪_{d=0}^{x} bucket[d]
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
    // 边界安全：x 可能超过 max_dist[center_fpga]
    if (x > max_dist[center_fpga]) {
        x = max_dist[center_fpga];
    }
    if (x < 0) x = 0;
    return S_cache[center_fpga][x];
}

// ================================================================
//  CircuitGraph 成员函数
// ================================================================

void CircuitGraph::init_move_nodes() {
    move_nodes.clear();
    for (int i = 0; i < node_num; ++i) {
        if (fixed_nodes.find(i) == fixed_nodes.end()) {
            move_nodes.insert(i);
        }
    }
}

void CircuitGraph::build_dist() {
    dist_all = CalcAllPairShortestPath(adj, node_num);
}

// ================================================================
//  CandidateSet 成员函数
// ================================================================

CandidateSet::CandidateSet(int fpga_num) : bits(0) {
    T_vec.assign(fpga_num, 0);
}

void CandidateSet::reset_all(int fpga_num) {
    bits = (fpga_num >= 64) ? ~0ULL : ((1ULL << fpga_num) - 1);
    T_vec.assign(fpga_num, 0);
    for (int f = 0; f < fpga_num; ++f) T_vec[f] = 1;
}

void CandidateSet::intersect_with(const unordered_set<int>& other, int fpga_num) {
    uint64_t other_bits = 0;
    for (int f : other) other_bits |= (1ULL << f);
    uint64_t old = bits;
    bits &= other_bits;
    uint64_t removed = old ^ bits;
    while (removed) {
        int f; _BitScanForward64((unsigned long*)&f, removed);
        if (f < (int)T_vec.size()) T_vec[f] = 0;
        removed &= removed - 1;
    }
}

bool CandidateSet::intersect_with_check(const unordered_set<int>& other, int fpga_num) {
    uint64_t old = bits;
    intersect_with(other, fpga_num);
    return bits != old;
}

void CandidateSet::intersect_with_bits(uint64_t other_bits) {
    bits &= other_bits;
}

void CandidateSet::remove_fpga(int f, int fpga_num) {
    bits &= ~(1ULL << f);
    if (f < (int)T_vec.size()) T_vec[f] = 0;
}

// get_only_fpga() is defined inline in types.h

// ================================================================
//  PartitionResult 成员函数
// ================================================================

void PartitionResult::init(int node_cnt, int fpga_cnt) {
    node2fpga.assign(node_cnt, -1);
    fpga2nodes.assign(fpga_cnt, vector<int>());
    cut_size = 0;
    violation_cnt = 0;
}

void PartitionResult::assign(int vi, int f) {
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
    if (vi < 0 || vi >= (int)node2fpga.size()) return;
    node2fpga[vi] = -1;
    auto& vec = fpga2nodes[f];
    vec.erase(remove(vec.begin(), vec.end(), vi), vec.end());
}
