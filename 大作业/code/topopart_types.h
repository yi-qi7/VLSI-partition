#ifndef TOPOPART_TYPES_H
#define TOPOPART_TYPES_H

/**
 * @file    topopart_types.h
 * @brief   TopoPart 全局数据结构定义（严格遵循 ICCAD2021 论文术语）
 *
 * 论文：TopoPart: a Multi-level Topology-Driven Partitioning Framework
 *       for Multi-FPGA Systems (ICCAD 2021)
 *
 * 核心约束（三大硬约束）：
 *   约束1 固定节点约束：固定电路节点必须分配至指定FPGA；
 *   约束2 资源约束：每个FPGA承载节点总数 ≤ 资源容量 r_i；
 *   约束3 拓扑约束：任意二引脚网两端节点，同FPGA OR 两块FPGA存在直连边，零拓扑违规。
 *
 * 核心数学定理（定理 III.1）：
 *   电路节点 vi, vj 最短电路距离 x；
 *   对应 FPGA v̂i, v̂j 硬件最短距离 y；
 *   vi 分配 v̂i、vj 分配 v̂j 无拓扑违规 ⇔ x ≥ y。
 */

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <intrin.h>  // __popcnt64, _BitScanForward64

using namespace std;

// ================================================================
//  1. FPGA 硬件拓扑图 Ĝ(Ê, V̂)
// ================================================================
struct FPGAGraph {
    int fpga_num;                       // |V̂| FPGA 总数量
    vector<vector<int>> adj;            // FPGA 邻接表，存储直连 FPGA 编号 [0..fpga_num-1]
    vector<int> resource_cap;           // r_i 每个 FPGA 资源上限（可承载最大节点数）
    vector<vector<int>> dist_all;       // dist(v̂i, v̂j)：任意两块 FPGA 最短距离矩阵
    vector<int> max_dist;               // maxDist(v̂i)：每个 FPGA 到所有其他 FPGA 的最大最短距离

    /**
     * Ŝ(v̂i, x)：距离 v̂i 不超过 x 的 FPGA 集合（预计算缓存）
     * S_cache[center_fpga][x] = 距离 center_fpga ≤ x 的所有 FPGA 集合
     * 实现为 vector<unordered_set<int>>，索引：S_cache[center_fpga * max_possible_x + x]
     * 为简化，使用二维 vector：S_cache[center][x]
     */
    vector<vector<unordered_set<int>>> S_cache;

    /** 初始化 S_cache：预计算所有 Ŝ(v̂i, x)，x ∈ [0, max_dist[v̂i]] */
    void build_S_cache();

    /** 获取 Ŝ(v̂i, x)：距离 center_fpga ≤ x 的 FPGA 集合引用 */
    const unordered_set<int>& get_S(int center_fpga, int x) const;
};

// ================================================================
//  2. 电路图 G(E, V) —— 全部 n 引脚网转为二引脚网 star 模型
// ================================================================
struct CircuitGraph {
    int node_num;                       // 电路总节点数
    vector<vector<int>> adj;            // 电路节点邻接表 Nbrs(vi)：相邻节点编号 [0..node_num-1]
    vector<vector<int>> dist_all;       // dist(vi, vj)：电路节点间最短距离矩阵

    unordered_set<int> fixed_nodes;     // V_γ：固定节点集合
    unordered_map<int, int> node2fpga;  // 固定节点映射：固定节点 id → 绑定 FPGA id [0..fpga_num-1]
    unordered_set<int> move_nodes;      // V_λ：可移动节点集合

    /** 初始化可移动节点集合（非固定节点） */
    void init_move_nodes();

    void build_dist();
};

// ================================================================
//  3. 候选 FPGA 集合 Cddt(vi)
// ================================================================
struct CandidateSet {
    // ★ Bitset 优化: K≤64 时用 uint64_t 替代 unordered_set
    //   交集/插入/大小统计全部按位运算，速度提升 10-100x
    uint64_t bits;                      // 位图: bit f = 1 表示 FPGA f 在候选集中

    /**
     * T_vec：计数向量，长度 = FPGA 总数
     * T_vec[f] 记录 FPGA f 被多少固定节点「推荐」过
     * 用于快速更新 / 恢复候选集：T_vec[f] > 0 表示 f 仍在候选集中
     */
    vector<int> T_vec;

    CandidateSet() : bits(0) {}
    explicit CandidateSet(int fpga_num);

    /** 将候选集重置为全部 FPGA */
    void reset_all(int fpga_num);

    /** 候选集与另一个集合求交 */
    void intersect_with(const unordered_set<int>& other, int fpga_num);

    /** 候选集与另一个集合求交（返回 true 表示发生变化） */
    bool intersect_with_check(const unordered_set<int>& other, int fpga_num);

    /** 候选集与位图求交（bitset 版本，更快） */
    void intersect_with_bits(uint64_t other_bits);

    /** 从候选集中移除单个 FPGA */
    void remove_fpga(int f, int fpga_num);

    /** 检查候选集是否为空 */
    bool empty() const { return bits == 0; }

    /** 检查候选集是否仅剩一个 FPGA */
    bool is_singleton() const { return bits != 0 && (bits & (bits - 1)) == 0; }

    /** 获取唯一候选 FPGA */
    int get_only_fpga() const {
        if (!is_singleton()) {
            throw runtime_error("CandidateSet::get_only_fpga(): not a singleton!");
        }
        unsigned long idx;
        _BitScanForward64(&idx, bits);
        return (int)idx;
    }

    /** 候选集大小 */
    int size() const { return (int)__popcnt64(bits); }

    /** 检查 FPGA f 是否在候选集中 */
    bool contains(int f) const { return (bits >> f) & 1; }

    /** 迭代器：遍历所有候选 FPGA */
    template<typename Func>
    void for_each(Func&& fn) const {
        uint64_t b = bits;
        while (b) {
            int f; _BitScanForward64((unsigned long*)&f, b);
            fn(f);
            b &= b - 1;
        }
    }
};

// ================================================================
//  4. 划分结果存储 P = {p_i}
// ================================================================
struct PartitionResult {
    vector<int> node2fpga;              // part(vi)：电路节点 → 分配 FPGA，-1 代表未分配
    vector<vector<int>> fpga2nodes;     // p_i：每个 FPGA 承载的所有电路节点列表
    long long cut_size;                 // 当前跨 FPGA 割边数量
    int violation_cnt;                  // 拓扑违规数量（TopoPart 全程保证 = 0）

    PartitionResult() : cut_size(0), violation_cnt(0) {}

    /** 初始化指定节点数的划分结果 */
    void init(int node_cnt, int fpga_cnt);

    /** 分配节点 vi 到 FPGA f */
    void assign(int vi, int f);

    /** 撤销节点 vi 的分配 */
    void unassign(int vi, int f);
};

// ================================================================
//  5. 多层粗化层级存储
// ================================================================
struct CoarsenLevel {
    CircuitGraph super_graph;               // 当前层超节点电路图（超图）
    vector<CandidateSet> super_cddt;        // 每个超节点的候选 FPGA 集合
    unordered_map<int, vector<int>> super2origin; // 超节点 id → 原始子节点 id 列表
    vector<int> super_weight;              // 每个超节点代表的原始节点数（用于加权负载均衡）
    PartitionResult level_part;             // 该层划分结果

    CoarsenLevel() = default;
};

// ================================================================
//  6. 算法全局配置参数
// ================================================================
struct TopoPartConfig {
    double beta_ratio = 2.0 / 3.0;      // β 初始合并阈值比例（β = beta_ratio * fpga_num）
    int    max_coarsen_retry = 5;        // 粗化失败最大重试次数（每次增大 β）
    double beta_grow_factor = 1.1;       // β 增大因子
    int    max_traceback_depth = 100;    // 初始划分最大回溯深度
    bool   verbose = true;               // 是否输出详细日志

    /** 根据 FPGA 数量计算初始 β 阈值 */
    int get_beta(int fpga_num) const {
        return max(1, (int)(beta_ratio * fpga_num));
    }
};

#endif // TOPOPART_TYPES_H
