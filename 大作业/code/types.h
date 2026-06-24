#ifndef TYPES_H
#define TYPES_H

/**
 * @file    types.h
 * @brief   全局数据结构定义（遵循 ICCAD 2021 论文术语）
 *
 * 三大硬约束：
 *   约束1 固定节点约束：固定电路节点必须分配至指定 FPGA；
 *   约束2 资源约束：每个 FPGA 承载节点总数 ≤ 资源容量 r_i；
 *   约束3 拓扑约束：任意二引脚网两端节点，同 FPGA 或所在 FPGA 存在直连边。
 *
 * 定理 III.1（拓扑无违规判别条件）：
 *   设电路节点 vi, vj 最短电路距离为 x，
 *   对应 FPGA v̂i, v̂j 硬件最短距离为 y，则：
 *   vi→v̂i、vj→v̂j 无拓扑违规 ⟺ x ≥ y。
 */

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <stdexcept>

// 跨平台位运算内联函数：MSVC 用 <intrin.h>，GCC/Clang 用 __builtin_ 系列
#ifdef _MSC_VER
#include <intrin.h>  // 提供 __popcnt64, _BitScanForward64
#else
// GCC/Clang 使用内置函数替代
#define __popcnt64 __builtin_popcountll
// 用 __builtin_ctzll 模拟 _BitScanForward64：返回最低置位索引（调用者保证 bits != 0）
#define _BitScanForward64(p, bits64) \
    (*(p) = (unsigned long)__builtin_ctzll(bits64), 1)
#endif

using namespace std;

// ================================================================
//  1. FPGA 硬件拓扑图 Ĝ(Ê, V̂)
// ================================================================
struct FPGAGraph {
    int fpga_num;                       // |V̂|：FPGA 总数量
    vector<vector<int>> adj;            // FPGA 邻接表 [0..fpga_num-1]
    vector<int> resource_cap;           // r_i：每个 FPGA 资源容量上限
    vector<vector<int>> dist_all;       // dist(v̂i, v̂j)：FPGA 全点对最短距离矩阵
    vector<int> max_dist;               // maxDist(v̂i)：每个 FPGA 到其他 FPGA 的最大距离

    /**
     * Ŝ(v̂i, x)：距离 FPGA v̂i 不超过 x 的所有 FPGA 集合（预计算缓存）
     * 使用二维 vector 结构：S_cache[center][x]
     */
    vector<vector<unordered_set<int>>> S_cache;

    /** 预计算所有 Ŝ(v̂i, x)，其中 x ∈ [0, max_dist[v̂i]] */
    void build_S_cache();

    /** 获取 Ŝ(v̂i, x)：返回距离 center_fpga ≤ x 的 FPGA 集合引用 */
    const unordered_set<int>& get_S(int center_fpga, int x) const;
};

// ================================================================
//  2. 电路图 G(E, V) —— 多引脚网以 star 模型转为二引脚网
// ================================================================
struct CircuitGraph {
    int node_num;                       // 电路总节点数
    vector<vector<int>> adj;            // 邻接表 Nbrs(vi)：[0..node_num-1]
    vector<vector<int>> dist_all;       // dist(vi, vj)：电路节点最短距离矩阵

    unordered_set<int> fixed_nodes;     // V_γ：固定节点集合
    unordered_map<int, int> node2fpga;  // 固定节点 id → 绑定 FPGA id 的映射
    unordered_set<int> move_nodes;      // V_λ：可移动节点集合

    /** 初始化可移动节点集合（从全集中排除固定节点） */
    void init_move_nodes();

    /** 预计算电路图全点对最短距离 */
    void build_dist();
};

// ================================================================
//  3. 候选 FPGA 集合 Cddt(vi)
// ================================================================
struct CandidateSet {
    // 位图优化：当 K ≤ 64 时用 uint64_t 替代 unordered_set，
    // 交集 / 插入 / 大小统计全部使用按位运算，速度提升 10~100 倍
    uint64_t bits;                      // 位图：bit f = 1 表示 FPGA f 在候选集中

    /**
     * T_vec：计数向量，长度 = FPGA 总数
     * T_vec[f] 记录 FPGA f 被多少固定节点"推荐"过，
     * 用于快速更新 / 恢复候选集：T_vec[f] > 0 表示 f 仍在候选集中
     */
    vector<int> T_vec;

    CandidateSet() : bits(0) {}
    explicit CandidateSet(int fpga_num);

    /** 将候选集重置为全部 FPGA */
    void reset_all(int fpga_num);

    /** 候选集与另一个集合求交集 */
    void intersect_with(const unordered_set<int>& other);

    /** 候选集与另一个集合求交，返回 true 表示发生改变 */
    bool intersect_with_check(const unordered_set<int>& other);

    /** 候选集与位图求交（bitset 版本，速度更快） */
    void intersect_with_bits(uint64_t other_bits);

    /** 从候选集中移除单个 FPGA */
    void remove_fpga(int f);

    /** 检查候选集是否为空 */
    bool empty() const { return bits == 0; }

    /** 检查候选集是否仅剩唯一 FPGA（单元素） */
    bool is_singleton() const { return bits != 0 && (bits & (bits - 1)) == 0; }

    /** 获取唯一的候选 FPGA（要求候选集为单元素） */
    int get_only_fpga() const {
        if (!is_singleton()) {
            throw runtime_error("CandidateSet::get_only_fpga(): not a singleton!");
        }
        unsigned long idx;
        _BitScanForward64(&idx, bits);
        return (int)idx;
    }

    /** 返回候选集大小（即候选 FPGA 个数） */
    int size() const { return (int)__popcnt64(bits); }

    /** 检查 FPGA f 是否在候选集中 */
    bool contains(int f) const { return (bits >> f) & 1; }

    /** 迭代器：遍历候选集中的所有 FPGA */
    template<typename Func>
    void for_each(Func&& fn) const {
        uint64_t b = bits;
        while (b) {
            unsigned long f_ul; _BitScanForward64(&f_ul, b);
            int f = (int)f_ul;
            fn(f);
            b &= b - 1;
        }
    }
};

// ================================================================
//  4. 划分结果存储 P = {p_i}
// ================================================================
struct PartitionResult {
    vector<int> node2fpga;              // part(vi)：节点 → FPGA 映射，-1 表示未分配
    vector<vector<int>> fpga2nodes;     // 每个 FPGA 所承载的节点列表
    long long cut_size;                 // 跨 FPGA 割边总数
    int violation_cnt;                  // 拓扑违规数（算法保证最终为零）

    PartitionResult() : cut_size(0), violation_cnt(0) {}

    /** 初始化指定规模的划分结果 */
    void init(int node_cnt, int fpga_cnt);

    /** 将节点 vi 分配到 FPGA f */
    void assign(int vi, int f);

    /** 撤销节点 vi 的分配 */
    void unassign(int vi, int f);
};

// ================================================================
//  5. 多层粗化层级存储
// ================================================================
struct CoarsenLevel {
    CircuitGraph super_graph;                   // 当前层超图
    vector<CandidateSet> super_cddt;            // 超节点候选 FPGA 集合
    unordered_map<int, vector<int>> super2origin; // 超节点 → 原始子节点列表映射
    vector<int> super_weight;                   // 超节点权重（所含原始节点数）
    PartitionResult level_part;                 // 该层划分结果

    CoarsenLevel() = default;
};

// ================================================================
//  6. 算法全局配置参数
// ================================================================
struct TopoPartConfig {
    double beta_ratio = 2.0 / 3.0;      // β 初始占比（β = beta_ratio × fpga_num）
    int    max_coarsen_retry = 5;        // 粗化失败时最大重试次数
    double beta_grow_factor = 1.1;       // 每次重试 β 增大的因子
    int    max_traceback_depth = 100;    // 初始划分最大回溯深度
    bool   verbose = true;               // 是否输出详细日志

    /** 根据 FPGA 数量计算初始 β 阈值 */
    int get_beta(int fpga_num) const {
        return max(1, (int)(beta_ratio * fpga_num));
    }
};

#endif // TYPES_H
