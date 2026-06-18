#ifndef TOPOPART_CUDA_CUH
#define TOPOPART_CUDA_CUH

/**
 * @file    topopart_cuda.cuh
 * @brief   TopoPart GPU 并行加速模块（CUDA C++17）
 *
 * 目标硬件：NVIDIA RTX 5060 (Blackwell, sm_120), CUDA 13.3
 *
 * GPU 并行化的核心思路：
 *   1. 全点对最短路径：每个源点的 Dijkstra → GPU 多源并行 BFS
 *   2. 候选集传播：每个固定节点的传播 → GPU 线程并行收缩候选集
 *   3. 粗化匹配：节点对的合并判定 → GPU 线程并行评估 + 原子操作选最优
 *   4. 初始划分：标签传播 → GPU 线程并行邻域投票
 *   5. 细化：边界节点 FM → GPU 线程并行评估迁移增益
 *   6. 拓扑检查：定理 III.1 验证 → GPU 线程并行检查每条边
 *
 * 数据布局（GPU 侧，CSR 格式）：
 *   - node_adj_off[0..N]：节点邻接偏移
 *   - node_adj_list[0..2E-1]：邻接表扁平化
 *   - fpga_dist[0..K*K-1]：FPGA 距离矩阵（行优先）
 *   - fpga_S_cache：Ŝ(v̂i, x) 预计算缓存（CSR 格式）
 *   - cddt_mask[0..N*K-1]：候选集位掩码（每个节点 K 位）
 *   - partition[0..N-1]：当前划分
 */

#include <cuda_runtime.h>
#include <vector>
#include <cstdint>

using namespace std;

// ================================================================
//  GPU 侧数据结构（Device 内存）
// ================================================================

struct GPUContext {
    // 图结构 (CSR 格式)
    int N;                          // 节点数
    int K;                          // FPGA 数
    int E;                          // 边数（无向边 × 2）

    int* d_node_adj_off;            // [N+1] 邻接表偏移
    int* d_node_adj_list;           // [2E] 邻接表
    int* d_circuit_dist;            // [N*N] 电路距离矩阵（行优先）

    // FPGA 拓扑
    int* d_fpga_dist;               // [K*K] FPGA 距离矩阵（行优先）
    int* d_fpga_resource_cap;       // [K] 资源容量

    // Ŝ(v̂i, x) 缓存（CSR 格式）
    int* d_S_cache_off;             // [K*(max_dist+2)] 每个 (center,x) 的偏移
    int* d_S_cache_data;            // 扁平化的 FPGA 集合
    int* d_max_dist;                // [K] maxDist(v̂i)

    // 候选集（位掩码：每个节点 K 个 int32，每 32 FPGA 一组）
    int  K_words;                   // 候选集位掩码宽度 = (K + 31) / 32
    int* d_cddt_mask;               // [N * K_words] 候选集位掩码

    // 候选集计数（用于 T_vec 等效）
    int* d_cddt_sizes;              // [N] 每个节点候选集大小

    // 固定节点
    int* d_fixed_mask;              // [N] 1=固定, 0=可移动
    int* d_fixed_fpga;              // [N] 固定绑定的 FPGA（仅 fixed_mask[i]=1 时有效）

    // 划分结果
    int* d_partition;               // [N] partition[i] = FPGA id
    int* d_fpga_loads;              // [K] 每个 FPGA 当前负载

    // 工作缓冲区
    int* d_work_buffer1;            // [N] 通用工作缓冲区
    int* d_work_buffer2;            // [N]
    int* d_work_int_NxK;            // [N*K]

    // 常量内存（通过 __constant__ 或 kernel 参数传递）
    // （CUDA 中通过 kernel 参数传递标量即可）

    GPUContext() : N(0), K(0), E(0), K_words(0),
                   d_node_adj_off(nullptr), d_node_adj_list(nullptr),
                   d_circuit_dist(nullptr), d_fpga_dist(nullptr),
                   d_fpga_resource_cap(nullptr), d_S_cache_off(nullptr),
                   d_S_cache_data(nullptr), d_max_dist(nullptr),
                   d_cddt_mask(nullptr), d_cddt_sizes(nullptr),
                   d_fixed_mask(nullptr), d_fixed_fpga(nullptr),
                   d_partition(nullptr), d_fpga_loads(nullptr),
                   d_work_buffer1(nullptr), d_work_buffer2(nullptr),
                   d_work_int_NxK(nullptr) {}
};

// ================================================================
//  GPU 内存管理 API
// ================================================================

/** 初始化 GPU 上下文并分配设备内存 */
bool gpu_init(GPUContext& ctx, int N, int K, int E);

/** 释放 GPU 上下文所有设备内存 */
void gpu_free(GPUContext& ctx);

/** 上传电路图数据到 GPU */
bool gpu_upload_graph(GPUContext& ctx,
                      const vector<int>& node_adj_off,
                      const vector<int>& node_adj_list,
                      const vector<int>& circuit_dist_flat,
                      const vector<int>& fixed_mask,
                      const vector<int>& fixed_fpga);

/** 上传 FPGA 拓扑数据到 GPU */
bool gpu_upload_fpga(GPUContext& ctx,
                     const vector<int>& fpga_dist_flat,
                     const vector<int>& resource_cap,
                     const vector<int>& max_dist,
                     const vector<int>& S_cache_off,
                     const vector<int>& S_cache_data);

/** 上传候选集数据到 GPU */
bool gpu_upload_cddt(GPUContext& ctx,
                     const vector<uint32_t>& cddt_mask_flat,
                     const vector<int>& cddt_sizes);

/** 下载划分结果从 GPU */
bool gpu_download_partition(const GPUContext& ctx, vector<int>& partition);

/** 下载候选集从 GPU */
bool gpu_download_cddt(const GPUContext& ctx,
                       vector<uint32_t>& cddt_mask_flat,
                       vector<int>& cddt_sizes);

// ================================================================
//  GPU Kernel 声明
// ================================================================

/**
 * Kernel 1: 并行全点对最短路径
 *   - 使用 BFS 从每个源点并行
 *   - 网格: (N, 1), 块: (256, 1)
 *   - 每个线程块处理一个源点
 */
__global__ void kernel_all_pair_bfs(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    int* __restrict__ dist_matrix,
    int N);

/**
 * Kernel 2: 并行候选集传播（单轮）
 *   对每个固定节点 vi，传播约束到其邻域
 *   网格: (num_fixed_nodes, N), 块: (32, 1)
 *   每个线程处理一对 (固定节点 vi, 可移动节点 vj) 的收缩
 */
__global__ void kernel_candidate_propagate(
    const int* __restrict__ fixed_nodes_list,
    int num_fixed,
    const int* __restrict__ fixed_fpga,
    const int* __restrict__ circuit_dist,       // [N*N]
    const int* __restrict__ fpga_dist,           // [K*K]
    const int* __restrict__ max_dist,
    const int* __restrict__ S_cache_off,
    const int* __restrict__ S_cache_data,
    uint32_t* __restrict__ cddt_mask,            // [N * K_words]
    int* __restrict__ cddt_sizes,
    int* __restrict__ new_fixed_flags,           // [N] 新固定标记
    int N, int K, int K_words);

/**
 * Kernel 3: 并行匹配评估（粗化阶段）
 *   每个线程处理一条边，评估合并质量（候选集交集大小）
 *   网格: (E, 1), 块: (256, 1)
 */
__global__ void kernel_evaluate_match(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const uint32_t* __restrict__ cddt_mask,
    int* __restrict__ match_scores,              // [E] 匹配得分
    int* __restrict__ match_partner,             // [N] 最佳匹配伙伴
    int N, int K_words, int beta);

/**
 * Kernel 4: 并行标签传播初始化
 *   随机初始化 + 固定节点覆盖
 *   网格: (N, 1), 块: (256, 1)
 */
__global__ void kernel_label_init(
    int* __restrict__ partition,
    const int* __restrict__ fixed_mask,
    const int* __restrict__ fixed_fpga,
    int N, int K, unsigned int seed);

/**
 * Kernel 5: 并行标签传播更新
 *   邻域投票：每个节点选择邻居中出现最多的 FPGA
 *   网格: (N, 1), 块: (256, 1)
 */
__global__ void kernel_label_update(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    const int* __restrict__ fixed_mask,
    const uint32_t* __restrict__ cddt_mask,
    int* __restrict__ new_partition,
    const int* __restrict__ fpga_loads,
    const int* __restrict__ resource_cap,
    int N, int K, int K_words);

/**
 * Kernel 6: 并行 FPGA 负载统计
 *   网格: (K, 1), 块: (256, 1)
 */
__global__ void kernel_compute_loads(
    const int* __restrict__ partition,
    int* __restrict__ fpga_loads,
    int N, int K);

/**
 * Kernel 7: 并行边界节点检测
 *   网格: (N, 1), 块: (256, 1)
 */
__global__ void kernel_detect_boundary(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    int* __restrict__ boundary_flags,
    int N);

/**
 * Kernel 8: 并行 FM 增益评估
 *   每个线程评估一个边界节点迁移到其候选 FPGA 的割边增益
 *   网格: (num_boundary, K), 块: (32, 4)
 */
__global__ void kernel_fm_evaluate_gains(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    const uint32_t* __restrict__ cddt_mask,
    const int* __restrict__ boundary_list,
    int num_boundary,
    const int* __restrict__ fpga_dist,
    const int* __restrict__ circuit_dist,
    int* __restrict__ gains,                     // [num_boundary * K]
    int N, int K, int K_words);

/**
 * Kernel 9: 并行拓扑违规检查
 *   检查每条边是否满足定理 III.1
 *   网格: (E, 1), 块: (256, 1)
 */
__global__ void kernel_check_topology(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    const int* __restrict__ fpga_dist,
    const int* __restrict__ circuit_dist,
    int* __restrict__ violation_flags,            // [E] 违规标记
    int N, int K);

// ================================================================
//  GPU 高级 API（封装多 Kernel 调用）
// ================================================================

/**
 * @brief GPU 加速的候选 FPGA 传播
 * @param ctx          GPU 上下文
 * @param fixed_seeds  固定节点列表
 * @param max_rounds   最大传播轮数
 * @return 成功标志
 */
bool gpu_candidate_propagation(GPUContext& ctx,
                                const vector<int>& fixed_seeds,
                                int max_rounds = 50);

/**
 * @brief GPU 加速的并行匹配粗化
 * @param ctx      GPU 上下文
 * @param beta     合并阈值
 * @param matches  输出：match[u] = v 表示 u,v 匹配
 */
bool gpu_coarsen_matching(GPUContext& ctx, int beta,
                           vector<int>& matches);

/**
 * @brief GPU 加速的标签传播初始划分
 * @param ctx        GPU 上下文
 * @param max_iters  最大迭代轮数
 */
bool gpu_label_propagation(GPUContext& ctx, int max_iters = 20);

/**
 * @brief GPU 加速的 FM 细化
 * @param ctx        GPU 上下文
 * @param max_passes 最大 FM 趟数
 */
bool gpu_fm_refinement(GPUContext& ctx, int max_passes = 10);

/**
 * @brief GPU 加速的拓扑违规检查
 * @param ctx              GPU 上下文
 * @param violation_count  输出：违规数
 */
bool gpu_check_topology(const GPUContext& ctx, int& violation_count);

#endif // TOPOPART_CUDA_CUH
