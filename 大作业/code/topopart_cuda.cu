/**
 * @file    topopart_cuda.cu
 * @brief   TopoPart GPU 并行加速 CUDA Kernel 实现
 *
 * 编译要求：CUDA 13.3+, sm_120 (Blackwell RTX 5060), C++17
 *
 * Kernel 设计原则：
 *   1. 粗粒度并行：每个线程处理独立数据单元（节点/边）
 *   2. 避免线程发散：使用统一执行路径，减少分支
 *   3. 共享内存优化：邻域投票使用 __shared__ 暂存
 *   4. 原子操作：用于全局计数器和负载更新
 *   5. Warp-level 原语：在可行场景使用 __shfl_sync 等
 */

#include "topopart_cuda.cuh"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>

// ================================================================
//  错误检查宏
// ================================================================

#define CUDA_CHECK(call) do {                                    \
    cudaError_t err = call;                                      \
    if (err != cudaSuccess) {                                    \
        fprintf(stderr, "CUDA Error at %s:%d — %s\n",           \
                __FILE__, __LINE__, cudaGetErrorString(err));    \
        exit(EXIT_FAILURE);                                      \
    }                                                            \
} while (0)

#define CUDA_CHECK_KERNEL() CUDA_CHECK(cudaGetLastError())

// ================================================================
//  GPU 内存管理实现
// ================================================================

bool gpu_init(GPUContext& ctx, int N, int K, int E) {
    ctx.N = N;
    ctx.K = K;
    ctx.E = E;
    ctx.K_words = (K + 31) / 32;

    try {
        // 图结构
        CUDA_CHECK(cudaMalloc(&ctx.d_node_adj_off, (N + 1) * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_node_adj_list, E * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_circuit_dist, (size_t)N * N * sizeof(int)));

        // FPGA 拓扑
        CUDA_CHECK(cudaMalloc(&ctx.d_fpga_dist, (size_t)K * K * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_fpga_resource_cap, K * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_max_dist, K * sizeof(int)));

        // 候选集
        CUDA_CHECK(cudaMalloc(&ctx.d_cddt_mask, (size_t)N * ctx.K_words * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&ctx.d_cddt_sizes, N * sizeof(int)));

        // 固定节点
        CUDA_CHECK(cudaMalloc(&ctx.d_fixed_mask, N * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_fixed_fpga, N * sizeof(int)));

        // 划分结果
        CUDA_CHECK(cudaMalloc(&ctx.d_partition, N * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_fpga_loads, K * sizeof(int)));

        // 工作缓冲区
        CUDA_CHECK(cudaMalloc(&ctx.d_work_buffer1, N * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&ctx.d_work_buffer2, N * sizeof(int)));

        return true;
    } catch (...) {
        gpu_free(ctx);
        return false;
    }
}

void gpu_free(GPUContext& ctx) {
    auto safeFree = [](int*& ptr) {
        if (ptr) { cudaFree(ptr); ptr = nullptr; }
    };
    safeFree(ctx.d_node_adj_off);
    safeFree(ctx.d_node_adj_list);
    safeFree(ctx.d_circuit_dist);
    safeFree(ctx.d_fpga_dist);
    safeFree(ctx.d_fpga_resource_cap);
    safeFree(ctx.d_max_dist);
    safeFree(ctx.d_cddt_mask);
    safeFree(ctx.d_cddt_sizes);
    safeFree(ctx.d_fixed_mask);
    safeFree(ctx.d_fixed_fpga);
    safeFree(ctx.d_partition);
    safeFree(ctx.d_fpga_loads);
    safeFree(ctx.d_work_buffer1);
    safeFree(ctx.d_work_buffer2);
    safeFree(ctx.d_work_int_NxK);
    safeFree(ctx.d_S_cache_off);
    safeFree(ctx.d_S_cache_data);

    ctx.N = ctx.K = ctx.E = ctx.K_words = 0;
}

bool gpu_upload_graph(GPUContext& ctx,
                      const vector<int>& node_adj_off,
                      const vector<int>& node_adj_list,
                      const vector<int>& circuit_dist_flat,
                      const vector<int>& fixed_mask,
                      const vector<int>& fixed_fpga) {
    CUDA_CHECK(cudaMemcpy(ctx.d_node_adj_off, node_adj_off.data(),
                          (ctx.N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_node_adj_list, node_adj_list.data(),
                          ctx.E * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_circuit_dist, circuit_dist_flat.data(),
                          (size_t)ctx.N * ctx.N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_fixed_mask, fixed_mask.data(),
                          ctx.N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_fixed_fpga, fixed_fpga.data(),
                          ctx.N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaDeviceSynchronize());
    return true;
}

bool gpu_upload_fpga(GPUContext& ctx,
                     const vector<int>& fpga_dist_flat,
                     const vector<int>& resource_cap,
                     const vector<int>& max_dist,
                     const vector<int>& S_cache_off,
                     const vector<int>& S_cache_data) {
    CUDA_CHECK(cudaMemcpy(ctx.d_fpga_dist, fpga_dist_flat.data(),
                          (size_t)ctx.K * ctx.K * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_fpga_resource_cap, resource_cap.data(),
                          ctx.K * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_max_dist, max_dist.data(),
                          ctx.K * sizeof(int), cudaMemcpyHostToDevice));

    if (!S_cache_off.empty()) {
        CUDA_CHECK(cudaMalloc(&ctx.d_S_cache_off, S_cache_off.size() * sizeof(int)));
        CUDA_CHECK(cudaMemcpy(ctx.d_S_cache_off, S_cache_off.data(),
                              S_cache_off.size() * sizeof(int), cudaMemcpyHostToDevice));
    }
    if (!S_cache_data.empty()) {
        CUDA_CHECK(cudaMalloc(&ctx.d_S_cache_data, S_cache_data.size() * sizeof(int)));
        CUDA_CHECK(cudaMemcpy(ctx.d_S_cache_data, S_cache_data.data(),
                              S_cache_data.size() * sizeof(int), cudaMemcpyHostToDevice));
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    return true;
}

bool gpu_upload_cddt(GPUContext& ctx,
                     const vector<uint32_t>& cddt_mask_flat,
                     const vector<int>& cddt_sizes) {
    CUDA_CHECK(cudaMemcpy(ctx.d_cddt_mask, cddt_mask_flat.data(),
                          (size_t)ctx.N * ctx.K_words * sizeof(uint32_t),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx.d_cddt_sizes, cddt_sizes.data(),
                          ctx.N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaDeviceSynchronize());
    return true;
}

bool gpu_download_partition(const GPUContext& ctx, vector<int>& partition) {
    partition.resize(ctx.N);
    CUDA_CHECK(cudaMemcpy(partition.data(), ctx.d_partition,
                          ctx.N * sizeof(int), cudaMemcpyDeviceToHost));
    return true;
}

bool gpu_download_cddt(const GPUContext& ctx,
                       vector<uint32_t>& cddt_mask_flat,
                       vector<int>& cddt_sizes) {
    cddt_mask_flat.resize((size_t)ctx.N * ctx.K_words);
    cddt_sizes.resize(ctx.N);
    CUDA_CHECK(cudaMemcpy(cddt_mask_flat.data(), ctx.d_cddt_mask,
                          (size_t)ctx.N * ctx.K_words * sizeof(uint32_t),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(cddt_sizes.data(), ctx.d_cddt_sizes,
                          ctx.N * sizeof(int), cudaMemcpyDeviceToHost));
    return true;
}

// ================================================================
//  Kernel 1: 并行全点对 BFS 最短路径
// ================================================================

__global__ void kernel_all_pair_bfs(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    int* __restrict__ dist_matrix,
    int N) {

    int src = blockIdx.x;  // 每个 block 处理一个源点

    if (src >= N) return;

    int tid = threadIdx.x;
    int* src_dist = dist_matrix + (size_t)src * N;  // 该源点的距离行

    // 初始化：INF
    if (tid < N) {
        src_dist[tid] = N + 999;  // 大于任何可能距离
    }
    __syncthreads();

    // 线程 0 设置源点距离
    if (tid == 0) {
        src_dist[src] = 0;
    }
    __syncthreads();

    // 迭代 BFS（最坏情况 N 轮，实际图直径远小于 N）
    for (int level = 0; level < N; ++level) {
        bool changed = false;

        // 每个线程处理一段节点
        for (int u = tid; u < N; u += blockDim.x) {
            if (src_dist[u] == level) {
                // 遍历邻居
                int start = node_adj_off[u];
                int end   = node_adj_off[u + 1];
                for (int e = start; e < end; ++e) {
                    int v = node_adj_list[e];
                    // 原子更新（如果当前距离更优）
                    int old = atomicMin(&src_dist[v], level + 1);
                    if (old > level + 1) {
                        changed = true;
                    }
                }
            }
        }

        __syncthreads();

        // 检查是否收敛
        __shared__ bool any_changed;
        if (tid == 0) any_changed = changed;
        __syncthreads();

        if (!any_changed) break;
    }
}

// ================================================================
//  Kernel 2: 并行候选集传播
// ================================================================

__global__ void kernel_candidate_propagate(
    const int* __restrict__ fixed_nodes_list,
    int num_fixed,
    const int* __restrict__ fixed_fpga,
    const int* __restrict__ circuit_dist,
    const int* __restrict__ fpga_dist,
    const int* __restrict__ max_dist,
    const int* __restrict__ S_cache_off,
    const int* __restrict__ S_cache_data,
    uint32_t* __restrict__ cddt_mask,
    int* __restrict__ cddt_sizes,
    int* __restrict__ new_fixed_flags,
    int N, int K, int K_words) {

    int fi = blockIdx.x;   // 固定节点索引
    int vj = threadIdx.x + blockIdx.y * blockDim.x;  // 待检查的可移动节点

    if (fi >= num_fixed || vj >= N) return;

    int vi = fixed_nodes_list[fi];
    int v_hat_i = fixed_fpga[vi];

    // 跳过自身和固定节点
    if (vj == vi) return;
    if (new_fixed_flags[vj]) return;  // 已固定

    // 电路距离 k = dist(vi, vj)
    int k = circuit_dist[(size_t)vi * N + vj];
    int max_d = max_dist[v_hat_i];

    // 只有 k < max_d 才进行收缩（定理 III.1 的传播条件）
    if (k >= max_d) return;

    // 获取 Ŝ(v̂i, k)：通过 S_cache 查找
    // S_cache_off[center * (max_d+1) + x] 给出偏移
    // 简化：直接在 FPGA 距离矩阵中逐位检查
    //
    // 对于 vj 的每个候选 FPGA f：
    //   如果 dist_fpga(v̂i, f) > k → f 不在 Ŝ(v̂i, k) → 从候选集移除

    uint32_t* vj_mask = cddt_mask + (size_t)vj * K_words;

    for (int f_word = 0; f_word < K_words; ++f_word) {
        uint32_t word = vj_mask[f_word];
        if (word == 0) continue;

        uint32_t new_word = word;
        int base = f_word * 32;

        // 逐位检查该 word 中的每个 FPGA
        for (int bit = 0; bit < 32 && (base + bit) < K; ++bit) {
            if (!(word & (1u << bit))) continue;

            int f = base + bit;
            int fpga_d = fpga_dist[(size_t)v_hat_i * K + f];

            if (fpga_d > k) {
                // 不在 Ŝ(v̂i, k) 中 → 清除该位
                new_word &= ~(1u << bit);
            }
        }

        if (new_word != word) {
            vj_mask[f_word] = new_word;
        }
    }

    // 重新计数候选集大小
    int new_size = 0;
    for (int f_word = 0; f_word < K_words; ++f_word) {
        new_size += __popc(vj_mask[f_word]);
    }

    if (new_size < cddt_sizes[vj]) {
        cddt_sizes[vj] = new_size;

        // 若收缩为单元素 → 标记为新固定节点
        if (new_size == 1 && new_size > 0) {
            new_fixed_flags[vj] = 1;

            // 找出唯一的 FPGA
            for (int f_word = 0; f_word < K_words; ++f_word) {
                uint32_t w = vj_mask[f_word];
                if (w != 0) {
                    int only_f = f_word * 32 + __ffs(w) - 1;
                    fixed_fpga[vj] = only_f;  // 注意：这里 fixed_fpga 实际存储绑定FPGA
                    break;
                }
            }
        }
    }
}

// ================================================================
//  Kernel 3: 并行匹配评估
// ================================================================

__global__ void kernel_evaluate_match(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const uint32_t* __restrict__ cddt_mask,
    int* __restrict__ match_scores,
    int* __restrict__ match_partner,
    int N, int K_words, int beta) {

    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= node_adj_off[N]) return;

    // 从扁平的邻接表中恢复 (u, v)
    // 由于每条无向边出现两次，需要去重逻辑
    // 简化：遍历所有节点，每个线程处理一个节点的第一个邻居

    // 查找 e 对应的节点 u
    int u = 0;
    while (u < N && e >= node_adj_off[u + 1]) ++u;
    if (u >= N) return;

    int local_e = e - node_adj_off[u];
    int v = node_adj_list[e];

    if (v < u) return;  // 仅处理 u < v 的边（去重）

    // 计算候选集交集大小
    const uint32_t* mask_u = cddt_mask + (size_t)u * K_words;
    const uint32_t* mask_v = cddt_mask + (size_t)v * K_words;

    int intersect_cnt = 0;
    for (int w = 0; w < K_words; ++w) {
        intersect_cnt += __popc(mask_u[w] & mask_v[w]);
    }

    // 条件检查：交集大小 ≥ β
    int score = (intersect_cnt >= beta) ? intersect_cnt : -1;
    match_scores[e] = score;

    // 更新最佳匹配（使用原子操作）
    if (score > 0) {
        atomicMax(&match_partner[u], score * N + v);  // 编码：(score, v)
    }
}

// ================================================================
//  Kernel 4: 并行标签传播初始化
// ================================================================

__global__ void kernel_label_init(
    int* __restrict__ partition,
    const int* __restrict__ fixed_mask,
    const int* __restrict__ fixed_fpga,
    int N, int K, unsigned int seed) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    if (fixed_mask[i]) {
        partition[i] = fixed_fpga[i];
    } else {
        // 简单哈希随机初始化
        unsigned int hash = (i * 2654435761u) ^ seed;
        partition[i] = hash % K;
    }
}

// ================================================================
//  Kernel 5: 并行标签传播更新
// ================================================================

__global__ void kernel_label_update(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    const int* __restrict__ fixed_mask,
    const uint32_t* __restrict__ cddt_mask,
    int* __restrict__ new_partition,
    const int* __restrict__ fpga_loads,
    const int* __restrict__ resource_cap,
    int N, int K, int K_words) {

    extern __shared__ int local_votes[];  // [blockDim.x][K] 实际分配

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    if (i >= N) return;

    // 固定节点不变
    if (fixed_mask[i]) {
        new_partition[i] = partition[i];
        return;
    }

    // 投票数组清零（shared memory）
    int* my_votes = local_votes + tid * K;
    for (int f = 0; f < K; ++f) {
        my_votes[f] = 0;
    }

    // 遍历邻居，收集 FPGA 投票
    int start = node_adj_off[i];
    int end   = node_adj_off[i + 1];
    for (int e = start; e < end; ++e) {
        int nb = node_adj_list[e];
        int nb_fpga = partition[nb];
        if (nb_fpga >= 0 && nb_fpga < K) {
            my_votes[nb_fpga]++;
        }
    }

    // 选择得票最高的合法 FPGA
    // 合法性：在 Cddt[i] 中 + 资源容量不超限
    int best_f = partition[i];  // 默认保持当前
    int best_votes = 0;

    const uint32_t* my_mask = cddt_mask + (size_t)i * K_words;

    for (int f = 0; f < K; ++f) {
        // 检查是否在候选集中
        int word_idx = f / 32;
        int bit = f % 32;
        if (!(my_mask[word_idx] & (1u << bit))) continue;

        // 资源检查：目标 FPGA 还有空间
        // （简化：此处用全局 load 近似，精确实现需 atomicAdd）
        if (resource_cap[f] > 0 && fpga_loads[f] >= resource_cap[f]) continue;

        if (my_votes[f] > best_votes) {
            best_votes = my_votes[f];
            best_f = f;
        }
    }

    new_partition[i] = best_f;
}

// ================================================================
//  Kernel 6: 并行 FPGA 负载统计
// ================================================================

__global__ void kernel_compute_loads(
    const int* __restrict__ partition,
    int* __restrict__ fpga_loads,
    int N, int K) {

    // 初始化
    int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f < K) {
        fpga_loads[f] = 0;
    }
    __syncthreads();

    // 每个线程统计一段节点的负载
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;

    int p = partition[tid];
    if (p >= 0 && p < K) {
        atomicAdd(&fpga_loads[p], 1);
    }
}

// ================================================================
//  Kernel 7: 并行边界节点检测
// ================================================================

__global__ void kernel_detect_boundary(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    int* __restrict__ boundary_flags,
    int N) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    int my_fpga = partition[i];
    if (my_fpga < 0) {
        boundary_flags[i] = 0;
        return;
    }

    int start = node_adj_off[i];
    int end   = node_adj_off[i + 1];
    for (int e = start; e < end; ++e) {
        int nb = node_adj_list[e];
        int nb_fpga = partition[nb];
        if (nb_fpga >= 0 && nb_fpga != my_fpga) {
            boundary_flags[i] = 1;
            return;
        }
    }
    boundary_flags[i] = 0;
}

// ================================================================
//  Kernel 8: 并行 FM 增益评估
// ================================================================

__global__ void kernel_fm_evaluate_gains(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    const uint32_t* __restrict__ cddt_mask,
    const int* __restrict__ boundary_list,
    int num_boundary,
    const int* __restrict__ fpga_dist,
    const int* __restrict__ circuit_dist,
    int* __restrict__ gains,
    int N, int K, int K_words) {

    int bi = blockIdx.x;           // 边界节点索引
    int f  = threadIdx.x;          // 目标 FPGA

    if (bi >= num_boundary || f >= K) return;

    int vi = boundary_list[bi];

    // 检查 f 是否在 vi 的候选集中
    const uint32_t* my_mask = cddt_mask + (size_t)vi * K_words;
    int word_idx = f / 32;
    int bit = f % 32;
    if (!(my_mask[word_idx] & (1u << bit))) {
        gains[(size_t)bi * K + f] = INT_MAX;  // 不可行
        return;
    }

    int old_f = partition[vi];
    if (old_f == f) {
        gains[(size_t)bi * K + f] = 0;
        return;
    }

    // 计算割边变化（简化：只计跨 FPGA 边数变化）
    int delta = 0;
    int start = node_adj_off[vi];
    int end   = node_adj_off[vi + 1];

    for (int e = start; e < end; ++e) {
        int nb = node_adj_list[e];
        int nb_f = partition[nb];
        if (nb_f < 0) continue;

        bool was_cut = (old_f != nb_f);
        bool new_cut = (f != nb_f);

        if (!was_cut && new_cut) delta += 1;
        else if (was_cut && !new_cut) delta -= 1;
    }

    // 拓扑约束检查：定理 III.1
    for (int e = start; e < end; ++e) {
        int nb = node_adj_list[e];
        int nb_f = partition[nb];
        if (nb_f < 0) continue;

        int circuit_d = circuit_dist[(size_t)vi * N + nb];
        int fpga_d = fpga_dist[(size_t)f * K + nb_f];

        if (circuit_d < fpga_d) {
            // 违反 x ≥ y → 不可行
            gains[(size_t)bi * K + f] = INT_MAX;
            return;
        }
    }

    gains[(size_t)bi * K + f] = delta;
}

// ================================================================
//  Kernel 9: 并行拓扑违规检查
// ================================================================

__global__ void kernel_check_topology(
    const int* __restrict__ node_adj_off,
    const int* __restrict__ node_adj_list,
    const int* __restrict__ partition,
    const int* __restrict__ fpga_dist,
    const int* __restrict__ circuit_dist,
    int* __restrict__ violation_flags,
    int N, int K) {

    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= node_adj_off[N]) return;

    // 解析边 (u, v)
    int u = 0;
    while (u < N && e >= node_adj_off[u + 1]) ++u;
    if (u >= N) return;
    int v = node_adj_list[e];

    if (v <= u) {
        violation_flags[e] = 0;
        return;  // 每条无向边只检查一次
    }

    int fu = partition[u];
    int fv = partition[v];

    if (fu < 0 || fv < 0 || fu == fv) {
        violation_flags[e] = 0;
        return;
    }

    // 定理 III.1：需要 circuit_dist ≥ fpga_dist
    // 论文硬约束：不同 FPGA 之间必须有直连边（dist == 1）
    int fpga_d = fpga_dist[(size_t)fu * K + fv];
    violation_flags[e] = (fpga_d > 1) ? 1 : 0;
}

// ================================================================
//  GPU 高级 API 实现
// ================================================================

bool gpu_candidate_propagation(GPUContext& ctx,
                                const vector<int>& fixed_seeds,
                                int max_rounds) {
    int N = ctx.N, K = ctx.K;
    int num_fixed = (int)fixed_seeds.size();

    if (num_fixed == 0) return true;

    // 上传固定种子列表
    int* d_fixed_list;
    CUDA_CHECK(cudaMalloc(&d_fixed_list, num_fixed * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_fixed_list, fixed_seeds.data(),
                          num_fixed * sizeof(int), cudaMemcpyHostToDevice));

    // 新固定标记数组
    CUDA_CHECK(cudaMemset(ctx.d_work_buffer1, 0, N * sizeof(int)));

    // 网格配置
    dim3 block(32, 1);
    dim3 grid(num_fixed, (N + 31) / 32);

    for (int round = 0; round < max_rounds; ++round) {
        kernel_candidate_propagate<<<grid, block>>>(
            d_fixed_list, num_fixed,
            ctx.d_fixed_fpga,
            ctx.d_circuit_dist,
            ctx.d_fpga_dist,
            ctx.d_max_dist,
            ctx.d_S_cache_off,
            ctx.d_S_cache_data,
            (uint32_t*)ctx.d_cddt_mask,
            ctx.d_cddt_sizes,
            ctx.d_work_buffer1,  // new_fixed_flags
            N, K, ctx.K_words);
        CUDA_CHECK_KERNEL();

        // 检查是否有新固定节点（简化：轮询收敛）
        // 完整实现应使用 device 端的规约判断收敛
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaFree(d_fixed_list));
    return true;
}

bool gpu_coarsen_matching(GPUContext& ctx, int beta,
                           vector<int>& matches) {
    int N = ctx.N;

    // 清零匹配分数和伙伴数组
    CUDA_CHECK(cudaMemset(ctx.d_work_buffer1, 0, ctx.E * sizeof(int)));
    CUDA_CHECK(cudaMemset(ctx.d_work_buffer2, -1, N * sizeof(int)));

    dim3 block(256);
    dim3 grid((ctx.E + 255) / 256);

    kernel_evaluate_match<<<grid, block>>>(
        ctx.d_node_adj_off,
        ctx.d_node_adj_list,
        (const uint32_t*)ctx.d_cddt_mask,
        ctx.d_work_buffer1,   // match_scores
        ctx.d_work_buffer2,   // match_partner
        N, ctx.K_words, beta);
    CUDA_CHECK_KERNEL();
    CUDA_CHECK(cudaDeviceSynchronize());

    // 下载匹配结果
    matches.resize(N, -1);
    CUDA_CHECK(cudaMemcpy(matches.data(), ctx.d_work_buffer2,
                          N * sizeof(int), cudaMemcpyDeviceToHost));

    // 解码：match_partner[u] = score * N + v → v = match_partner[u] % N
    for (int u = 0; u < N; ++u) {
        if (matches[u] >= 0) {
            matches[u] = matches[u] % N;
        }
    }

    return true;
}

bool gpu_label_propagation(GPUContext& ctx, int max_iters) {
    int N = ctx.N, K = ctx.K;

    // 随机初始化
    {
        dim3 block(256);
        dim3 grid((N + 255) / 256);
        kernel_label_init<<<grid, block>>>(
            ctx.d_partition,
            ctx.d_fixed_mask,
            ctx.d_fixed_fpga,
            N, K, 42u);
        CUDA_CHECK_KERNEL();
    }

    // 分配 new_partition 缓冲区
    int* d_new_partition;
    CUDA_CHECK(cudaMalloc(&d_new_partition, N * sizeof(int)));

    // 共享内存大小：K * sizeof(int) * blockDim.x
    size_t shared_bytes = (size_t)K * 256 * sizeof(int);

    for (int iter = 0; iter < max_iters; ++iter) {
        // 更新负载
        {
            dim3 block(256);
            dim3 grid((K + 255) / 256);
            kernel_compute_loads<<<grid, block>>>(
                ctx.d_partition, ctx.d_fpga_loads, N, K);
            CUDA_CHECK_KERNEL();
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        // 标签传播更新
        {
            dim3 block(256);
            dim3 grid((N + 255) / 256);
            kernel_label_update<<<grid, block, shared_bytes>>>(
                ctx.d_node_adj_off,
                ctx.d_node_adj_list,
                ctx.d_partition,
                ctx.d_fixed_mask,
                (const uint32_t*)ctx.d_cddt_mask,
                d_new_partition,
                ctx.d_fpga_loads,
                ctx.d_fpga_resource_cap,
                N, K, ctx.K_words);
            CUDA_CHECK_KERNEL();
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        // 交换 buffer
        swap(ctx.d_partition, d_new_partition);
    }

    CUDA_CHECK(cudaFree(d_new_partition));
    return true;
}

bool gpu_fm_refinement(GPUContext& ctx, int max_passes) {
    int N = ctx.N, K = ctx.K;

    // 边界标记数组
    int* d_boundary;
    CUDA_CHECK(cudaMalloc(&d_boundary, N * sizeof(int)));

    // 增益数组（动态分配）
    int* d_gains;
    CUDA_CHECK(cudaMalloc(&d_gains, (size_t)N * K * sizeof(int)));

    for (int pass = 0; pass < max_passes; ++pass) {
        // 更新负载
        {
            dim3 block(256);
            dim3 grid((K + 255) / 256);
            kernel_compute_loads<<<grid, block>>>(
                ctx.d_partition, ctx.d_fpga_loads, N, K);
            CUDA_CHECK_KERNEL();
        }

        // 检测边界
        {
            dim3 block(256);
            dim3 grid((N + 255) / 256);
            kernel_detect_boundary<<<grid, block>>>(
                ctx.d_node_adj_off,
                ctx.d_node_adj_list,
                ctx.d_partition,
                d_boundary, N);
            CUDA_CHECK_KERNEL();
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        // 简化：此处仅评估增益，实际移动需要在 CPU 侧决策
        // 完整实现需要 device 端原子操作协调 FM pass
        {
            dim3 block(min(K, 32), 4);
            dim3 grid(N, 1);
            kernel_fm_evaluate_gains<<<grid, block>>>(
                ctx.d_node_adj_off,
                ctx.d_node_adj_list,
                ctx.d_partition,
                (const uint32_t*)ctx.d_cddt_mask,
                d_boundary, N,
                ctx.d_fpga_dist,
                ctx.d_circuit_dist,
                d_gains,
                N, K, ctx.K_words);
            CUDA_CHECK_KERNEL();
        }
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    CUDA_CHECK(cudaFree(d_boundary));
    CUDA_CHECK(cudaFree(d_gains));
    return true;
}

bool gpu_check_topology(const GPUContext& ctx, int& violation_count) {
    int N = ctx.N, K = ctx.K;
    int total_edges = ctx.E;

    int* d_violations;
    CUDA_CHECK(cudaMalloc(&d_violations, total_edges * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_violations, 0, total_edges * sizeof(int)));

    dim3 block(256);
    dim3 grid((total_edges + 255) / 256);

    kernel_check_topology<<<grid, block>>>(
        ctx.d_node_adj_off,
        ctx.d_node_adj_list,
        ctx.d_partition,
        ctx.d_fpga_dist,
        ctx.d_circuit_dist,
        d_violations,
        N, K);
    CUDA_CHECK_KERNEL();
    CUDA_CHECK(cudaDeviceSynchronize());

    // 规约求和
    vector<int> violations(total_edges);
    CUDA_CHECK(cudaMemcpy(violations.data(), d_violations,
                          total_edges * sizeof(int), cudaMemcpyDeviceToHost));

    violation_count = 0;
    for (int v : violations) violation_count += v;

    CUDA_CHECK(cudaFree(d_violations));
    return true;
}
