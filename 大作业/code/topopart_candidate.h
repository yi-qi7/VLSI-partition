#ifndef TOPOPART_CANDIDATE_H
#define TOPOPART_CANDIDATE_H

/**
 * @file    topopart_candidate.h
 * @brief   TopoPart 模块2：候选 FPGA 传播 CandidateFPGAPropagation（Algorithm 1）
 *
 * 论文 Algorithm 1: Candidate FPGA Set Propagation
 *
 * 输入：原始电路图 G、FPGA 拓扑图 Ĝ
 * 输出：每个电路节点初始化完成的候选集 Cddt[] 数组
 *
 * 核心思想：
 *   从固定节点出发，利用定理 III.1 的「距离下界」逐层收缩
 *   可移动节点的候选 FPGA 集合，直至所有节点候选集稳定。
 *
 * 定理 III.1（判别依据）：
 *   电路节点 vi, vj 最短电路距离 x；
 *   FPGA v̂i, v̂j 硬件最短距离 y；
 *   vi 分配 v̂i、vj 分配 v̂j 无拓扑违规 ⇔ x ≥ y。
 *
 * 传播逻辑：
 *   对固定节点 vi（已绑定 v̂i），考察距离 < maxDist(v̂i) 的可移动节点 vj：
 *     令 k = dist(vi, vj)（电路距离）
 *     则 vj 只能分配在 Ŝ(v̂i, k) 内（距离 v̂i ≤ k 的 FPGA），
 *     否则将违反定理 III.1。
 *   即：Cddt[vj] = Cddt[vj] ∩ Ŝ(v̂i, k)
 */

#include "topopart_types.h"
#include <queue>
#include <utility>

/**
 * @brief 候选 FPGA 传播主函数
 * @param g   原始电路图（需预先计算 dist_all）
 * @param fg  FPGA 拓扑图（需预先计算 dist_all / max_dist / S_cache）
 * @return 每个电路节点的初始化候选集 Cddt 数组（长度 = g.node_num）
 *
 * @throws runtime_error 若任意节点候选集为空（无可行拓扑划分解）
 *
 * 内部执行步骤（严格按照 Algorithm 1）：
 *   1. 预计算 FPGA 图全点对距离 dist_all、每个 FPGA max_dist；预缓存所有 Ŝ(v̂i, x)；
 *      —— 由调用者在调用前通过 BuildFPGADist(fg) 完成
 *   2. 初始化候选集：
 *      - 固定节点仅保留绑定 FPGA
 *      - 可移动节点候选集 = 全部 FPGA
 *   3. 队列 queue<pair<int,int>> q：存储 (固定节点 id, 绑定 FPGA id)，所有固定节点入队
 *   4. 循环弹出队列元素 (vi, v̂i)：
 *      a. d = fg.max_dist[v̂i]；构建 S(vi, d)：电路上距离 vi < d 的所有可移动节点 vj
 *      b. k = g.dist_all[vi][vj]（电路距离）；取出 Ŝ(v̂i, k)（FPGA 距离 ≤ k 的 FPGA 集）
 *      c. Cddt[vj].fpga_set = Cddt[vj].fpga_set ∩ Ŝ(v̂i, k)；同步更新 T_vec 计数
 *      d. 分支判定：
 *         - 若 Cddt[vj].fpga_set.size() == 1：vj 转为新固定节点，(vj, 唯一FPGA) 入队
 *         - 若 Cddt[vj].fpga_set.empty()：抛出异常，无可行拓扑划分解
 *   5. 返回完整候选集数组
 */
vector<CandidateSet> CandidateFPGAPropagation(CircuitGraph& g, FPGAGraph& fg);

#endif // TOPOPART_CANDIDATE_H
