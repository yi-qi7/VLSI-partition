#ifndef TOPOPART_COARSEN_H
#define TOPOPART_COARSEN_H

/**
 * @file    topopart_coarsen.h
 * @brief   TopoPart 模块3：候选集感知多层粗化 Coarsening
 *
 * 关键参数：
 *   β = 2/3 * fpga_num（初始合并阈值）
 *
 * 合并规则（双条件同时满足才允许合并两节点/超节点）：
 *   条件1：两节点电路互连权重高（邻接关系强）；
 *   条件2：两节点候选集交集大小 ≥ β；
 *
 * 超节点候选集规则：
 *   super_cddt = cddt_a ∩ cddt_b
 *
 * 迭代粗化流程：
 *   循环匹配节点生成超节点，每轮将图规模压缩，
 *   直到无法再合并（图规模不再缩小）或达到最小规模。
 *
 * 容错逻辑：
 *   若后续 InitialPartition 无法找到可行解，
 *   自动增大 β 阈值重新执行 Coarsening。
 */

#include "topopart_types.h"
#include <list>

/**
 * @brief 候选集感知多层粗化
 * @param origin_g    原始电路图
 * @param node_cddt   原始节点的候选 FPGA 集合数组
 * @param fg          FPGA 拓扑图
 * @param beta_init   初始 β 阈值（默认 = 2/3 * fpga_num）
 * @return 多层粗化层级列表 levels，levels.back() 为最粗粒度超图
 *
 * 层级结构：
 *   levels[0]：原始层（细粒度）
 *   levels[1]：第 1 轮粗化
 *   ...
 *   levels[L-1]：最粗粒度超图
 *
 * 每层粗化步骤：
 *   1. 遍历当前层所有节点对，检查合并双条件
 *   2. 贪心匹配：优先合并候选集交集大、互连权重高的节点对
 *   3. 构建超节点图：更新邻接关系、候选集、超节点→子节点映射
 *   4. 若本轮无合并 → 粗化终止
 */
vector<CoarsenLevel> Coarsening(const CircuitGraph& origin_g,
                                 const vector<CandidateSet>& node_cddt,
                                 const FPGAGraph& fg,
                                 int beta_init);

/**
 * @brief 重试粗化（增大 β 阈值）
 * @param origin_g    原始电路图
 * @param node_cddt   原始节点候选集
 * @param fg          FPGA 拓扑图
 * @param beta_init   初始 β 阈值
 * @param max_retry   最大重试次数
 * @param grow_factor 每次重试 β 增大因子
 * @return 多层粗化层级列表
 *
 * 若 InitialPartition 对当前粗化结果无法找到可行解，
 * 调用此函数增大 β（更严格合并条件 → 更少合并 → 更大解空间）。
 */
vector<CoarsenLevel> CoarseningWithRetry(const CircuitGraph& origin_g,
                                          const vector<CandidateSet>& node_cddt,
                                          const FPGAGraph& fg,
                                          const TopoPartConfig& cfg);

#endif // TOPOPART_COARSEN_H
