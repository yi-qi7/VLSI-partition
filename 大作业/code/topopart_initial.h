#ifndef TOPOPART_INITIAL_H
#define TOPOPART_INITIAL_H

/**
 * @file    topopart_initial.h
 * @brief   TopoPart 模块4：拓扑驱动贪心初始划分 InitialPartition（Algorithm 2，带回溯 Traceback）
 *
 * 论文 Algorithm 2: Topology-Driven Greedy Initial Partitioning with Traceback
 *
 * 输入：最粗层超图层级、FPGA 资源拓扑图
 * 输出：粗粒度零违规初始划分结果
 *
 * 内部流程：
 *   1. 初始化：固定超节点直接绑定指定 FPGA；可移动节点标记 node2fpga = -1
 *   2. 双优先队列：
 *      Q：priority_queue<(候选集大小, 节点 id), greater> —— 候选越少越优先分配
 *      R：vector<priority_queue<(割边增量, FPGA id), greater>> —— 割边增量小优先
 *   3. 主分配循环：
 *      a. 弹出 Q 顶部节点 vj'，取出 R[j'] 最优 FPGA v̂j；临时分配 part(vj') = v̂j
 *      b. 遍历所有未分配邻节点 vk'：移除 vk' 候选集中与 v̂j 无直连的 FPGA，更新 T_vec
 *      c. 冲突检测：若任意 vk' 候选集为空 → Traceback = true，进入回溯
 *   4. 回溯恢复逻辑（Traceback）：
 *      a. 撤销 vj' 分配，将 v̂j 永久从 Cddt[vj'] 中删除
 *      b. 恢复所有邻节点 vk' 的 T_vec 计数、候选集
 *      c. vj' 重新加入分配队列 Q，重新选择其他合法 FPGA
 *   5. 分配完成后校验：资源不超限、拓扑违规数 = 0，返回划分结果
 */

#include "topopart_types.h"
#include <queue>
#include <stack>
#include <string>

/**
 * @brief 初始划分（双模式）
 * @param mode  "paper" = Algorithm 2+traceback | "grow" = region-growing+spread
 */
PartitionResult InitialPartition(CoarsenLevel& coarsest_level,
                                  const FPGAGraph& fg,
                                  int max_traceback = 100,
                                  const std::string& mode = "paper");

/**
 * @brief 割边增量评估：计算将节点 vi 分配到 FPGA f 后的割边变化
 * @param level    当前粗化层级
 * @param fg       FPGA 拓扑图
 * @param part     当前划分状态
 * @param vi       待分配节点 id
 * @param target_f 目标 FPGA id
 * @return 割边增量（正值表示割边增多，负值表示割边减少）
 */
long long EvaluateCutDelta(const CoarsenLevel& level,
                            const FPGAGraph& fg,
                            const PartitionResult& part,
                            int vi, int target_f);

#endif // TOPOPART_INITIAL_H
