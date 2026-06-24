#ifndef INITIAL_H
#define INITIAL_H

/**
 * @file    initial.h
 * @brief   模块4：拓扑驱动贪心初始划分（算法2，带回溯机制）
 *
 * 论文 Algorithm 2: Topology-Driven Greedy Initial Partitioning with Traceback
 *
 * 输入：最粗层超图层级、FPGA 资源拓扑图
 * 输出：粗粒度零违规初始划分结果
 *
 * 内部流程：
 *   1. 初始化：固定超节点直接绑定指定 FPGA；可移动节点标记 node2fpga = -1
 *   2. 优先队列：
 *      Q：priority_queue<(候选集大小, 节点 id), greater> —— 候选越少越优先分配
 *   3. 主分配循环：
 *      a. 弹出 Q 顶部节点，尝试在候选集中寻找最优 FPGA
 *      b. 遍历所有未分配邻节点，收缩其候选集以保持拓扑合法性
 *      c. 冲突检测：若任意邻节点候选集为空 → 触发回溯
 *   4. 回溯恢复逻辑：
 *      a. 撤销当前节点分配，将失败 FPGA 永久从候选集中删除
 *      b. 恢复所有邻节点的候选集备份
 *      c. 当前节点重新入队，选择其他合法 FPGA
 *   5. 分配完成后校验：资源不超限、拓扑违规数为零，返回划分结果
 */

#include "types.h"
#include <queue>
#include <stack>
#include <string>

/**
 * @brief 初始划分（双模式）
 * @param coarsest_level  最粗层超图层级
 * @param fg              FPGA 拓扑图
 * @param max_traceback   最大回溯深度（默认 100）
 * @param mode            "paper" = 算法2+回溯 | "grow" = 区域生长+扩散
 * @return 粗粒度初始划分结果
 */
PartitionResult InitialPartition(CoarsenLevel& coarsest_level,
                                  const FPGAGraph& fg,
                                  int max_traceback = 100,
                                  const std::string& mode = "paper");

/**
 * @brief 割边增量评估：计算将节点 vi 分配到 FPGA target_f 后的割边变化
 * @param level    当前粗化层级
 * @param part     当前划分状态
 * @param vi       待分配节点 id
 * @param target_f 目标 FPGA id
 * @return 割边增量（正值表示割边增多，负值表示割边减少）
 */
long long EvaluateCutDelta(const CoarsenLevel& level,
                            const PartitionResult& part,
                            int vi, int target_f);

#endif // INITIAL_H
