#ifndef REFINE_H
#define REFINE_H

/**
 * @file    refine.h
 * @brief   模块5：解粗化与约束感知细化（FM）
 *
 * 论文 IV-D 节：Uncoarsening and Constraint-Aware Refinement
 * 从最粗层向最细层迭代，在每层执行 FM 风格的边界节点迁移。
 *
 * 强制约束：
 *   - 资源约束：每次迁移前检查 FPGA 容量
 *   - 拓扑约束：每次迁移前按定理 III.1 检查
 *   - 固定节点：固定节点永不移动
 */

#include "types.h"

/** 从最粗层到最细层的解粗化与细化 */
PartitionResult UncoarsenRefine(vector<CoarsenLevel>& levels, const FPGAGraph& fg);

/** 单层 FM 风格细化迭代 */
bool RefineSingleLevel(CoarsenLevel& level, const FPGAGraph& fg, int max_iter = 10);

#endif // REFINE_H
