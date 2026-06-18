#ifndef TOPOPART_REFINE_H
#define TOPOPART_REFINE_H

/**
 * @file    topopart_refine.h
 * @brief   TopoPart 模块5：逐层解粗化与约束感知细化 UncoarsenRefine
 *
 * 论文 §IV-D: Uncoarsening and Constraint-Aware Refinement
 *
 * 输入：全部粗化层级列表、FPGA 拓扑图
 * 输出：原始细粒度电路最优合法划分
 *
 * 分层迭代逻辑（从最粗层向原始层逐层还原）：
 *   1. 当前层继承上一层划分结果作为初始解；
 *      拆分超节点 → 子节点继承对应 FPGA 分配
 *   2. 筛选边界节点：存在跨 FPGA 连线的节点，
 *      仅对边界节点做迁移优化
 *   3. 迭代细化循环：
 *      a. 遍历每个边界节点 vi，仅遍历 Cddt[vi] 内合法 FPGA
 *      b. 计算 vi 迁移至每个候选 FPGA 后的割边变化；
 *         选择割边下降最大的 FPGA 执行迁移
 *      c. 迁移全程校验资源上限、拓扑约束，禁止违规移动
 *   4. 单层终止：一轮遍历无割边下降 → 进入下一层解粗化
 *   5. 全部层级还原完成 → 返回原始电路划分结果
 */

#include "topopart_types.h"

/**
 * @brief 逐层解粗化与约束感知细化
 * @param levels  全部粗化层级列表（levels[0] = 原始层，levels[L-1] = 最粗层）
 * @param fg      FPGA 拓扑图
 * @return 原始细粒度电路最优合法划分结果
 *
 * 约束保证：
 *   - 资源约束：每次移动前检查目标 FPGA 容量
 *   - 拓扑约束：每次移动前检查拓扑违规（定理 III.1 验证）
 *   - 固定节点约束：固定节点不允许移动
 */
PartitionResult UncoarsenRefine(vector<CoarsenLevel>& levels, const FPGAGraph& fg);

/**
 * @brief 单层细化迭代
 * @param level     当前层粗化层级
 * @param fg        FPGA 拓扑图
 * @param max_iter  最大迭代轮数（一轮 = 遍历所有边界节点一次）
 * @return 本轮是否有割边下降（true 表示有改进）
 *
 * 内部步骤：
 *   1. 标记边界节点（存在跨 FPGA 邻接的节点）
 *   2. 对每个边界节点，遍历其 Cddt 内的候选 FPGA
 *   3. 选择割边下降最大且合法的 FPGA 迁移
 *   4. 若割边下降 > 0 则执行迁移
 */
bool RefineSingleLevel(CoarsenLevel& level, const FPGAGraph& fg, int max_iter = 10);

#endif // TOPOPART_REFINE_H
