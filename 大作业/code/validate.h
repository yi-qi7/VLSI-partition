#ifndef VALIDATE_H
#define VALIDATE_H

/**
 * @file    validate.h
 * @brief   划分合法性校验 + 输入输出模块
 *
 * 校验接口 CheckPartitionValid：
 *   自动校验拓扑 / 资源 / 固定节点三大约束，输出违规数量
 *
 * 输入输出接口：
 *   支持读取自定义网表、FPGA 拓扑、固定节点配置
 */

#include "types.h"
#include <string>

// ================================================================
//  校验
// ================================================================

/**
 * @brief 校验结果结构体
 */
struct ValidationReport {
    bool valid;                         // 整体是否合法（三大约束全部满足）
    int fixed_violations;              // 固定节点违规数
    int resource_violations;           // 资源超限违规数（超限 FPGA 数）
    int topology_violations;           // 拓扑违规数（不满足直连条件的跨 FPGA 网数）
    string detail;                     // 详细信息

    ValidationReport() : valid(true), fixed_violations(0),
                         resource_violations(0), topology_violations(0) {}
};

/**
 * @brief 校验划分结果合法性（三大约束）
 * @param res  划分结果
 * @param g    电路图（需要固定节点信息）
 * @param fg   FPGA 拓扑图（需要资源容量、直连矩阵）
 * @return 校验报告 ValidationReport
 *
 * 校验项：
 *   约束1 固定节点约束：固定节点 partition[i] == node2fpga[i]
 *   约束2 资源约束：每个 FPGA 承载节点数 ≤ resource_cap[f]
 *   约束3 拓扑约束：每条网涉及的 FPGA 集合构成团（两两直连）
 */
ValidationReport CheckPartitionValid(const PartitionResult& res,
                                      const CircuitGraph& g,
                                      const FPGAGraph& fg);

// ================================================================
//  输入
// ================================================================

/**
 * @brief 读取 FPGA 拓扑文件
 * @param filename  文件路径
 * @param fg        输出 FPGA 图
 *
 * 文件格式（纯文本）：
 *   第1行：fpga_num  edge_num
 *   后续每行：u v（一条无向边，0-indexed）
 *   后续：每行一个 int 表示该 FPGA 资源容量（共 fpga_num 行，可选）
 */
bool ReadFPGATopology(const string& filename, FPGAGraph& fg);

/**
 * @brief 读取电路网表文件（多引脚网格式，内部转为 star 模型的二引脚网）
 * @param filename  文件路径
 * @param g         输出电路图
 *
 * 文件格式（纯文本）：
 *   第1行：node_num  net_num（可分行）
 *   后续每行：一组节点 id（一条多引脚网）
 *   后续：每行一组固定到该 FPGA 的节点 id（可选）
 */
bool ReadCircuitNetlist(const string& filename, CircuitGraph& g);

/**
 * @brief 读取固定节点配置文件
 * @param filename  文件路径
 * @param g         电路图（修改内部 fixed_nodes / node2fpga）
 *
 * 文件格式（纯文本）：
 *   每行：node_id  fpga_id
 */
bool ReadFixedNodes(const string& filename, CircuitGraph& g);

// ================================================================
//  输出
// ================================================================

/**
 * @brief 导出划分结果到文件
 * @param filename  输出文件路径
 * @param res       划分结果
 *
 * 文件格式：
 *   每行一个整数，第 i 行 = 节点 i 分配的 FPGA id
 */
void ExportPartition(const string& filename, const PartitionResult& res);

/**
 * @brief 打印划分结果汇总统计
 * @param res  划分结果
 * @param fg   FPGA 拓扑图（提供资源容量信息）
 */
void PrintPartitionStats(const PartitionResult& res, const FPGAGraph& fg);

#endif // VALIDATE_H
