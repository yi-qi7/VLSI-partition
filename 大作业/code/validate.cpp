/**
 * @file    validate.cpp
 * @brief   校验、文件I/O与统计功能实现
 */

#include "validate.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

// ================================================================
//  校验：三大约束检查
// ================================================================

ValidationReport CheckPartitionValid(const PartitionResult& res,
                                      const CircuitGraph& g,
                                      const FPGAGraph& fg) {
    ValidationReport report;
    int N = (int)res.node2fpga.size();
    int K = fg.fpga_num;

    // ----------------------------------------------------------
    // 约束1：固定节点约束
    //   每个固定节点 partition[i] == node2fpga[i]
    // ----------------------------------------------------------
    for (int vi : g.fixed_nodes) {
        if (vi >= N) continue;  // 防御越界
        auto it = g.node2fpga.find(vi);
        if (it == g.node2fpga.end()) continue;

        int required_f = it->second;
        int actual_f = res.node2fpga[vi];

        if (actual_f != required_f) {
            ++report.fixed_violations;
            report.detail += "Fixed violation: node " + to_string(vi) +
                             " should be FPGA " + to_string(required_f) +
                             ", actual FPGA " + to_string(actual_f) + "\n";
        }
    }

    // ----------------------------------------------------------
    // 约束2：资源约束
    //   每个 FPGA 承载节点数 ≤ resource_cap[f]
    // ----------------------------------------------------------
    for (int f = 0; f < K; ++f) {
        int load = (int)res.fpga2nodes[f].size();
        int cap = fg.resource_cap[f];
        if (load > cap) {
            ++report.resource_violations;
            report.detail += "Resource overflow: FPGA " + to_string(f) +
                             " load=" + to_string(load) + ","
                             "capacity=" + to_string(cap) + "\n";
        }
    }

    // ----------------------------------------------------------
    // 约束3：拓扑约束（定理 III.1）
    //   对每条边 (vi, vj)，验证：
    //     若两个 FPGA 不同，则它们之间必须有直接边（dist=1）
    //   注意：这里检查的是 star 模型下的二引脚网
    // ----------------------------------------------------------
    // 拓扑约束验证：逐边检查，每条无向边仅计数一次（vi < vj）
    if (!fg.dist_all.empty()) {
        for (int vi = 0; vi < N; ++vi) {
            int f_vi = res.node2fpga[vi];
            if (f_vi < 0) continue;

            for (int vj : g.adj[vi]) {
                if (vj >= N) continue;
                if (vi >= vj) continue;   // 每条无向边仅计数一次
                int f_vj = res.node2fpga[vj];
                if (f_vj < 0) continue;

                if (f_vi == f_vj) continue;  // 同 FPGA，OK

                // 不同 FPGA：需要满足直连条件（dist_fpga == 1）
                int fpga_dist = fg.dist_all[f_vi][f_vj];
                if (fpga_dist > 1) {
                    ++report.topology_violations;
                    // 仅报告前几条违规详情
                    if (report.topology_violations <= 5) {
                        report.detail += "Topology violation: edge (" + to_string(vi) + "," +
                                         to_string(vj) + ") across FPGA " +
                                         to_string(f_vi) + "<->" + to_string(f_vj) +
                                         ", FPGA distance = " + to_string(fpga_dist) + "\n";
                    }
                }
            }
        }
    }

    // 汇总判定
    report.valid = (report.fixed_violations == 0 &&
                    report.resource_violations == 0 &&
                    report.topology_violations == 0);

    return report;
}

// ================================================================
//  输入：读取 FPGA 拓扑文件
// ================================================================

bool ReadFPGATopology(const string& filename, FPGAGraph& fg) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "[Error] Cannot open FPGA topology file: " << filename << endl;
        return false;
    }

    int K, M;
    file >> K >> M;
    file.ignore();

    fg.fpga_num = K;
    fg.adj.assign(K, vector<int>());
    fg.resource_cap.assign(K, 0);  // 默认 0 = 未设置（由调用者设置）

    // 读取 FPGA 间连线
    for (int i = 0; i < M; ++i) {
        int u, v;
        file >> u >> v;
        fg.adj[u].push_back(v);
        fg.adj[v].push_back(u);
    }

    // 尝试读取资源容量（可选字段）
    string line;
    getline(file, line);  // 跳过换行
    for (int i = 0; i < K; ++i) {
        if (getline(file, line)) {
            istringstream iss(line);
            int cap;
            if (iss >> cap) {
                fg.resource_cap[i] = cap;
            }
        }
    }

    file.close();

    // 自动计算距离
    BuildFPGADist(fg);

    cout << "[ReadFPGATopology] Loaded: " << K << " FPGAs, "
         << M << " edges" << endl;
    return true;
}

// ================================================================
//  输入：读取电路网表文件
// ================================================================

bool ReadCircuitNetlist(const string& filename, CircuitGraph& g) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "[Error] Cannot open netlist file: " << filename << endl;
        return false;
    }

    string line;

    // 跳过注释行，读取第一行数据
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        break;
    }

    // 解析第一行：node_num（可能与 net_num 同行或分行）
    istringstream iss(line);
    int node_cnt = 0, net_cnt = 0;
    iss >> node_cnt;
    // 尝试从同一行读取 net_cnt
    if (!(iss >> net_cnt)) {
        // net_cnt 在下一行（如 case1 格式：节点数和网数分行）
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            istringstream iss2(line);
            if (iss2 >> net_cnt) break;
        }
    }

    if (node_cnt <= 0 || net_cnt <= 0) {
        cerr << "[Error] Invalid netlist header: nodes=" << node_cnt
             << " nets=" << net_cnt << endl;
        return false;
    }

    g.node_num = node_cnt;
    g.adj.assign(node_cnt, vector<int>());

    // 读取每条网
    // Bug fix: 使用 for 循环代替 while(getline && condition)
    //   原 while 的 getline 在 net_loaded==net_cnt 时会多读一行，
    //   导致第一个 FPGA 的固定节点行被吞掉
    int net_loaded = 0;
    for (int ni = 0; ni < net_cnt; ++ni) {
        // 跳过空行和注释
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            break;
        }
        istringstream net_iss(line);
        vector<int> net_nodes;
        int node;
        while (net_iss >> node) {
            if (node >= 0 && node < node_cnt) {
                net_nodes.push_back(node);
            }
        }
        if (net_nodes.size() >= 2) {
            int center = net_nodes[0];
            for (size_t i = 1; i < net_nodes.size(); ++i) {
                int other = net_nodes[i];
                g.adj[center].push_back(other);
                g.adj[other].push_back(center);
            }
        }
        ++net_loaded;
    }

    // ---- 读取固定节点（第四部分） ----
    // 剩余行：每行对应一个 FPGA，行中节点固定分配到该 FPGA
    // 第 i 行剩余行 → FPGA i（从 0 开始编号）
    int fixed_cnt = 0;
    int fpga_idx = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        istringstream fix_iss(line);
        int node;
        while (fix_iss >> node) {
            if (node >= 0 && node < node_cnt) {
                g.fixed_nodes.insert(node);
                g.node2fpga[node] = fpga_idx;
                ++fixed_cnt;
            }
        }
        ++fpga_idx;
    }

    file.close();

    // 初始化可移动节点
    g.init_move_nodes();

    cout << "[ReadCircuitNetlist] Loaded: " << node_cnt << " nodes, "
         << net_loaded << " nets";
    if (fixed_cnt > 0)
        cout << ", " << fixed_cnt << " fixed nodes (" << fpga_idx << " FPGAs)";
    cout << endl;
    return true;
}

// ================================================================
//  输入：读取固定节点配置文件
// ================================================================

bool ReadFixedNodes(const string& filename, CircuitGraph& g) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "[Warning] Fixed node file not found: " << filename
             << " -- no fixed constraints." << endl;
        return false;
    }

    string line;
    int cnt = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        int node_id, fpga_id;
        if (iss >> node_id >> fpga_id) {
            g.fixed_nodes.insert(node_id);
            g.node2fpga[node_id] = fpga_id;
            ++cnt;
        }
    }

    file.close();
    g.init_move_nodes();

    cout << "[ReadFixedNodes] Loaded: " << cnt << " fixed nodes" << endl;
    return true;
}

// ================================================================
//  输出：导出划分结果
// ================================================================

void ExportPartition(const string& filename, const PartitionResult& res) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "[Error] Cannot write partition file: " << filename << endl;
        return;
    }

    // 每行一个整数，第 i 行 = 节点 i 分配的 FPGA id
    for (int f : res.node2fpga) {
        file << f << "\n";
    }

    file.close();
    cout << "[ExportPartition] Partition exported: " << filename
         << " (" << res.node2fpga.size() << " lines)" << endl;
}

// ================================================================
//  输出：打印划分结果统计
// ================================================================

void PrintPartitionStats(const PartitionResult& res, const FPGAGraph& fg) {
    int K = fg.fpga_num;

    cout << "\n========== Partition Stats ==========" << endl;
    cout << "  Total nodes: " << res.node2fpga.size() << endl;
    cout << "  Cut size: " << res.cut_size << endl;
    cout << "  Topology violations: " << res.violation_cnt << endl;
    cout << "  FPGA load distribution:" << endl;

    int total_load = 0;
    for (int f = 0; f < K; ++f) {
        int load = (int)res.fpga2nodes[f].size();
        total_load += load;
        double util = (fg.resource_cap[f] > 0)
                      ? (100.0 * load / fg.resource_cap[f]) : 0.0;
        cout << "    FPGA " << setw(3) << f << ": "
             << setw(6) << load << " nodes";
        if (fg.resource_cap[f] > 0) {
            cout << " (cap " << fg.resource_cap[f]
                 << ", util " << fixed << setprecision(1) << util << "%)";
        }
        cout << endl;
    }
    cout << "  Total load: " << total_load << " nodes" << endl;
    cout << "==================================\n" << endl;
}
