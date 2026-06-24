/**
 * @file    main.cpp
 * @brief   多 FPGA 拓扑驱动划分框架主入口（ICCAD 2021）
 *
 * 算法流水线（固定顺序，不可调换）：
 *   步骤0：读取输入 + 预计算全点对距离
 *   步骤1：候选 FPGA 传播（算法1）
 *   步骤2：候选集感知多层粗化
 *   步骤3：拓扑驱动贪心初始划分（算法2 + 回溯）
 *   步骤4：逐层解粗化 + 约束感知细化（FM）
 *   步骤5：三大约束校验 + 输出结果
 *
 * 用法：
 *   main.exe <电路网表> <FPGA拓扑> [固定节点文件] [-o 输出] [-v]
 *   --mode paper|grow   初始划分模式（默认 paper）
 *   --beta N            手动设置 β 阈值（默认 2/3 × FPGA数）
 */

#include "types.h"
#include "utils.h"
#include "candidate.h"
#include "coarsen.h"
#include "initial.h"
#include "refine.h"
#include "validate.h"

#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <climits>

using namespace std;
using namespace std::chrono;

// ================================================================
//  打印命令行使用说明
// ================================================================

void PrintUsage(const char* prog) {
    cout << "\n===== Multi-level Topology-Driven Partitioning for Multi-FPGA Systems (ICCAD 2021) =====\n\n"
         << "Usage:\n"
         << "  " << prog << " <circuit_file> <fpga_topo_file> [options]\n\n"
         << "Arguments:\n"
         << "  circuit_file    Circuit netlist file (2-pin net format)\n"
         << "  fpga_topo_file  FPGA topology file\n"
         << "  fixed_file      Fixed node config file (optional, 3rd positional arg)\n\n"
         << "Options:\n"
         << "  -o <file>       Output partition file (default: ../result/<bench>.part)\n"
         << "  -v              Verbose logging mode\n"
         << "  --mode <M>      Init mode: paper (Algorithm 2+traceback) | grow (region-growing)\n"
         << "  --beta <N>      Manually set beta threshold (default: 2/3 * fpga_num)\n"
         << "  --help          Show this help\n"
         << endl;
}

// ================================================================
//  主函数入口
// ================================================================

int main(int argc, char* argv[]) {
    // 设置控制台为 UTF-8 编码（防止制表符乱码）
    system("chcp 65001 > nul");

    // 解析命令行参数
    if (argc < 3) {
        PrintUsage(argv[0]);
        return -1;
    }

    string circuit_file;
    string fpga_topo_file;
    string fixed_file;
    string output_file;
    bool verbose = false;
    int manual_beta = -1;
    string init_mode = "paper";  // "paper" or "grow"

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "--beta" && i + 1 < argc) {
            manual_beta = atoi(argv[++i]);
        } else if (arg == "--mode" && i + 1 < argc) {
            init_mode = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (circuit_file.empty()) {
            circuit_file = arg;
        } else if (fpga_topo_file.empty()) {
            fpga_topo_file = arg;
        } else if (fixed_file.empty() && arg[0] != '-') {
            fixed_file = arg;
        }
    }

    if (circuit_file.empty() || fpga_topo_file.empty()) {
        cerr << "[Error] Missing required arguments!\n" << endl;
        PrintUsage(argv[0]);
        return -1;
    }

    // 构造默认输出路径（相对于可执行文件目录的 ../result/）
    if (output_file.empty()) {
        string base = circuit_file;
        size_t pos = base.find_last_of("/\\");
        if (pos != string::npos) base = base.substr(pos + 1);
        pos = base.find_last_of('.');
        if (pos != string::npos) base = base.substr(0, pos);
        // 获取可执行文件所在目录
        string exe_dir = argv[0];
        pos = exe_dir.find_last_of("/\\");
        if (pos != string::npos) exe_dir = exe_dir.substr(0, pos);
        output_file = exe_dir + "/../result/" + base + ".part";
    }

    cout << "\n===== Multi-FPGA Topology-Driven Partitioning (ICCAD 2021) =====\n" << endl;
    cout << "Circuit netlist: " << circuit_file << endl;
    cout << "FPGA topology: " << fpga_topo_file << endl;
    if (!fixed_file.empty())
        cout << "Fixed nodes: " << fixed_file << endl;
    cout << "Output file: " << output_file << endl;
    cout << endl;

    // ============================================================
    //  步骤0：读取输入数据
    // ============================================================
    cout << "========== Step 0: Read Input Data ==========" << endl;

    FPGAGraph fg;
    CircuitGraph g;
    TopoPartConfig cfg;
    cfg.verbose = verbose;

    // 读取 FPGA 拓扑文件
    if (!ReadFPGATopology(fpga_topo_file, fg)) {
        cerr << "[Fatal] Failed to read FPGA topology!" << endl;
        return -1;
    }

    // 读取电路网表文件（含固定节点）
    if (!ReadCircuitNetlist(circuit_file, g)) {
        cerr << "[Fatal] Failed to read circuit netlist!" << endl;
        return -1;
    }

    // 设置默认资源容量（每个 FPGA 容量 = ⌈总节点数 / FPGA数 × 1.25⌉）
    bool need_default_cap = fg.resource_cap.empty();
    if (!need_default_cap) {
        need_default_cap = true;
        for (int cap : fg.resource_cap) {
            if (cap > 1) { need_default_cap = false; break; }
        }
    }
    if (need_default_cap) {
        int default_cap = (int)std::ceil((double)g.node_num / fg.fpga_num * 1.25);
        fg.resource_cap.assign(fg.fpga_num, default_cap);
        cout << "[Info] Default resource capacity: per FPGA " << default_cap << " nodes" << endl;
    }

    // 读取额外的固定节点配置文件（可选）
    if (!fixed_file.empty()) {
        ReadFixedNodes(fixed_file, g);
    }

    // 预计算电路图全点对距离（候选传播阶段必需）
    BuildCircuitDist(g);

    // ============================================================
    //  步骤1：候选 FPGA 传播（算法1）
    // ============================================================
    auto t1_start = high_resolution_clock::now();

    vector<CandidateSet> cddt;
    try {
        cddt = CandidateFPGAPropagation(g, fg);
    } catch (const exception& e) {
        cerr << "\n[Fatal] Candidate propagation failed: " << e.what() << endl;
        return -1;
    }

    auto t1_end = high_resolution_clock::now();
    auto t1_ms = duration_cast<milliseconds>(t1_end - t1_start).count();

    // ============================================================
    //  步骤2：候选集感知多层粗化
    // ============================================================
    auto t2_start = high_resolution_clock::now();

    // ============================================================
    //  步骤2+3：论文 β-重试循环
    //   若初始划分存在违规 → 增大 β 重新粗化
    // ============================================================
    int beta = (manual_beta > 0) ? manual_beta : cfg.get_beta(fg.fpga_num);

    vector<CoarsenLevel> levels;
    PartitionResult initial_result;
    bool zero_violation_found = false;

    for (int retry = 0; retry < cfg.max_coarsen_retry && !zero_violation_found; ++retry) {
        if (retry > 0) {
            beta = max(beta + 1, (int)(beta * cfg.beta_grow_factor));
            // β 不能超过 FPGA 总数，否则无法合并任何节点
            if (beta > fg.fpga_num) beta = fg.fpga_num;
            cout << "\n  [Retry " << retry << "] Increasing beta to " << beta << "..." << endl;
        }
        if (verbose) cout << "  beta threshold = " << beta << endl;

        try {
            levels = Coarsening(g, cddt, fg, beta);
        } catch (const exception& e) {
            cerr << "\n[Fatal] Coarsening failed: " << e.what() << endl;
            return -1;
        }

        PartitionResult best_trial_result;
        long long best_score = LLONG_MAX;
        int num_trials = (init_mode == "grow") ? 3 : 1;

        for (int trial = 0; trial < num_trials; ++trial) {
            vector<CoarsenLevel> trial_levels = levels;
            srand(trial * 12345 + 67890 + retry * 100000);

            PartitionResult trial_result;
            try {
                trial_result = InitialPartition(trial_levels.back(), fg,
                                                 cfg.max_traceback_depth, init_mode);
            } catch (const exception& e) {
                cout << "  Exception: " << e.what() << endl;
                if (init_mode != "grow") {
                    trial_result = InitialPartition(trial_levels.back(), fg,
                                                     cfg.max_traceback_depth, "grow");
                } else continue;
            }
            trial_levels.back().level_part = trial_result;

            int viol = 0, over = 0; long long cut = 0;
            auto& part = trial_result;
            auto& gr = trial_levels.back().super_graph;
            for (int vi = 0; vi < gr.node_num; ++vi) {
                for (int vj : gr.adj[vi]) {
                    int fi = part.node2fpga[vi], fj = part.node2fpga[vj];
                    if (fi >= 0 && fj >= 0 && fi != fj) ++cut;
                    if (fi >= 0 && fj >= 0 && fg.dist_all[fi][fj] > 1) ++viol;
                }
            }
            cut /= 2; viol /= 2;
            for (int f = 0; f < fg.fpga_num; ++f) {
                int ld = 0;
                for (int nid : part.fpga2nodes[f])
                    ld += (nid < (int)trial_levels.back().super_weight.size())
                          ? trial_levels.back().super_weight[nid] : 1;
                if (ld > fg.resource_cap[f]) ++over;
            }
            long long score = (long long)viol * 100000 + over * 1000 + cut / 100;
            if (verbose || num_trials > 1)
                cout << "  Trial " << trial << ": viol=" << viol << " over=" << over
                     << " cut=" << cut << " score=" << score << endl;

            if (score < best_score) {
                best_score = score;
                best_trial_result = std::move(trial_result);
                levels = std::move(trial_levels);
            }
        }

        initial_result = std::move(best_trial_result);
        levels.back().level_part = initial_result;

        if (initial_result.violation_cnt == 0) {
            cout << "  ZERO violations with beta=" << beta << "!" << endl;
            zero_violation_found = true;
        } else {
            cout << "  " << initial_result.violation_cnt << " violations with beta=" << beta << endl;
        }
    }

    auto t2_end = high_resolution_clock::now();
    auto t23_ms = duration_cast<milliseconds>(t2_end - t2_start).count();

    // ============================================================
    //  步骤4：逐层解粗化与约束感知细化
    // ============================================================
    auto t4_start = high_resolution_clock::now();

    PartitionResult final_result;
    try {
        final_result = UncoarsenRefine(levels, fg);
    } catch (const exception& e) {
        cerr << "\n[Fatal] Uncoarsen/refine failed: " << e.what() << endl;
        return -1;
    }

    auto t4_end = high_resolution_clock::now();
    auto t4_ms = duration_cast<milliseconds>(t4_end - t4_start).count();

    // ============================================================
    //  步骤5：三大约束校验 + 输出结果
    // ============================================================
    cout << "\n========== Step 5: Validation & Output ==========" << endl;

    // 校验三大约束（固定节点 / 资源容量 / 拓扑距离）
    ValidationReport report = CheckPartitionValid(final_result, g, fg);

    cout << "\n----- Validation Results -----" << endl;
    cout << "  Fixed node constraint: " << report.fixed_violations
         << " violations  " << (report.fixed_violations == 0 ? "PASS" : "FAIL") << endl;
    cout << "  Resource constraint:  " << report.resource_violations
         << " violations  " << (report.resource_violations == 0 ? "PASS" : "FAIL") << endl;
    cout << "  Topology constraint:  " << report.topology_violations
         << " violations  " << (report.topology_violations == 0 ? "PASS" : "FAIL") << endl;
    cout << "  Overall validity:     "
         << (report.valid ? "PASS" : "FAIL") << endl;
    cout << "------------------------------\n";

    if (!report.detail.empty() && verbose) {
        cout << "\nDetails:\n" << report.detail << endl;
    }

    // 输出各 FPGA 负载分布统计
    PrintPartitionStats(final_result, fg);

    // 导出划分结果到文件
    ExportPartition(output_file, final_result);

    // ============================================================
    //  输出各阶段耗时统计
    // ============================================================
    cout << "\n----- Performance Statistics -----" << endl;
    cout << "  Step1 Candidate Prop: " << t1_ms << " ms" << endl;
    cout << "  Step2+3 Coarsen+Init: " << t23_ms << " ms" << endl;
    cout << "  Step4 UncoarsenRefine:" << t4_ms << " ms" << endl;
    auto total_ms = t1_ms + t23_ms + t4_ms;
    cout << "  Total time:           " << total_ms << " ms" << endl;
    cout << "------------------------------\n" << endl;

    return report.valid ? 0 : 1;
}
