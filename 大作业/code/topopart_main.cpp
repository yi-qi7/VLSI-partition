/**
 * @file    topopart_main.cpp
 * @brief   TopoPart 主函数测试框架
 *
 * 完整算法流水线（固定顺序，不可调换）：
 *   输入：电路网图 G + FPGA 硬件拓扑图 Ĝ + 固定节点映射表
 *   → 步骤0：读取输入 + 预计算全点对距离
 *   → 步骤1：候选 FPGA 传播 CandidateFPGAPropagation（Algorithm 1）
 *   → 步骤2：候选集感知多层粗化 Coarsening（带 β 阈值自适应回退）
 *   → 步骤3：拓扑驱动贪心初始划分 InitialPartition（Algorithm 2 + Traceback）
 *   → 步骤4：逐层解粗化 + 约束感知细化 UncoarsenRefine
 *   → 步骤5：校验三大约束 + 输出结果
 *   输出：合法划分方案 P
 *
 * 用法（命令行）：
 *   topopart.exe <circuit_file> <fpga_topo_file> [fixed_file] [-o output_file] [-v]
 *
 * 参数：
 *   circuit_file    : 电路网表文件（二引脚网格式）
 *   fpga_topo_file  : FPGA 拓扑文件
 *   fixed_file      : 固定节点配置文件（可选）
 *   -o output_file  : 输出划分文件路径（默认 ../result/xxx.part）
 *   -v              : 详细日志模式
 *   --cuda          : 启用 GPU 加速
 *   --beta N        : 手动设置 β 阈值（默认 2/3 * fpga_num）
 */

#include "topopart_types.h"
#include "topopart_utils.h"
#include "topopart_candidate.h"
#include "topopart_coarsen.h"
#include "topopart_initial.h"
#include "topopart_refine.h"
#include "topopart_validate.h"

#ifdef TOPOPART_USE_CUDA
#include "topopart_cuda.cuh"
#endif

#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <climits>

using namespace std;
using namespace std::chrono;

// ================================================================
//  打印使用说明
// ================================================================

void PrintUsage(const char* prog) {
    cout << "\n╔══════════════════════════════════════════════════════════════╗\n"
         << "║  TopoPart: Multi-level Topology-Driven Partitioning         ║\n"
         << "║  for Multi-FPGA Systems (ICCAD 2021)                       ║\n"
         << "╚══════════════════════════════════════════════════════════════╝\n\n"
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
         << "  --cuda          Enable GPU parallel acceleration (requires CUDA)\n"
         << "  --beta <N>      Manually set beta threshold (default: 2/3 * fpga_num)\n"
         << "  --help          Show this help\n"
         << endl;
}

// ================================================================
//  主函数
// ================================================================

int main(int argc, char* argv[]) {
    // ----------------------------------------------------------
    // 设置控制台为 UTF-8 编码（防止 box-drawing 字符乱码）
    // ----------------------------------------------------------
    system("chcp 65001 > nul");

    // ----------------------------------------------------------
    // 解析命令行参数
    // ----------------------------------------------------------
    if (argc < 3) {
        PrintUsage(argv[0]);
        return -1;
    }

    string circuit_file;
    string fpga_topo_file;
    string fixed_file;
    string output_file;
    bool verbose = false;
    bool use_cuda = false;
    int manual_beta = -1;
    string init_mode = "paper";  // "paper" or "grow"

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "--cuda") {
            use_cuda = true;
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

    // 默认输出路径 (相对于 exe 所在目录的 ../result/)
    if (output_file.empty()) {
        string base = circuit_file;
        size_t pos = base.find_last_of("/\\");
        if (pos != string::npos) base = base.substr(pos + 1);
        pos = base.find_last_of('.');
        if (pos != string::npos) base = base.substr(0, pos);
        // 获取 exe 所在目录
        string exe_dir = argv[0];
        pos = exe_dir.find_last_of("/\\");
        if (pos != string::npos) exe_dir = exe_dir.substr(0, pos);
        output_file = exe_dir + "/../result/" + base + ".topopart.part";
    }

    cout << "\n╔══════════════════════════════════════════════════════════════╗\n"
         << "║  TopoPart — ICCAD 2021 Multi-FPGA Topology-Driven Partitioning ║\n"
         << "╚══════════════════════════════════════════════════════════════╝\n"
         << endl;
    cout << "Circuit netlist: " << circuit_file << endl;
    cout << "FPGA topology: " << fpga_topo_file << endl;
    if (!fixed_file.empty())
        cout << "Fixed nodes: " << fixed_file << endl;
    cout << "Output file: " << output_file << endl;
    if (use_cuda) cout << "GPU acceleration: Enabled" << endl;
    cout << endl;

    // ============================================================
    //  步骤0：读取输入数据
    // ============================================================
    cout << "========== Step 0: Read Input Data ==========" << endl;

    FPGAGraph fg;
    CircuitGraph g;
    TopoPartConfig cfg;
    cfg.verbose = verbose;

    // 读取 FPGA 拓扑
    if (!ReadFPGATopology(fpga_topo_file, fg)) {
        cerr << "[Fatal] Failed to read FPGA topology!" << endl;
        return -1;
    }

    // 设置资源容量默认值（若未在文件中指定）
    if (fg.resource_cap.empty() || fg.resource_cap[0] <= 0) {
        // 默认：平均分配 + 10% 余量
        int avg_load = 0;
        // 先读取电路以获取节点数
    }

    // 读取电路网表
    if (!ReadCircuitNetlist(circuit_file, g)) {
        cerr << "[Fatal] Failed to read circuit netlist!" << endl;
        return -1;
    }

    // 设置默认资源容量（每个 FPGA 容量 = ceil(总节点数 / FPGA数 * 1.25)）
    bool need_default_cap = fg.resource_cap.empty();
    if (!need_default_cap) {
        need_default_cap = true;
        for (int cap : fg.resource_cap) {
            if (cap > 1) { need_default_cap = false; break; }
        }
    }
    if (need_default_cap) {
        int default_cap = (int)ceil((double)g.node_num / fg.fpga_num * 1.25);
        fg.resource_cap.assign(fg.fpga_num, default_cap);
        cout << "[Info] Default resource capacity: per FPGA " << default_cap << " nodes" << endl;
    }

    // 读取固定节点配置
    if (!fixed_file.empty()) {
        ReadFixedNodes(fixed_file, g);
    }

    // 预计算电路图全点对距离（候选传播必需）
    BuildCircuitDist(g);

    // ============================================================
    //  步骤1：候选 FPGA 传播 CandidateFPGAPropagation
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
    //  步骤2：候选集感知多层粗化 Coarsening
    // ============================================================
    auto t2_start = high_resolution_clock::now();

    int beta = (manual_beta > 0) ? manual_beta : cfg.get_beta(fg.fpga_num);
    if (verbose) cout << "  beta threshold = " << beta << endl;

    vector<CoarsenLevel> levels;
    try {
        levels = Coarsening(g, cddt, fg, beta);
    } catch (const exception& e) {
        cerr << "\n[Fatal] Coarsening failed: " << e.what() << endl;
        return -1;
    }

    auto t2_end = high_resolution_clock::now();
    auto t2_ms = duration_cast<milliseconds>(t2_end - t2_start).count();

    // ============================================================
    //  步骤3：初始划分（paper 模式: Algorithm 2+traceback / grow 模式: 区域生长）
    // ============================================================
    auto t3_start = high_resolution_clock::now();

    PartitionResult initial_result;
    int num_trials = (init_mode == "grow") ? 3 : 1;  // grow mode benefits from multi-seed

    long long best_score = LLONG_MAX;
    for (int trial = 0; trial < num_trials; ++trial) {
        vector<CoarsenLevel> trial_levels = levels;
        srand(trial * 12345 + 67890);

        PartitionResult trial_result;
        try {
            trial_result = InitialPartition(trial_levels.back(), fg,
                                             cfg.max_traceback_depth, init_mode);
        } catch (const exception& e) {
            cout << "  Paper mode failed: " << e.what() << endl;
            cout << "  Falling back to grow mode..." << endl;
            trial_result = InitialPartition(trial_levels.back(), fg,
                                             cfg.max_traceback_depth, "grow");
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
        long long score = (long long)viol * 10000 + over * 1000 + cut / 100;
        if (verbose || num_trials > 1)
            cout << "  Trial " << trial << ": viol=" << viol << " over=" << over
                 << " cut=" << cut << " score=" << score << endl;

        if (score < best_score) {
            best_score = score;
            initial_result = std::move(trial_result);
            levels = std::move(trial_levels);
        }
    }
    levels.back().level_part = initial_result;

    auto t3_end = high_resolution_clock::now();
    auto t3_ms = duration_cast<milliseconds>(t3_end - t3_start).count();

    // ============================================================
    //  步骤4：逐层解粗化与细化 UncoarsenRefine
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
    //  步骤5：校验 + 输出
    // ============================================================
    cout << "\n========== Step 5: Validation & Output ==========" << endl;

    // 校验三大约束
    ValidationReport report = CheckPartitionValid(final_result, g, fg);

    cout << "\n╔══════════════════════════════════════════════════════════════╗\n";
    cout << "║                Validation Results                            ║\n";
    cout << "╠══════════════════════════════════════════════════════════════╣\n";
    cout << "║  Fixed node constraint: " << setw(5) << report.fixed_violations
         << " violations  " << (report.fixed_violations == 0 ? "✓" : "✗")
         << "                    ║\n";
    cout << "║  Resource constraint:  " << setw(5) << report.resource_violations
         << " violations  " << (report.resource_violations == 0 ? "✓" : "✗")
         << "                    ║\n";
    cout << "║  Topology constraint:  " << setw(5) << report.topology_violations
         << " violations  " << (report.topology_violations == 0 ? "✓" : "✗")
         << "                    ║\n";
    cout << "╠══════════════════════════════════════════════════════════════╣\n";
    cout << "║  Overall validity:    "
         << (report.valid ? "✓ PASS" : "✗ FAIL")
         << "                                  ║\n";
    cout << "╚══════════════════════════════════════════════════════════════╝\n";

    if (!report.detail.empty() && verbose) {
        cout << "\nDetails:\n" << report.detail << endl;
    }

    // 输出统计信息
    PrintPartitionStats(final_result, fg);

    // 导出划分结果
    ExportPartition(output_file, final_result);

    // ============================================================
    //  性能总结
    // ============================================================
    cout << "\n╔══════════════════════════════════════════════════════════════╗\n";
    cout << "║                Performance Statistics                        ║\n";
    cout << "╠══════════════════════════════════════════════════════════════╣\n";
    cout << "║  Step1 Candidate Prop: " << setw(8) << t1_ms << " ms"
         << "                          ║\n";
    cout << "║  Step2 Multi-Coarsen:  " << setw(8) << t2_ms << " ms"
         << "                          ║\n";
    cout << "║  Step3 Initial Part:   " << setw(8) << t3_ms << " ms"
         << "                          ║\n";
    cout << "║  Step4 UncoarsenRefine:" << setw(8) << t4_ms << " ms"
         << "                          ║\n";
    auto total_ms = t1_ms + t2_ms + t3_ms + t4_ms;
    cout << "╠══════════════════════════════════════════════════════════════╣\n";
    cout << "║  Total time:           " << setw(8) << total_ms << " ms"
         << "                        ║\n";
    cout << "╚══════════════════════════════════════════════════════════════╝\n"
         << endl;

    return report.valid ? 0 : 1;
}
