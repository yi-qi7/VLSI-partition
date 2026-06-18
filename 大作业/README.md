# TopoPart — ICCAD 2021 多 FPGA 拓扑驱动划分框架

> **多级拓扑驱动多 FPGA 系统划分**
>
> 纯 CPU 实现 · C++17 · MSVC 编译 · 双模式初始划分 · Bitset 加速

---

## 一、项目结构

```
大作业/
├── code/
│   ├── topopart_main.cpp            # CLI 入口、5 步流水线、--mode 开关
│   ├── topopart_types.h / .cpp      # 核心数据结构 (bitset CandidateSet)
│   ├── topopart_utils.h / .cpp      # Dijkstra、按需 BFS、FPGA 距离预计算
│   ├── topopart_candidate.h / .cpp  # 模块二：候选 FPGA 传播
│   ├── topopart_coarsen.h / .cpp    # 模块三：候选集感知多层粗化 (自适应收敛)
│   ├── topopart_initial.h / .cpp    # 模块四：初始划分 (paper/grow 双模式)
│   ├── topopart_refine.h / .cpp     # 模块五：解粗化 + KL-FM + 智能重均衡
│   ├── topopart_validate.h / .cpp   # 三大约束校验 + 文件输入输出
│   ├── topopart_cuda.cu / .cuh      # GPU 加速 (可选)
│   └── Makefile                     # MSVC cl.exe 编译
├── FPGA Graph/                      # MFS1(8 FPGA) / MFS2(43 FPGA, 214 边)
├── Generated Benchmarks/            # case1~case8 (最大 300K 节点, 749K 网)
├── ISPD_benchmark/                  # ibm01~ibm18 (.hgr 格式)
├── Titan23 Benchmarks/              # 22 个大规模电路
├── build/                           # topopart.exe + .obj
└── result/                          # 输出划分文件
```

---

## 二、算法流水线

```
输入: 电路网表 G + FPGA 拓扑图 Ĝ + (可选) 固定节点
  │
  ├─ 步骤 0：读取输入 + 预计算
  │    · ReadFPGATopology：FPGA 全点对最短距离 (Dijkstra) + Ŝ 缓存
  │    · ReadCircuitNetlist：电路邻接表构建
  │    · BuildCircuitDist：按需 BFS (O(1) 内存, 不存储 N×N 矩阵)
  │
  ├─ 模块二：候选 FPGA 传播 CandidateFPGAPropagation
  │    · BFS 风格传播：固定节点 → 按定理 III.1 收缩可移动节点候选集
  │
  ├─ 模块三：候选集感知多层粗化 Coarsening (自适应收敛)
  │    · 双条件合并：(1) 电路边相连 (2) 候选集交集 ≥ β
  │    · Bitset 加速求交 (__popcnt64 单指令)
  │    · 自适应终止：连续 3 轮压缩率 < 1.02 或图节点过少
  │    · 输出：约 13 层 (300K 节点 → 30K 超节点)
  │
  ├─ 模块四：初始划分 InitialPartition (双模式, --mode 切换)
  │    ┌─ paper 模式 (默认)：ICCAD 2021 Algorithm 2 + 回溯 Traceback
  │    │   · 优先队列 Q (候选集小者优先)
  │    │   · 贪心分配 + 约束传播
  │    │   · 冲突回溯 (批量备份/恢复邻居候选集)
  │    │   · 回溯超限自动回退到 grow 模式
  │    └─ grow 模式 (--mode grow)：区域生长 + 拓扑保持扩散
  │        · BFS 区域生长：相连节点同 FPGA (零违规保证)
  │        · 自适应拓扑保持扩散 (收敛即停)
  │        · 多种子随机重启 (3 次选最优)
  │
  ├─ 模块五：逐层解粗化 UncoarsenRefine
  │    · 逐层继承上层划分 + RefineSingleLevel
  │    · 细层 (L0~L3) KL-FM 割边优化
  │    · L0：3 轮交替 (均衡 → KL-FM → 智能重均衡)
  │    · 智能重均衡：超载→欠载迁移 (容忍 ≤2 割边增加)
  │
  └─ 步骤 5：三大约束校验 + 输出
       · 固定节点约束 / 资源容量约束 / 拓扑距离约束
```

---

## 三、核心约束 (定理 III.1)

| 约束 | 检查方式 |
|------|----------|
| 固定节点 | `partition[vi] == required_fpga` |
| 资源容量 | `|fpga_nodes[f]| ≤ resource_cap[f]`（默认 `⌈N/K × 1.25⌉`） |
| 拓扑距离 | 每条边 `(vi,vj)`：`dist_fpga(part[vi], part[vj]) ≤ 1` |

**定理 III.1**：vi → v̂i, vj → v̂j 无拓扑违规 ⟺ `dist_circuit(vi, vj) ≥ dist_fpga(v̂i, v̂j)`

---

## 四、关键优化

| 优化 | 说明 | 效果 |
|------|------|------|
| **Bitset 候选集** | `uint64_t` 替代 `unordered_set`，`__popcnt64` 求大小 | 粗化 4.3× 加速 |
| **自适应收敛** | 所有循环：连续 5 轮无改进 → 停止 | 避免无效迭代 |
| **按需 BFS** | `BFSFromSource` 替代 N×N 距离矩阵 | 内存 ~360GB → ~2MB |
| **多级 KL-FM** | L0~L3 解粗化时 KL-FM + L0 交替循环 | 割边优化 |
| **智能重均衡** | 超载→欠载零违规迁移 (容忍 ≤2 割边增加) | 均衡 + 割边权衡 |
| **多种子重启** | grow 模式 3 次随机种子, 评分选最优 | 解质量 |
| **双模式划分** | paper (Algorithm 2+回溯) / grow (区域生长+扩散) | 灵活适配 |

---

## 五、编译与运行

```powershell
cd code
mingw32-make              # CPU 版本编译
mingw32-make clean        # 清理构建产物
```

```powershell
# paper 模式 (默认, Algorithm 2 + 回溯)
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2"

# grow 模式 (区域生长, 大图更快)
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2" --mode grow

# 详细日志
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2" -v
```

### 命令行参数

| 参数 | 说明 |
|------|------|
| `<circuit_file>` | 电路网表文件 (二引脚网格式) |
| `<fpga_topo_file>` | FPGA 拓扑文件 |
| `--mode paper\|grow` | 初始划分模式 (默认 paper) |
| `-o <file>` | 输出文件路径 |
| `-v` | 详细日志 |
| `--beta N` | 手动设置 β 粗化阈值 |
| `--cuda` | 启用 GPU 加速 |

---

## 六、输入格式

**电路网表** (例：case1)：
```
300000        ← 节点数
749091        ← 网数
212496 250542 ← 网：节点 u 节点 v
85805 88841
...
```

**FPGA 拓扑** (例：MFS2)：
```
43 214        ← FPGA 数量, 边数
0 5           ← 边：FPGA u FPGA v
0 12
...
```

---

## 七、性能数据 (case1 + MFS2, 300K 节点, 43 FPGA)

| 指标 | 数值 |
|------|------|
| 拓扑违规 | **0** ✓ |
| 资源违规 | 5 (拓扑约束固有限制) |
| 割边数 | ~375K / 749K 边 (~50%) |
| 总耗时 | ~23s (grow 模式) |

---

## 八、模块职责速查

| 文件 | 职责 |
|------|------|
| `topopart_main.cpp` | CLI 解析, 5 步流水线调度, 性能统计 |
| `topopart_types.h/.cpp` | FPGAGraph, CircuitGraph, CandidateSet (bitset), CoarsenLevel |
| `topopart_utils.h/.cpp` | Dijkstra 全点对最短路径, BFSFromSource, 距离预计算 |
| `topopart_candidate.cpp` | 候选 FPGA 传播 (Algorithm 1) |
| `topopart_coarsen.cpp` | 多层粗化 + 自适应终止 |
| `topopart_initial.cpp` | 双模式初始划分 (paper/grow) |
| `topopart_refine.cpp` | 解粗化 + 多级 KL-FM + 智能重均衡 |
| `topopart_validate.cpp` | 三大约束校验, 文件 I/O, 统计输出 |

---

## 九、参考文献

1. Murray et al., *"TopoPart: a Multi-level Topology-Driven Partitioning Framework for Multi-FPGA Systems,"* ICCAD 2021.
2. Fiduccia & Mattheyses, *"A Linear-Time Heuristic for Improving Network Partitions,"* DAC 1982.
3. Kernighan & Lin, *"An Efficient Heuristic Procedure for Partitioning Graphs,"* BSTJ 1970.
4. Karypis & Kumar, *"Multilevel k-way Hypergraph Partitioning,"* DAC 1999.

---

*最后更新：2026-06-19 · 2026 VLSI Design 课程大作业*
