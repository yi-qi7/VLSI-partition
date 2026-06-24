# TopoPart — ICCAD 2021 多 FPGA 拓扑驱动划分框架

> **多级拓扑驱动多 FPGA 系统划分**
>
> 纯 CPU 实现 · C++17 · MinGW g++ 编译 · 双模式初始划分 · Bitset 加速

---

## 一、项目结构

```
大作业/
├── code/
│   ├── main.cpp                 # CLI 入口，5 步流水线调度，--mode 开关
│   ├── types.h / types.cpp      # 核心数据结构（FPGAGraph, CircuitGraph, CandidateSet, CoarsenLevel, PartitionResult）
│   ├── utils.h / utils.cpp      # Dijkstra 全点对最短路径，按需 BFS（BFSFromSource），FPGA/电路距离预计算
│   ├── candidate.h / candidate.cpp  # 模块2：候选 FPGA 传播（Algorithm 1）
│   ├── coarsen.h / coarsen.cpp      # 模块3：候选集感知多层粗化（自适应收敛 + 带重试）
│   ├── initial.h / initial.cpp      # 模块4：初始划分（paper/grow 双模式）
│   ├── refine.h / refine.cpp        # 模块5：解粗化 + 多级 KL-FM + 拓扑保持负载均衡 + 智能重均衡
│   ├── validate.h / validate.cpp    # 三大约束校验 + 文件输入输出 + 统计输出
│   └── Makefile                     # MinGW g++ 编译
├── FPGA Graph/                      # MFS1（8 FPGA）/ MFS2（43 FPGA，214 边）
├── Generated Benchmarks/            # case1~case8（最大 300K 节点，749K 网）
├── ISPD_benchmark/                  # ibm01~ibm18（.hgr 格式）
├── Titan23 Benchmarks/              # 22 个大规模电路（含 .k8.part 参考解）
├── build/                           # 编译产物（topopart.exe + .o 文件）
├── result/                          # 输出划分文件
└── test.bat                         # 一键编译 + 测试脚本
```

---

## 二、算法流水线

```
输入：电路网表 G(E,V) + FPGA 拓扑图 Ĝ(Ê,V̂) +（可选）固定节点 V_γ
  │
  ├─ 步骤 0：读取输入 + 预计算
  │    · ReadFPGATopology：FPGA 全点对最短距离（Dijkstra）+ Ŝ 缓存（S_cache）
  │    · ReadCircuitNetlist：Star 模型多引脚网转换 + 固定节点解析
  │    · BuildCircuitDist：按需 BFS（O(1) 内存，不存储 N×N 矩阵）
  │
  ├─ 模块2：候选 FPGA 传播 CandidateFPGAPropagation（Algorithm 1）
  │    · BFS 风格传播：固定节点 → 按定理 III.1 收缩可移动节点候选集
  │    · 每个节点 vi 维护位图候选集 Cddt[vi]（uint64_t，支持 64 FPGA）
  │    · 候选集收缩为单元素 → 节点转为新固定节点，继续传播
  │
  ├─ 模块3：候选集感知多层粗化 Coarsening（自适应收敛 + 带重试）
  │    · 双条件合并：①电路边相连 ②候选集交集 |Cddt[u] ∩ Cddt[v]| ≥ β
  │    · 按节点度降序贪心匹配，优先合并高度节点（合并收益更大）
  │    · Bitset 加速求交（__popcnt64 单指令）
  │    · 自适应终止：连续 3 轮压缩率 < 1.02 或超节点数 ≤ min(5, 2K) 或 50 轮上限
  │    · 输出：约 13 层（300K 节点 → ~30K 超节点）
  │
  ├─ 模块4：初始划分 InitialPartition（双模式，--mode 切换）
  │    ┌─ paper 模式（默认）：ICCAD 2021 Algorithm 2 + 回溯 Traceback
  │    │   · 优先队列 Q（候选集小者优先分配）
  │    │   · 贪心分配：优先违规最少 → 负载最轻的 FPGA
  │    │   · 约束传播：分配后收缩邻居候选集（交集 Ŝ(v̂i, 1)）
  │    │   · 冲突回溯：邻居候选集变空 → 撤销分配 + 恢复备份 + 移除失败 FPGA
  │    │   · 回溯超限（默认 100 次）自动回退到 grow 模式
  │    └─ grow 模式（--mode grow）：区域生长 + 拓扑保持扩散
  │        · BFS 区域生长：连通分量统一分配到同一 FPGA（零违规保证）
  │        · 自适应拓扑保持扩散：超载 FPGA → 轻载 FPGA（仅接受零违规迁移）
  │        · 自适应收敛：连续 5 轮无改进 → 停止
  │        · 多种子随机重启（3 次选最优）
  │
  ├─ 模块5：逐层解粗化 UncoarsenRefine
  │    · 继承上层划分 → 拆分超节点（子节点继承超节点 FPGA）
  │    · 每层 RefineSingleLevel：FM 风格边界节点迁移
  │      - 拓扑硬约束：new_violations ≤ old_violations（不增加违规）
  │      - 复合增益：cut_delta - topo_improve×100 - balance_gain×2
  │    · 细层（L0~L3）KL-FM 割边优化：增益排序 + 锁定机制 + 最佳前缀回滚
  │    · L0 层 3 轮交替循环：
  │      ┌─ 阶段 A：拓扑保持负载均衡（超载 FPGA 节点 → 轻载 FPGA，零违规迁移）
  │      ├─ 阶段 B：KL-FM 割边优化（Kerighan-Lin 风格，允许临时恶化，最佳前缀回滚）
  │      ├─ 阶段 C：智能重均衡（容忍 ≤2 割边增加，超载→欠载迁移）
  │      └─ 阶段 D：最终快速 KL-FM 清理（5 pass）
  │
  └─ 步骤 5：三大约束校验 + 输出
       · CheckPartitionValid：固定节点 / 资源容量 / 拓扑距离三大约束
       · ExportPartition：输出每行一个 FPGA id 的划分文件
       · PrintPartitionStats：各 FPGA 负载分布 + 利用率统计
```

---

## 三、核心约束（定理 III.1）

| 约束 | 检查方式 | 代码位置 |
|------|----------|----------|
| 固定节点 | `partition[vi] == required_fpga` | [`validate.cpp:31`](code/validate.cpp:31) |
| 资源容量 | `|fpga_nodes[f]| ≤ resource_cap[f]` | [`validate.cpp:51`](code/validate.cpp:51) |
| 拓扑距离 | 每条边 `(vi,vj)`：`dist_fpga(part[vi], part[vj]) ≤ 1` | [`validate.cpp:69`](code/validate.cpp:69) |

**定理 III.1**：vi → v̂i, vj → v̂j 无拓扑违规 ⟺ `dist_circuit(vi, vj) ≥ dist_fpga(v̂i, v̂j)`

候选传播收缩逻辑（[`candidate.cpp:95-131`](code/candidate.cpp:95)）：
- 已知 vi 固定绑定 v̂i，电路距离 k = dist(vi, vj)
- 收缩 `Cddt[vj] = Cddt[vj] ∩ Ŝ(v̂i, k)`
- 若 Cddt 为空 → 无可行解；若收缩为单元素 → 转为新固定节点

---

## 四、关键优化技术

| 优化 | 说明 | 代码位置 |
|------|------|----------|
| **Bitset 候选集** | `uint64_t` 替代 `unordered_set`，`__popcnt64` 求大小 | [`types.h:84-147`](code/types.h:84) |
| **按需 BFS** | `BFSFromSource` 模板函数替代 N×N 距离矩阵 | [`utils.h:62-85`](code/utils.h:62) |
| **自适应收敛** | 所有循环：连续 N 轮无改进 → 停止（避免无效迭代） | 各模块循环内 |
| **多级 KL-FM** | L0~L3 解粗化时 KL-FM + L0 3 轮交替循环 | [`refine.cpp:234-468`](code/refine.cpp:234) |
| **智能重均衡** | 超载→欠载零违规迁移（容忍 ≤2 割边增加） | [`refine.cpp:405-443`](code/refine.cpp:405) |
| **拓扑硬约束守卫** | 迁移前检查 `new_violations ≤ old_violations`，拒绝违规增加 | [`refine.cpp:108-115`](code/refine.cpp:108) |
| **冲突回溯机制** | Paper 模式：邻居候选集备份/恢复，失败 FPGA 永久移除 | [`initial.cpp:163-187`](code/initial.cpp:163) |
| **多种子重启** | Grow 模式 3 次随机种子，评分选最优（viol×10000 + over×1000 + cut/100） | [`main.cpp:219-263`](code/main.cpp:219) |
| **双模式划分** | paper（Algorithm 2 + 回溯）/ grow（区域生长 + 扩散） | [`initial.cpp:62-70`](code/initial.cpp:62) |
| **S_cache 预缓存** | FPGA Ŝ(v̂i, x) 前缀累积预计算，O(1) 查询 | [`types.cpp:15-42`](code/types.cpp:15) |

---

## 五、编译与运行

### 编译

```powershell
cd code
mingw32-make clean    # 清理构建产物
mingw32-make all      # 编译（生成 ../build/topopart.exe）
```

或使用一键脚本（在 `大作业/` 目录下）：

```powershell
.\test.bat            # 默认测试 case1 + MFS2
.\test.bat case2 MFS2 # 指定 case 和 FPGA 拓扑
```

### 运行

```powershell
# paper 模式（默认，Algorithm 2 + 回溯）
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2"

# grow 模式（区域生长，大图更快更稳定）
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2" --mode grow

# 详细日志
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2" -v

# 手动设置 β 阈值
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2" --beta 28

# 指定输出文件
.\build\topopart.exe "Generated Benchmarks/case1" "FPGA Graph/MFS2" -o result/output.part
```

### 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `<circuit_file>` | 电路网表文件（二引脚网格式） | 必填 |
| `<fpga_topo_file>` | FPGA 拓扑文件 | 必填 |
| `[fixed_file]` | 固定节点配置文件（可选第3位置参数） | 无 |
| `--mode paper\|grow` | 初始划分模式 | paper |
| `-o <file>` | 输出文件路径 | `../result/<bench>.topopart.part` |
| `-v` / `--verbose` | 详细日志 | 关闭 |
| `--beta N` | 手动设置 β 粗化阈值 | `⌊2/3 × K⌋` |

---

## 六、输入格式

**电路网表**（例：case1）：
```
300000        ← 节点数
749091        ← 网数
212496 250542 ← 网：节点 u 节点 v（Star 模型：首节点为中心）
85805 88841
...
```

**FPGA 拓扑**（例：MFS2）：
```
43 214        ← FPGA 数量，边数
0 5           ← 边：FPGA u  FPGA v
0 12
...
```

---

## 七、模块职责速查

| 文件 | 职责 | 对应论文章节 |
|------|------|------------|
| [`main.cpp`](code/main.cpp) | CLI 解析，5 步流水线调度，性能统计 | — |
| [`types.h/.cpp`](code/types.h) | FPGAGraph, CircuitGraph, CandidateSet（bitset）, CoarsenLevel, PartitionResult | III-A/B |
| [`utils.h/.cpp`](code/utils.h) | Dijkstra, BFSFromSource（按需 BFS）, BuildFPGADist | III-B |
| [`candidate.h/.cpp`](code/candidate.cpp) | 候选 FPGA 传播（Algorithm 1） | IV-B |
| [`coarsen.h/.cpp`](code/coarsen.cpp) | 多层粗化 + 自适应终止 + 带重试 | IV-C |
| [`initial.h/.cpp`](code/initial.cpp) | 双模式初始划分（paper Algorithm 2 + grow 区域生长） | IV-D |
| [`refine.h/.cpp`](code/refine.cpp) | 解粗化 + 多级 KL-FM + 拓扑保持负载均衡 + 智能重均衡 | IV-D |
| [`validate.h/.cpp`](code/validate.cpp) | 三大约束校验（CheckPartitionValid），文件 I/O，统计输出 | — |

---

## 八、关键算法参数

| 参数 | 含义 | 默认值 | 代码位置 |
|------|------|--------|----------|
| `beta_ratio` | β 初始占比（β = ratio × K） | 2/3 | [`types.h:187`](code/types.h:187) |
| `max_coarsen_retry` | 粗化失败最大重试次数 | 5 | [`types.h:188`](code/types.h:188) |
| `beta_grow_factor` | 每次重试 β 增大因子 | 1.1 | [`types.h:189`](code/types.h:189) |
| `max_traceback_depth` | Paper 模式最大回溯深度 | 100 | [`types.h:190`](code/types.h:190) |
| `max_iter`（RefineSingleLevel） | 单层 FM 最大迭代轮数 | 20（最粗层）/ 10（中间层） | [`refine.cpp:173-232`](code/refine.cpp:173) |
| L0 交替轮数 | L0 均衡→KL-FM→重均衡循环 | 3 | [`refine.cpp:280`](code/refine.cpp:280) |
| Grow 模式随机重启 | 多种子评分选最优 | 3 | [`main.cpp:216`](code/main.cpp:216) |
| 自适应收敛 stall | 连续无改进轮数阈值 | 5（扩散）/ 3（其他） | 各模块 |

---

## 九、参考文献

1. Murray et al., *"TopoPart: a Multi-level Topology-Driven Partitioning Framework for Multi-FPGA Systems,"* ICCAD 2021.
2. Fiduccia & Mattheyses, *"A Linear-Time Heuristic for Improving Network Partitions,"* DAC 1982.
3. Kernighan & Lin, *"An Efficient Heuristic Procedure for Partitioning Graphs,"* BSTJ 1970.
4. Karypis & Kumar, *"Multilevel k-way Hypergraph Partitioning,"* DAC 1999.

---

*最后更新：2026-06-24 · 2026 VLSI Design 课程大作业*
