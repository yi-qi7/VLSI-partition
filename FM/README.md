# 基于 FM 算法的超图划分

一个高性能的超图划分实现，采用 Fiduccia-Mattheyses (FM) 算法，支持多起点优化和 OpenMP 并行化。

## 目录

- [简介](#简介)
- [特性](#特性)
- [项目结构](#项目结构)
- [环境要求](#环境要求)
- [如何运行](#如何运行)
- [算法详解](#算法详解)
- [参数配置](#参数配置)
- [输出结果](#输出结果)

## 简介

本项目实现了 Fiduccia-Mattheyses (FM) 算法用于超图二划分。目标是将超图的顶点划分为两个均衡的子集，同时最小化被切割的超边数量（cut size）。

FM 算法是一种经典的迭代改进算法，广泛应用于 VLSI 设计自动化领域。本实现包含多项优化措施以获得更好的解质量。

## 特性

- **多起点 FM**：从不同初始划分运行多个独立的 FM 实例
- **双重初始化策略**：结合随机初始化和 BFS 初始化，增加多样性
- **OpenMP 并行化**：并行执行多个 FM 运行实例
- **增量增益更新**：每次移动后高效地重新计算增益
- **早停机制**：检测到无改进时提前终止
- **可视化进度条**：实时显示运行进度
- **批量处理**：一键测试所有数据集

## 项目结构

```
.
├── main.cpp              # 主程序入口
├── solution.cpp          # FM 算法实现
├── solution.h            # Solution 类头文件
├── Graph.cpp             # 图数据结构
├── Graph.h               # Graph 类头文件
├── Node.cpp              # 节点数据结构
├── Node.h                # Node 类头文件
├── Net.cpp               # 超边（Net）数据结构
├── Net.h                 # Net 类头文件
├── evaluate.cpp          # Cut size 评估函数
├── evaluate.h            # 评估函数头文件
├── Makefile              # 编译配置
├── run_benchmarks.sh     # 批量测试脚本
└── ISPD_benchmark/       # 基准测试数据集 (ibm01-ibm18)
```

## 环境要求

- **编译器**：支持 C++17 的 g++
- **OpenMP**：用于并行执行
- **操作系统**：Linux 

## 如何运行

### 步骤 1：编译项目

```bash
make all
```

### 步骤 2：准备数据集

解压基准测试数据集：

```bash
unzip ISPD_benchmark.zip
```

### 步骤 3：运行单个测试

```bash
./main ISPD_benchmark/ibm01.hgr
```

### 步骤 4：批量测试（推荐）

```bash
chmod +x run_benchmarks.sh
./run_benchmarks.sh
```

这将自动测试全部 18 个基准数据集（ibm01 - ibm18）并保存结果。

## 算法详解

### 核心 FM 算法

FM 算法以 pass 为单位进行迭代。每个 pass 包含：

1. **初始化**：计算所有节点的增益值
2. **节点移动**：重复移动增益最高且满足平衡约束的未锁定节点
3. **回滚**：恢复到 pass 过程中遇到的最佳划分
4. **终止**：当没有改进时停止

### 关键优化

#### 1. 多起点策略

```cpp
#define MULTI_START 20
```

不运行单次 FM，而是从不同初始划分运行 20 次，选择最佳结果。这大大降低了陷入较差局部最优的可能性。

#### 2. 双重初始化方法

**随机初始化** (`initialize_partition`)：
- 随机打乱所有节点
- 将前半部分分配到分区 X，后半部分分配到分区 Y

**BFS 初始化** (`initialize_partition_bfs`)：
- 从随机节点开始
- 使用 BFS 遍历超图
- 将访问的节点分配到分区 X 直到平衡
- 这种方法倾向于保持连接的节点在一起

算法交替使用这两种方法以增加多样性。

#### 3. OpenMP 并行化

```cpp
#pragma omp parallel
{
    #pragma omp for schedule(dynamic)
    for(int run = 0; run < MULTI_START; run++) {
        // 每个线程运行独立的 FM 实例
    }
}
```

多个 FM 运行并行执行，有效利用多核 CPU。

#### 4. 增量增益更新

每次移动后不重新计算所有增益，而是使用增量更新：

```cpp
void Solution::update_gains_incremental(Node *moved_node, ...);
```

只有受移动影响的节点才更新增益，基于关键节点概念：
- **源分区关键节点**：源分区中唯一的节点
- **目标分区关键节点**：目标分区中唯一的节点

这将复杂度从 O(|V| × |E|) 降低到 O(被移动节点的度)。

#### 5. 早停机制

```cpp
if(new_cut >= prev_cut && pass > 5) {
    break;
}
```

如果 cut size 在 5 个 pass 后仍无改进，算法提前终止以节省计算时间。

#### 6. 平衡约束

```cpp
#define EPSILON 0.02
```

划分必须保持 2% 容差内的平衡：
- 分区 X 比例：48% ~ 52%
- 分区 Y 比例：48% ~ 52%

## 参数配置

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `EPSILON` | 0.02 | 平衡容差（±2%） |
| `MAX_PASSES` | 50 | 每次 FM 运行的最大 pass 数 |
| `MULTI_START` | 20 | 独立 FM 运行的次数 |

可在 `solution.cpp` 中调整：

```cpp
#define EPSILON 0.02
#define MAX_PASSES 50
#define MULTI_START 20
```

## 输出结果

### 单次运行输出

```
Num nodes: 12752
Num nets: 14111

Progress: [=========================>                        ] 10/20 runs

Partition completed.
Size of partition X: 6376
Size of partition Y: 6376
Ratio X: 0.5
Ratio Y: 0.5
Best cut size found: 55
Partition result saved to: ibm01_partition.txt
Cut size: 55
```

### 批量运行输出文件

**benchmark_results.txt** - 完整结果：
```
Benchmark       Cut Size        Partition Ratio
----------------------------------------
ibm01.hgr       55              X: 0.50    Y: 0.50
ibm02.hgr       120             X: 0.50    Y: 0.50
...
```

**cut_sizes.txt** - 简洁格式，便于分析：
```
ibm01 55
ibm02 120
ibm03 85
...
```

**{benchmark}_partition.txt** - 每个数据集的划分结果：
```
0
1
0
1
...
```
（0 = 分区 X，1 = 分区 Y）

## 基准测试数据集

ISPD98 基准测试套件包含 18 个超图（ibm01 - ibm18），广泛用于评估划分算法。

| 数据集 | 节点数 | 超边数 |
|--------|--------|--------|
| ibm01 | 12,752 | 14,111 |
| ibm02 | 19,642 | 19,568 |
| ... | ... | ... |
| ibm18 | 210,613 | 201,920 |

## 许可证

本项目仅供教育目的使用。

## 参考文献

- Fiduccia, C. M., & Mattheyses, R. M. (1982). A linear-time heuristic for improving network partitions. *DAC '82*.
