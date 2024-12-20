# 生物信息算法

## 一、序列比对算法
RUST-BIO-TYPES

1. 全局比对(Global)
```
序列A: ATCG--AT
序列B: ATCGCCAT
```
- 将两个序列从头到尾完全比对
- 必须考虑序列两端的所有字符
- 适用于长度相近且相似度高的序列
- 例如:比对同源基因或蛋白质序列

2. 半全局比对(Semiglobal)

```
序列A:   TCGAT
序列B: ATCGATT
```
- 允许在序列两端自由引入空位(gap)
- 不对序列两端的空位进行惩罚
- 适用于一个序列是另一个序列的子序列的情况
- 例如:短读序列与参考基因组的比对

3. 局部比对(Local)

```
序列A: AAATCGAAA
      |||
序列B: CCCTCGTTTT
```
- 只寻找序列中相似度最高的片段
- 忽略不相似的区域
- 适用于序列中包含保守区域的情况
- 例如:寻找motif或功能域


## 二、基因组装算法

基因组组装是生物信息学中的一个重要领域，其目的是将短的DNA读段（reads）拼接成更长的连续序列（contigs），最终组装成整个基因组。基因组组装算法主要分为两类：基于重叠图的算法（OLC算法）和基于de Bruijn图的算法。

### 基于重叠图的算法（OLC算法）
OLC（Overlap-Layout-Consensus）算法是基因组组装中的一种经典方法，特别适用于较长读段（reads）的测序数据，如Sanger测序、PacBio和Nanopore数据。OLC算法的基本思想是通过找到读段之间的重叠区域，将它们拼接起来，最终形成完整的基因组序列。以下是OLC算法的详细步骤和原理。
#### 步骤
##### 1. 重叠（Overlap）
重叠步骤的目标是找到所有读段之间的重叠区域。这个过程通常需要大量的计算资源，因为需要将每个读段与所有其他读段进行比对。为了加速这一过程，常用的策略包括：

- 使用索引和哈希表快速查找潜在的重叠区域。
- 仅比对具有相似种子的读段。
##### 2. 布局（Layout）
布局步骤涉及构建一个重叠图，其中每个节点表示一个读段，边表示重叠关系。重叠图的构建和优化可以使用图理论中的多种算法，例如：

- Hamilton路径：找到覆盖所有节点（读段）的路径。
- 最大权匹配：优先选择具有最长重叠区域的边。
- 通过分析和处理重叠图，可以确定读段的最佳排列顺序。

##### 3. 共识（Consensus）
共识步骤的目的是生成一个高质量的基因组序列。这个过程通常包括：

- 多重序列比对：将所有重叠区域的读段进行比对，生成一个多重比对矩阵。
- 共识生成：根据多重比对矩阵，生成共识序列。共识生成算法会考虑每个位置上的碱基频率，选择出现频率最高的碱基作为共识碱基。

#### 优点：
- 适用于较长读长的读段，能够处理复杂的重复序列。
#### 缺点：
- 计算复杂度高，随着读段数量的增加，计算需求呈指数增长。

#### 例子

- 重叠（Overlap）：
找出所有读段之间的重叠区域。每个读段会和其他读段进行比对，识别出它们的重叠部分。常用的工具如BLAST、Minimap2等可以用于此步骤。

```
Read 1: ATCGTAC
Read 2: CGTACGGA
重叠部分：CGTAC
````
- 布局（Layout）：
根据重叠信息构建一个图，每个节点代表一个读段，边表示重叠关系。然后，通过图的分析和处理，将读段按照正确的顺序排列起来。
```
Nodes: [Read 1, Read 2]
Edges: [Read 1 -> Read 2 (重叠部分: CGTAC)]
```
- 共识（Consensus）：
通过遍历图，生成共识序列。共识序列是基于重叠区域的多个读段的序列，以确定最终的基因组序列。这个步骤通常使用多重比对和共识生成算法，如MUSCLE或MAFFT。

#### 算法实现

```
from collections import defaultdict

def find_overlaps(reads, min_overlap=3):
    overlaps = defaultdict(list)
    for i, read1 in enumerate(reads):
        for j, read2 in enumerate(reads):
            if i != j:
                overlap = find_overlap(read1, read2, min_overlap)
                if overlap:
                    overlaps[read1].append((read2, overlap))
    return overlaps

def find_overlap(read1, read2, min_overlap):
    start = max(0, len(read1) - len(read2))
    for i in range(start, len(read1)):
        if read1[i:] == read2[:len(read1) - i]:
            if len(read1) - i >= min_overlap:
                return read1[i:]
    return None

def assemble_genome(overlaps):
    # 简单的贪心算法，从第一个读段开始，逐步拼接
    start_read = list(overlaps.keys())[0]
    assembled = start_read
    current_read = start_read
    while overlaps[current_read]:
        next_read, overlap = overlaps[current_read].pop(0)
        assembled += next_read[len(overlap):]
        current_read = next_read
    return assembled

def main():
    reads = ["ATCGTAC", "CGTACGGA", "GGAATC"]
    overlaps = find_overlaps(reads, min_overlap=3)
    
    assembled_genome = assemble_genome(overlaps)
    print(f"Assembled Genome: {assembled_genome}")

if __name__ == "__main__":
    main()
```

### 基于de Bruijn图的算法
de Bruijn图算法是基因组组装中常用的一种方法，特别适用于处理短读段（reads）数据，如Illumina测序数据。它通过将读段分割成较短的k-mers（长度为k的子序列），并构建一个图来表示这些k-mers之间的关系，从而实现基因组的组装。以下是详细的步骤和原理。


#### 步骤






##### 1. 构建k-mers：
将所有读段分割成长度为k的子序列（k-mers）。每个读段会产生多个k-mers。例如，对一个读段ACGTGCA，如果k=3，则产生的k-mers为ACG，CGT，GTG，TGC，GCA。

##### 2. 构建de Bruijn图：
- 节点：每个k-mer作为一个节点。
- 边：如果两个k-mers在原始读段中相邻，则在它们之间添加有向边。例如，ACG和CGT是相邻的，因此在图中添加一条从ACG指向CGT的边。

##### 3. 简化图结构：
- 移除冗余边和节点：合并重复的边和节点，简化图结构。
- 平滑处理：去除低覆盖度的k-mers和错误的连接，减少噪音。

##### 4. 遍历图：
- 在简化后的de Bruijn图中找到Eulerian路径或回路（遍历所有边且每条边只遍历一次的路径或回路），生成组装的基因组序列。

#### 优点：
- 适用于短读段数据，计算效率高。
- 能够处理高覆盖度和重复序列。
#### 缺点：
- k-mer的选择对组装结果影响较大，选择不当可能导致图结构复杂化。
- 对于长重复序列和低覆盖度区域，处理较为困难。

#### 例子
- 假设我们有以下读段：

```
Read 1: ACGT
Read 2: CGTG
Read 3: GTGC
Read 4: TGCA
```
- 如果k=3，则生成的k-mers为：

```ACG, CGT, GTG, TGC, GCA```

- 构建de Bruijn图：

```
Nodes: ACG, CGT, GTG, TGC, GCA
Edges: ACG -> CGT, CGT -> GTG, GTG -> TGC, TGC -> GCA
```

#### 算法实现
```
from collections import defaultdict, deque

def build_de_bruijn_graph(kmers):
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

def find_eulerian_path(graph):
    def find_cycle(start):
        cycle = []
        stack = [start]
        while stack:
            u = stack[-1]
            if graph[u]:
                v = graph[u].pop()
                stack.append(v)
            else:
                cycle.append(stack.pop())
        return cycle[::-1]

    in_degrees = defaultdict(int)
    out_degrees = defaultdict(int)
    for u in graph:
        out_degrees[u] = len(graph[u])
        for v in graph[u]:
            in_degrees[v] += 1

    start = next((u for u in graph if out_degrees[u] > in_degrees[u]), None)
    start = start or next(iter(graph))

    cycle = find_cycle(start)
    path = deque(cycle)
    while any(graph.values()):
        for i, node in enumerate(path):
            if graph[node]:
                sub_cycle = find_cycle(node)
                path.rotate(-i)
                path.popleft()
                path.extend(sub_cycle)
                path.rotate(i)
                break
    return list(path)

def main():
    reads = ["ACGT", "CGTG", "GTGC", "TGCA"]
    k = 3
    kmers = [read[i:i+k] for read in reads for i in range(len(read) - k + 1)]
    graph = build_de_bruijn_graph(kmers)
    path = find_eulerian_path(graph)
    
    # 输出组装结果
    assembled_genome = path[0]
    for node in path[1:]:
        assembled_genome += node[-1]
    print(f"Assembled Genome: {assembled_genome}")

if __name__ == "__main__":
    main()
```



### 实际应用中的工具
以下是一些常用的基因组组装工具：

- Velvet：基于de Bruijn图的短读段组装工具，适用于Illumina数据。
- SPAdes：改进的de Bruijn图算法，能够处理多种测序平台的数据。
- Canu：基于OLC算法的长读段组装工具，适用于PacBio和Nanopore数据。
- Flye：另一款适用于长读段数据的组装工具，支持PacBio和Nanopore。




