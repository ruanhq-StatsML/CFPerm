# Graph localization（不用 community detection / 不用 ego-net）

> 前提不变：**图谱特征**进 \(X\)；归因门仍是 cmean / 缩支撑  
> （[`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md)）。  
> 本文问：在已有边级/点级漂移分之后，能不能用 **别的 graph-localization**  
> 把支撑收成一块 **连通（或近似连通）的图区域**——  
> **不用**全局 community detection，也 **不用** ego-centric 邻域当核。

---

## 0. 问题改写成一句

已有分数（来自协议，非图算法）：

```text
s(v) 或 s(e)   ←  cmean r_v / PO g_v / 边级 shift
G=(V,E)        ←  交互图：user–item 二分、covisit、转化路径图 …
```

求子集 \(C\subset V\)（或边集）使得：

- 分数尖（捕获 \(\gamma\) 高）  
- 在 \(G\) 上有结构约束（连通 / 低导率 / 扩散质量高）  
- **不是**全图社区划分，**不是**「种子的 1-hop ego」

这才叫 graph-localization；community / ego 是我们显式排除的两条。

---

## 1. 排除什么、为什么

| 排除 | 原因 |
|---|---|
| **Community detection**（Louvain 等） | 目标是全局划分模块度，不是「最大化漂移分的支撑」；社区边界≠漂移边界 |
| **Ego-centric net** | 固定半径邻域，与分数无关；高分点的 ego 常又脏又大，不是 localization |

图谱特征（covisit 计数、credit…）照旧用；这里多的是 **可选的第二阶段：在 \(G\) 上收核**。

---

## 2. 可用的 graph-localization 族（调研口径）

### 2.1 Seeded local expansion（有种子，无 ego 硬切）

种子 = 高 \(r_v\) / 高 \(g_v\) 的少量点（来自 cmean 行谱），然后：

| 方法 | 想法 | 输出 |
|---|---|---|
| **Personalized PageRank (PPR)** | 从种子重启的扩散质量 | 质量排序 → sweep 取前缀 |
| **PageRank-Nibble / Nibble** | 截断 PPR + conductance sweep | 低导率局部簇（体积约与簇大小成正比的时间） |
| **Heat kernel diffusion** | 热核代替 PPR | 同上，带宽可调 |

与 ego 的差别：边界由 **扩散质量 / 导率** 决定，不是固定 1-hop。  
与 community 的差别：**局部**、服务种子，不做全图划分。

**接到协议：**  
`D_drill` 给出种子集合 \(S_0=\mathrm{Top}(r)\) → PPR-Nibble 扩成连通核 \(C\) → 用 cmean 的 \(R(C)\) 做 D_check。  
若 \(R\) 不降，撤回（结构扩出来的是假核）。

### 2.2 Score-driven anomalous connected subgraph（可无指定社区）

经典：**GraphScan / subset scan with connectivity**（Neill et al.）：

- 输入：图 + 点上的异常分（我们用漂移分）+ 可加/满足 LTSS 的 score  
- 输出：分数最高的 **连通** 子图（不规则形状）  
- 相对 ULS：可扩展且（在设定下）保证最高分连通子集；相对穷举 FlexScan：可算更大 \(k\)

变体：additive graph scan、prize-collecting Steiner、密度约束子图。

**接到协议：**  
点分 \(s(v)=r_v\) 或边分聚合到点 → GraphScan 得 \(C^\star\) → 当作 \(K\)，再 \(R\) / \(\pi\)。  
**不**先跑 Louvain。

### 2.3 Upper Level Set (ULS) / 分数阈值 + 连通分量

```text
按 s(v) 降序；阈值 τ 以上的点导出诱导子图；取最大/最优连通分量
```

便宜、好实现；不保证全局最优连通高分集，但是强基线。  
和「只取 TopK id、不管边」比：多了 **连通** 约束，仍不是 community。

### 2.4 Graph-smoothed scores（TV / fused lasso / 拉普拉斯平滑）

\[
\min_z\ \sum_v (z_v-s_v)^2 + \lambda\sum_{(u,v)\in E}|z_u-z_v|
\]

平滑后再阈值化 / 取高水平集。  
作用：把零散高分点收成图上连贯区域；\(\lambda\) 控「结构 vs 分数」。  
仍不是社区发现（没有全局 \(k\) 个社区的目标）。

### 2.5 局部谱 / 种子谱定位

种子相关的局部特征向量、或 spectral embedding 后只在种子邻域子图上切。  
介于 Nibble 与异常子图之间；实现重于 ULS，轻于全图谱聚类。

---

## 3. 推荐怎么「弄一下看看」（实验阶梯）

在 **同一套 cmean 分数** 上比四条 localization 后端（结构约束递增）：

| # | 后端 | 结构？ | 实现难度 | 优先 |
|---|---|---|---|---|
| 0 | 现在：Top\(r_v\) / entity key（无图连通） | 无 | 已有 | 基线 |
| 1 | **ULS**：高分点诱导子图的连通分量 | 连通 | 低 | **先做** |
| 2 | **PPR-Nibble**：Top\(r\) 为种子扩局部簇 | 低导率局部 | 中 | 第二 |
| 3 | **GraphScan 风格**：最大分连通子图 | 最优连通（近似/精确） | 中–高 | 有 1–2 后再上 |

评估仍用协议量（非社区指标）：

\[
\alpha=\frac{|C|}{|S|},\quad
\gamma=1-R(C),\quad
\mathrm{eff}=\gamma/\alpha,\quad
\pi\text{-稳定性}
\]

外加可选：\(C\) 的边密度 / 直径（描述用，不当因果）。

图 \(G\) 建议：

- user–item 二分（TencentGR 交互）  
- 或 item–item covisit（已有图谱特征同源）  
- 点 = 当前 `drill_key` 的 value（商户图 / 用户图），边 = 共现或业务邻接  

---

## 4. 和主协议怎么叠（口径）

```text
图谱特征 → X
cmean δ, D, r_v     → 分数 + 是否值得收支撑（行谱）
若 D_drill：
   可选 graph-loc(G, s=r) → C
   D_check: R(C)；不过则退回 Top(r) 或停
FSDS 仍只在最终 K* 上一次
```

Justify 补一句：

> Graph localization enforces a *connectivity* (or conductance) constraint on the retrieval support. It does not estimate communities or treatment effects; seeds/scores come from period conditional-mean contrasts.

禁止话术：「检测到异常社区」「ego 影响范围」。

---

## 5. 还有其他吗

| 可玩 | 默认？ |
|---|---|
| ULS / PPR-Nibble / GraphScan | 建议按阶梯试 |
| Graph-TV 平滑再 ULS | 可选 |
| 多层边（点击+covisit）作 \(G\) 的加权并 | 可选；仍非 community |
| Community detection 当 localization | **否** |
| Ego-1hop 当核 | **否** |
| GNN 端到端归因 | **否**（与「图谱特征」路线冲突） |

---

**收束：**  
能。不用 community、不用 ego，用 **分数驱动的图局部化**：ULS → PPR-Nibble → GraphScan。  
分数仍来自 cmean；图只提供连通/导率约束；验收仍用 \(R,\gamma,\pi\)。  
这和「只用图谱特征、主门是 cmean」完全兼容——多的是可选的结构收核后端。
