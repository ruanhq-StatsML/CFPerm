# 多层图算法速览 × 多模态归因逻辑

> 基础方法论：[`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md)、
> [`Feature_Entity_Joint_CMean.tex`](./Feature_Entity_Joint_CMean.tex)。  
> 本文：**(1)** 多层图算法调研（怎么用、怎么不用）；**(2)** 多模态归因如何挂在同一套
> conditional-mean / 缩支撑协议上。

---

## A. 多层图算法赶紧调研（面向我们的取舍）

### A.1 对象是什么

| 术语 | 含义 | 和我们的对应 |
|---|---|---|
| **Multiplex** | 同一批节点，多层边（不同关系） | 同 user/item，多关系：点击 / covisit / 转化路径 |
| **Multi-layer / multi-aspect** | 层可不同节点集或不同语义 | 商户层、用户层、订单层；或 video/audio/text 块 |
| **Temporal multi-slice** | 时间片当层 | W1 / W2（我们用 period，不当 community 层） |
| **Heterogeneous graph** | 多类型节点+边 | user–item–merchant 异构，常被做成**图谱特征**而非跑 GNN |

经典综述口径（Kim et al. SIGMOD Rec.; Huang et al. DMKD; Magnani et al. CSUR; 近期 multilayer network science reviews）：

社区发现三条策略：

1. **Flattening** — 层合成单图再聚类（丢层语义）  
2. **Layer-by-layer + consensus** — 每层社区再集成  
3. **Direct multilayer** — 多层模块度 / Infomap / SBM / 张量分解 / 耦合随机游走  

表示学习：矩阵分解、多层随机游走嵌入、多层 GNN / 对比学习（层内消息 + 层间融合）。

### A.2 算法族速查（我们关心的性质）

| 族 | 代表想法 | 输出 | 算力 | 对归因的用途 |
|---|---|---|---|---|
| 多层模块度 / GenLouvain | 层内密 + 层间耦合 | 社区划分 | 中–高 | 可选：**造 entity key**（社区 id），不是归因门本身 |
| Infomap / 随机游走 | 多层流向压缩 | 社区 / flow | 中 | 同上；或边权特征 |
| 张量 / NMF | 邻接栈分解 | 因子、嵌入 | 高 | 嵌入进 \(X\) 的一列块 |
| 多层 GNN | 层内聚合+层间融合 | 节点嵌入 | 高 | 嵌入当模态块；**不替代** cmean drill |
| Flatten + 单层算法 | 简单 | 粗社区 | 低 | 基线；易混层 |

**关键结论（对我们）：**  
多层图算法擅长 **结构划分 / 嵌入**。  
我们的归因门擅长 **两窗分布差的支撑收缩**。  

**Scope 更新（更硬）：**  
- 默认 **不用** community detection，**不用** ego-centric network 分析（超 scope）。  
- **不用** ULS / PPR-Nibble / GraphScan 等图上收核（已标超纲 / deferred）。  
- 图谱特征进 \(X\) + cmean 门 **已足够闭环**。  
- 社区/GNN/图 localization 若以后单开，不进默认路径。

**不要**用 Louvain / GNN / 图收核直接当「要不要下钻」的判决器——归因门只认 \(\delta,R\)。

### A.3 和「图谱特征 only、不做图算法」的关系

TencentGR 路径已锁定：用 covisit / credit / activity 等 **图谱导出特征**，不跑社区发现或 GNN 归因。  
这与调研不矛盾：

- 调研告诉我们多层结构**可以**提供 key 或特征块；  
- 默认协议 **不依赖**它们也能闭环（表格特征 + merchant/user/order key 即可）；  
- 若上多层算法，只作为 **上游特征/键生成**，下游仍走 §B。

### A.4 时间多层

把 W1/W2 当 multi-slice 层去做社区，和我们把 \(W\) 当 **period 标签** 比均值，是两条线：

| | 时间片社区 | 我们的 period cmean |
|---|---|---|
| 问 | 结构社区是否跨时间稳定 | 特征均值差落在哪块 mass |
| 输出 | 社区标签 | \(K^\star\) + \(J^\star\) + FSDS |
| 因果 | 不谈 | 不谈 |

可并列：社区标签当 key；**Drill 仍看该 key 下的 \(r_v\)**。

---

## B. 多模态归因逻辑（同一套方法论）

「模态」在协议里不是特殊公民，只是 **\(X\) 的列块** 或 **并列 entity key**。

### B.1 两种多模态接法

**接法 M1 — 特征块（MSR-VTT 式）**

\[
X=\bigl[X^{(\mathrm{video})}\,\|\,X^{(\mathrm{audio})}\,\|\,X^{(\mathrm{text})}\bigr]
\]

- 列 guidance：\(\delta\) 的 share **按块求和** → 模态能量份额（可对齐 RF-VIMP / PO-LOGO，但是 cmean 版更便宜）  
- 行 drill：entity = video_id / user_id / … 照旧  
- 联合 F→E：\(J^\star\) 可先限制在高能量模态块内，再算 \(r_v\)

**接法 M2 — 关系层 / 图多层当 key 或子块**

- 每层边统计 → 特征子块 \(X^{(\ell)}\)  
- 或每层社区 id → 并列 `entity_key`  
- 用 \(\mathrm{eff}^{(\ell)}=\gamma/\alpha\) 比哪一层更好收支撑

TencentGR 默认更像 **M1（图谱特征拼进 \(X\)）+ 嵌套 key（merchant→user→order）**；  
MSR-VTT 脚本是 **M1 + 模态 LOGO 份额**（PO/RF/MMD），尚未接行谱早停——可接到同一决策表。

### B.2 闭环（多模态版，与单模态同一公式）

```text
1. 标准化（W1 scaler）→ X（可含多模态块 / 图谱特征）
2. 提名（可选 MMD TopK / 业务父集）→ 父支撑 S
3. 同一套 cmean：
     δ, share     → 列：哪模态/哪维；本层 FSDS？
     D, r_v       → 行：沿 drill_key 下不下钻
     F→E / E→F    → 联合顺序（须声明）
4. 多 key / 多模态块 → 比 eff、π；分栏或交集
5. 早停得 K* → 一次 FSDS
6. Justify：非因果、非 subgroup（J1–J8）
```

### B.3 和仓库里已有多模态板的对齐

| 已有 | 角色 | 接到本协议 |
|---|---|---|
| `run_msrvtt_mm_attribution.py` | 模态份额 RF / PO-LOGO / MMD-LOGO | = **列/块 guidance** 的重型版；可加 cmean 块份额作轻量门 |
| TencentGR standardize→MMD→FSDS | 图谱特征 + item/merchant 定位 | = 主路径；补 **行谱早停**（勿强制三层钻穿） |
| CF parsimonious smoke | 并列透镜 | **不进** Drill 乘积 |
| 多层图算法（若引入） | 社区/嵌入上游 | 只喂 \(X\) 或 key，不喂门 |

### B.4 Justify（多模态时更容易写歪）

多模态报告常写成「video 效应最大」。替换口径：

- 「video **块**上的 period mean-shift / 份额最大」  
- 「在 key=`video_id` 下行谱尖 → 漂移载体是少数样本」  
- LOGO \(\Delta\) = 去掉该块后风险变化，是 **贡献/份额**，不是 ATE  

多层图社区同理：社区是检索键，不是处理臂。

---

## C. 还有其他的吗（短清单）

| 项 | 默认？ | 说明 |
|---|---|---|
| 图谱特征进 \(X\) | 是 | 已落地 |
| cmean 列/行 + 早停 | 是 | 方法论闭环 |
| 模态块份额（cmean 或 LOGO） | 建议 | 多模态列 guidance |
| 并列 entity key / 多层社区 key | 可选 | 比 eff |
| F→E 联合 | 建议 | 降噪 |
| 真多层 GNN / GenLouvain 归因门 | **否** | 目标函数不对 |
| CF τ 方差当异质 | **否** | 已否决 |
| Seeded PPR / ULS / GraphScan / conductance | **否（本阶段超纲）** | deferred；见 [`GRAPH_LOCALIZATION_NO_COMMUNITY.md`](./GRAPH_LOCALIZATION_NO_COMMUNITY.md) |
| Ego-centric / 全局社区发现 | **否** | 超 scope |

---

## D. 参考文献（调研入口）

- Kim & Lee, *Community Detection in Multi-Layer Graphs: A Survey*, SIGMOD Record.  
- Huang et al., *A survey of community detection methods in multilayer networks*, DMKD.  
- Magnani et al., *Community Detection in Multiplex Networks*, ACM CSUR.  
- Recent reviews on multilayer network science / multilayer graph embedding（社区、动力学、表示学习）。  

公式原型：[`Feature_Entity_Joint_CMean.tex`](./Feature_Entity_Joint_CMean.tex)。

---

**收束：**  
多层图算法 = 可选的 **结构上游**（特征块 / entity key）。  
多模态归因 = 同一套 **conditional-mean 缩支撑**：模态是列块，实体是行键，联合是 F→E/E→F，停止看行谱与业务顶。  
不把图社区或 GNN 直接当成下钻判决——那才和这套方法论对齐。
