# 轻量 Graph Localization（不用 community / ego-network）

> **Scope 钉死：**  
> - ✅ 图谱特征进 \(X\)（已落地）  
> - ✅ 可选：**轻量 graph localization** 辅助缩支撑 / 验收核  
> - ❌ Community detection（Louvain / Infomap / 多层模块度 / SBM …）  
> - ❌ Ego-centric network analysis（角色、ego-motif、结构性洞 …）— 明显超 scope  
> - ❌ 训练 GNN 当归因器  
>
> 基础门仍是 [`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md) 的 **cmean \(\delta/D/r/R\)**。  
> 本文只谈：在已有边（user–item、covisit）上，还能不能用**别的 localization** 方法帮一把。

---

## 1. 「Graph localization」在这里指什么

不是「发现社区」，是：

> 给定种子或分数，在图上找出 **漂移质量集中的局部支撑**（一小撮节点/边），
> 计算量跟种子集走，不跟全图社区划分走。

经典 CS 里接近：**seeded local partitioning**（ACL 风格 PPR）、热核扩散、奖赏子图 / 剥皮——  
都是 **local**，不是 global community detection，也不是 ego 社会学分析。

和我们的 \(\alpha,\gamma\) 同构：缩支撑、看捕获率；只是多了一条 **图邻接约束**。

---

## 2. 能用的几条（按推荐顺序）

图 \(G\)：现成即可——TencentGR 上优先 **user–item 二部图** 或 **item–item covisit**（特征工程里已有 covisit，不必新造世界观）。

### 方法 L1 — Seeded diffusion（主候选）

1. 用 cmean 行谱得种子 \(Z=\mathrm{Top}(r_v)\)（或高 \(g_v\) 实体）  
2. 在 \(G\) 上做 **Personalized PageRank / 热核**，teleport/种子 = \(Z\)  
3. 取 sweep 前缀质量最好的集合 \(K_{\mathrm{diff}}\)（或分数阈值）  
4. 用同一套 \(R=\|\delta_{S\setminus K}\|^2/\|\delta_S\|^2\) 验收  

| | |
|---|---|
| 要什么 | 种子在图上「连着的」漂移载体，补纯 Top(\(r\)) 的碎点 |
| 不要什么 | 全图社区；ego 展开到 \(k\)-hop 角色特征 |
| 算力 | 近似 PPR / push 与 \(|Z|\)、边数相关，可控 |
| Justify | 仍是 period 漂移支撑的检索；扩散是几何约束，不是效应 |

伪码：

```text
Z = Top(r) ∩ mass_ok ∩ π        # cmean 种子
p = ApproximatePPR(G, seeds=Z, alpha)
K = sweep(p) that minimizes conductance or maximizes γ(cmean)
accept if R(K) ≤ R* and mass ok
```

### 方法 L2 — Kernel conductance 验收（更轻，建议常开）

不重新找核，只验收 cmean 提出的 \(K\)：

\[
\phi(K)=\frac{e(K,S\setminus K)}{\min(\mathrm{vol}(K),\mathrm{vol}(S\setminus K))}
\]

| | |
|---|---|
| \(\phi\) 低 | 核在图上相对「成块」→ 与业务叙事更合 |
| \(\phi\) 高 | 核是分数尖但图上碎 → 可标 fragile，或改用 L1 扩散收紧 |

**不当主 Drill 门**；当 D_check 的图侧旁证（与 \(R\) 并列）。

### 方法 L3 — Prize-guided peel（图约束剥皮）

节点奖赏 \(p_v=r_v\)（或 \(r_v\cdot n_v\)），在 \(G[S]\) 上反复去掉 **奖赏/度数比差** 的点（或只按奖赏剥，边仅用于禁止孤立爆炸）：

- 与文档里「剥离曲线」同族，但可要求剩余集诱导子图连通 / 低 conductance  
- 仍读 \(\gamma\) 累积；不是社区目标

### 方法 L4 — 1-hop 邻域只进特征（其实还是图谱特征）

\[
X_i \leftarrow X_i \,\|\, \mathrm{mean}\{X_j:j\sim i\}
\]

W1/W2 仍走 cmean。这是 **特征扩维**，不是 localization 算法——可做，但别叫 graph localization。

---

## 3. 明确不用的（超 scope）

| 方法 | 为何刷掉 |
|---|---|
| Louvain / Leiden / Infomap / 多层模块度 / SBM | 全局社区；目标≠漂移捕获 |
| Ego-network 角色、结构洞、ego-motif 谱 | 分析重、叙事社会学化，离归因门远 |
| 全图 betweenness / 昂贵中心性 | 成本与收益都不对 |
| 训练 GNN 做归因 | 另一套学习问题 |
| 「先社区再在社区内 cmean」当默认 | 偷偷把 community 变主路径 |

若业务以后强需求「类目/圈子」——用 **已有类目字段当 entity key**，或离线一次性社区 id 当 key；  
**默认协议不内嵌社区发现。**

---

## 4. 怎么嵌进现有闭环（仍非因果）

```text
图谱特征 → X
cmean: δ（列）/ r（行）→ 候选核 K_cmean
可选:
  L2  φ(K_cmean) 验收
  L1  以 Top(r) 为种子扩散 → K_diff，用 R 比 γ
  L3  奖赏剥皮曲线
取 γ/α/π 更好者（分栏报告，勿 silent 并）
停：行平 | 业务顶 | R-reject
FSDS 一次在 K*
```

话术：

- ✅ 「在 covisit/UI 图上，以 cmean 种子做局部扩散得到的支撑」  
- ❌ 「检测到的社区里效应更大」  
- ❌ 「ego 网络显示该用户处于结构洞」

---

## 5. 和「只用图谱特征」是否冲突？

**不冲突。**  

- 默认：**只用图谱特征 + cmean** 已闭环。  
- L1–L3 是 **可选插件**：多一次图上的局部运算，仍服务同一 \(\delta/R\) 验收。  
- 若图边质量差 / 覆盖低：关掉 L1–L3，退回纯 cmean——协议仍成立。

建议实验顺序（若试）：

1. 基线：纯 cmean Top(\(r\)) + \(R\)  
2. +L2 conductance 报告  
3. +L1 seeded PPR，比 \(\gamma,\phi,\pi\)  
4. 不上社区、不上 ego  

---

## 6. 还有其他的吗？

同级、仍在 scope 内：

- 二部图上只沿 user 侧或 item 侧扩散（半边 localization）  
- 扩散限制在父核 \(S\) 诱导子图（与 drill 父集一致）  
- holdout 窗上复现 \(K\) 的 overlap  

仍超 scope：多层社区共识、ego 角色、表示学习 GNN。

---

**收束：**  
能。不用 community detection，不用 ego-centric。  
用 **cmean 种子 + 轻量局部扩散 / conductance 验收 / 奖赏剥皮** 做 graph localization，验收仍看 \(\delta,R\)。  
默认可关；开了也只是缩支撑的图约束插件，不是第二套归因哲学。
