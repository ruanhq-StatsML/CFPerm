# Leiden 接到多层次下钻：elaborated proposal（审阅稿）

> 目的：把 **Leiden 社区发现** 讲清楚——为什么可能比 dissolve rollup 更合适、  
> 在我们锁定的 **cmean / MMD / FSDS** 协议里它能当什么、不能当什么、怎么落地试。  
> 讲武德：`W`=period；社区是 **检索键 / 结构块**，**不是** ATE / 效应子群。

相关：[`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md) ·  
[`MULTILAYER_GRAPH_AND_MM_ATTR.md`](./MULTILAYER_GRAPH_AND_MM_ATTR.md) ·  
[`GRAPH_LOCALIZATION_NO_COMMUNITY.md`](./GRAPH_LOCALIZATION_NO_COMMUNITY.md)（原「社区超纲」）·  
[`DRILL_DISSOLVE_LANDABILITY.md`](./DRILL_DISSOLVE_LANDABILITY.md)（dissolve 试探结论）

---

## 0. 一句话立场

| | dissolve rollup | **Leiden（本提案）** | 官方 drill（现状） |
|---|---|---|---|
| 在干什么 | 把叶子分 **加总** 到已知父 key | 从边上 **学** 出多分辨率社区键 | 在 **给定** key（商户/用户/订单）上算 δ/r/MMD |
| 优点 | 极便宜 | 键本身带结构；可跨固定 taxonomy | 口径干净、已落地 |
| 缺点 | 用户层对不上 MMD；父 key 必须先验存在 | 曾标超纲；图构造敏感；易被读成「发现因果群」 | 固定层级可能切错颗粒 |
| 建议角色 | 商户 prior（可选） | **上游 entity_key 生成** | **归因门仍只认 δ/r/R + stop** |

**主张：** Leiden 比 dissolve 更适合「多层次下钻」——因为它直接产出 **可下钻的层级键**；  
但 **判决是否下钻 / 是否停止** 仍必须走 cmean 行谱与业务顶，**不**用模块度当 Drill 门。

---

## 1. 为什么觉得 Leiden 比 dissolve 更好（对你关心的点）

### 1.1 Dissolve 卡在哪

试探结论（PR #73）：

- 商户层 `∑ shift_l2` 与 top-down MMD **有重叠**（Jaccard@5≈0.67）  
- 用户层几乎 **对不上**（目标函数不同：质量加总 ≠ 实体条件分布差）  
- 父 key（merchant_id）还依赖 item_feat 映射，映射率低时大量「假商户」

Dissolve 假设：**层级已经写在表里**，只是把分数「溶」上去。  
多层次下钻若卡在「层级切错 / 层级太粗 / 加密商户不可靠」，dissolve **救不了键**。

### 1.2 Leiden 补的是「键」

Leiden（Traag et al., 2019）在 Louvain 上修了 **连通性 / 局部最优** 问题，常用于：

- 得到 **质量更好、更稳** 的社区划分  
- **多分辨率**（调分辨率参数 \(\gamma\)）→ 粗社区 / 细社区 ≈ 自然「层」  
- 可在 **同一图** 上产出嵌套或并列的 `community_id` 当作 `drill_key`

这对我们的协议恰好是 **M2 接法**（[`MULTILAYER_GRAPH_AND_MM_ATTR.md`](./MULTILAYER_GRAPH_AND_MM_ATTR.md) §B.1）：

> 社区 id → 并列 / 嵌套 `entity_key` → 仍用 \(r_v\) / MMD / stop 表决定下不下钻。

也就是说：Leiden **换钥匙**，不换锁。

### 1.3 和「固定 merchant→user→order」比

| | 固定三层 | Leiden 多分辨率 |
|---|---|---|
| 键从哪来 | 业务 schema / 特征列 | 图结构优化 |
| 层数 | 预登记 3 层 | \(\gamma\) 扫描 → 2–K 层 |
| 坏键风险 | 映射缺失、颗粒不对 | 边定义错、分辨率乱扫 |
| 和 cmean 对齐 | 已有三步脚本 | 社区当 key 后 **同一套** δ/r |
| 讲武德风险 | 低（业务键） | 中（易被说成「发现人群效应」）→ 必须写 J1–J8 |

Leiden **不是**自动更因果，而是 **可能更贴漂移载体的结构切分**——尤其当「商户」只是弱代理时。

---

## 2. Leiden 算法（够用的直觉，不推全文）

### 2.1 目标（模块度家族）

在加权无向图 \(G=(V,E,w)\) 上，找划分 \(\mathcal{C}=\{C_k\}\) 使 **模块度** \(Q\)（或 CPM）大：

- 社区内边密、社区间边疏  
- Louvain：贪心聚合 + 局部移动；可能产生 **内部不连通** 社区  
- **Leiden**：保证社区内部连通，并有更强的局部最优性质；实践中常更稳、更快收敛

### 2.2 多分辨率 =「层」

分辨率参数 \(\gamma\)（越大社区越碎）：

```text
γ 小  →  粗社区  ≈  「L1 板块」
γ 中  →  中社区  ≈  「L2」
γ 大  →  细社区  ≈  「L3 / 近似订单簇」
```

这不是业务上的 merchant/user/order，而是 **结构上的粗→细**。  
下钻时：在粗社区 \(C\) 内再限制子图跑更细 \(\gamma'\)，或直接用预先算好的嵌套划分。

### 2.3 和 Louvain / Infomap / SBM 的取舍（我们场景）

| 方法 | 为何先试 Leiden |
|---|---|
| Louvain | 经典但连通性瑕疵；Leiden 几乎是严格改进 |
| Infomap | 流压缩强，实现/调参更重 |
| SBM | 模型漂亮，算力与辨识更重 |
| GNN 聚类 | 超纲且目标函数更漂 |

**默认提案：Leiden + 分辨率扫描**；Infomap 作对照，不进默认门。

---

## 3. 接到我们协议的正确位置（硬边界）

```text
[上游 · 可选]  构图 → Leiden(γ) → community_id 写入边/实体表
       ↓
[主协议 · 锁定] 图谱特征 X → 标准化 → cmean δ/D/r → stop → K* → FSDS
                 （可选：社区内再 MMD / PO-VIMP）
```

| 允许 | 禁止 |
|---|---|
| 社区 id 当 `drill_key` / 并列 key | 用 \(Q\) / 模块度增益当「要不要下钻」 |
| 社区标签进 \(X\) 的一列（one-hot / embed） | 「该社区转化效应更大」 |
| 比 \(\mathrm{eff}\)、\(\pi\)、残余 \(R\) | 把 Leiden 输出叫 causal subgroup |
| 与固定 merchant key **并列** 比哪把钥匙更好收支撑 | 取代 FSDS / PO-help |

一句话复读（来自多层调研收束）：

> 多层图算法 = 可选的 **结构上游**；归因门只认 \(\delta,R\)。

---

## 4. TencentGR 上怎么构图（落地前必须钉死）

社区质量 **几乎完全取决于图**，不取决于「叫 Leiden 还是 Louvain」。

### 4.1 候选图（由易到难）

| ID | 节点 | 边 | 权重想法 | 备注 |
|---|---|---|---|---|
| G1 | item | covisit / 共现 | 共现计数或 PMI | 已有图谱特征同源；实现近 |
| G2 | user | 共同点击 item / 二部投影 | 投影边权 | 用户下钻直接 |
| G3 | user–item 二部图 | 点击边 | 1 或频次 | Leiden 变体 / 投影后再 Leiden |
| G4 | merchant（若映射可靠） | 经 user 或 item 的投影 | — | 映射率低时慎用 |
| G5 | **时段切片** | W1 / W2 分别构图或差分边 | \(w_{W2}-w_{W1}\) 截断 | 更贴 period-shift；更敏感 |

**建议第一枪：G1（item covisit）+ 可选 G2（user 投影）**，在 **W1∪W2 并集** 或 **W2 侧** 构图；用 period 只作后续 cmean 的 \(W\)，先别把差分边当默认（方差大）。

### 4.2 节点–边–社区如何挂回 grid

每条 localize grid 边（user,item,ts）打上：

- `comm_item_γc` — item 所在粗社区  
- `comm_item_γf` — item 所在细社区  
- `comm_user_γc` — user 社区（若跑 G2）

然后：

```text
drill_key ∈ { merchant_id, comm_item_γc, comm_item_γf, comm_user_γc, user_id, ... }
```

对每个 key 算同一套：实体级 \(r_v\) 或 MMD² → TopK → stop 表 → 比 \(\mathrm{eff}\)。

### 4.3 多分辨率下钻的两种操作型定义

**定义 A — 并列钥匙（推荐先做）**

- 预先算好 \(\gamma\in\{\gamma_c,\gamma_m,\gamma_f\}\) 三套标签  
- 像现在比 merchant vs user：比哪把 key 的行谱更尖、π 更稳  
- **不下钻依赖**：只是换 key，不强制嵌套

**定义 B — 嵌套下钻（更像「层次」）**

1. 在全图用 \(\gamma_c\) 得粗社区，按 \(r\)/`MMD` 选 Top 粗社区 \(C^\star\)  
2. **诱导子图** \(G[C^\star]\) 上用 \(\gamma_f>\gamma_c\) 再 Leiden  
3. 在细社区上重复 stop；不行就停在 \(C^\star\)

B 更贴「多层次」，但也更像 subgroup 流程 → **Justify 必须写满**（见 §6）。

---

## 5. 和现有组件怎么拼（配方）

### 5.1 最小可行实验（MVE）

```text
1. 从 data/tencent_subset 建 G1（item covisit，cap 边数）
2. leiden(G1, γ ∈ {粗,中,细}) → 三列 community id
3. join 到 W1/W2 localized grid
4. 对 key ∈ {merchant_id, comm_γc, comm_γf}：
     - 实体 MMD 或 mean-shift 排名
     - TopK 支撑上跑官方 FSDS / po_help_fsds
     - 记：Jaccard(社区Top vs 商户Top)、W2 HGB、seed σ、π 稳定
5. 判决表：
     - 若 comm key 的 eff / π / W2 不差于 merchant → Leiden key 可并列落地
     - 若只是重排噪声 → 维持「社区超纲」，不进默认
```

### 5.2 成功标准（预先写死，防 fishing）

| 指标 | 过线（建议） |
|---|---|
| 检索 | 粗社区 TopK 与「业务可解释名单」有可述重叠（不必高 Jaccard） |
| 稳定 | 多种子 / 边抽样下社区标签 κ 或 NMI 中等以上；选中社区 π 不崩 |
| 下游 | 社区支撑上 FSDS / PO-help 的 W2 mean **不低于** merchant 三步基线（seed 均值） |
| 讲武德 | 报告零 ATE 用语；社区称「结构块 / 检索键」 |

任一不过 → Leiden **只进附录**，不进默认闭环。

### 5.3 失败模式（预期会踩）

1. **星型 / 头部 item** 吞噬社区 → 需度惩罚或投影阈值  
2. **W1/W2 图差太大** → 并集图社区在单窗上无意义；考虑窗内 Leiden + 匹配  
3. **社区数爆炸** → 细 γ 变成近乎 item 级，退化为 dissolve 的叶子  
4. **rare-pos** → 社区支撑上仍可能 0 正例；下游 AUC 无信息（与 dissolve 相同局限）  
5. **叙事滑坡** → 「Leiden 找到高转化社区」——必须在模板里禁掉

---

## 6. Justify 模板（比 dissolve 更需要）

社区外观比 merchant_id **更像** 因果子群。每份结果页抄这段：

1. **任务**：W1/W2 特征分布漂移的结构键检索，不是处理效应。  
2. **Leiden 输出**：图上的连通密集块 id，不是 treatment arm。  
3. **下钻门**：只使用 \(\delta/D/r/R\)（及可选 MMD），**不用** \(Q\)。  
4. **\(K^\star\)**：仍是 cmean 收完后的支撑，不是「Leiden 核」。  
5. **FSDS / PO-VIMP**：监督排序先验，不是机制证明。  
6. **比较**：Leiden key vs merchant key 是 **哪把检索钥匙更好**，不是哪个人群效应更大。

---

## 7. 和「超纲」声明怎么和解

[`GRAPH_LOCALIZATION_NO_COMMUNITY.md`](./GRAPH_LOCALIZATION_NO_COMMUNITY.md) 写过：社区发现 ❌ 超 scope。

本次提案 **不推翻**「社区不当归因门」，只申请：

| 原条款 | 调整建议 |
|---|---|
| 默认闭环不做 community | **维持** |
| 社区完全不做 | 改为：**可选上游实验**（本 elaboration → 若 MVE 过线再开 `run_leiden_key_*.py`） |
| ULS/PPR/GraphScan | **仍 deferred**（与 Leiden 不同族） |

即：**Leiden = 受控试验的钥匙工厂**；默认路径仍是图谱特征 + cmean。

---

## 8. 实施草图（你点头后再写代码）

```text
scripts/tencent_gr/
  build_covisit_graph.py          # G1 边表
  run_leiden_community_keys.py    # γ 扫描 → parquet 标签
  run_leiden_vs_merchant_drill.py # 并列 key：MMD/cmean/FSDS/π
docs/tencent_gr/
  LEIDEN_FOR_DRILL_ELABORATION.md # 本文
```

依赖候选：`python-igraph` + `leidenalg`（或 `networkx` + `cdlib`；优先 igraph/leidenalg，成熟稳定）。

不算进默认 CI；先 smoke 在 `data/tencent_subset` + 现有 localized grid。

---

## 9. 给你审的结论页

1. **Dissolve** 适合「已有父 key 的分数上卷」；对「多层次下钻」的 **键发现** 帮助有限。  
2. **Leiden** 更对口：多分辨率社区 ≈ 可下钻层级键；应作为 **上游 key**，归因门仍是 cmean/stop/FSDS。  
3. **第一枪** 应用 item-covisit 图 + 并列 key 对比 merchant；嵌套子图下钻放第二枪。  
4. **成功标准** 预先写死；不过线就继续保持「社区不进默认」。  
5. 全程讲武德：社区 ≠ 效应子群。

若认可 §5 MVE 与 §5.2 过线标准，下一步再开实现 PR；若只想要「嵌套定义 B」或换 G2 用户图，可以指定我改实验设计后再动代码。
