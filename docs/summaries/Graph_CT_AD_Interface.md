# Continuous-time graph AD：可做事项 + 数据接口

> Drill 固定：**MMD + PO + cmean**（rank-average），不改门。  
> 监测骨架：与 Grad-OnlineRFPerm 同构的 **OnlineRFPerm 标量流**。  
> 图 = 每个时间片上的 \((u,i,t)\) 边集（+ 可选社区键）。

相关：[`Graph_Localize_Algorithm_Elaboration.md`](Graph_Localize_Algorithm_Elaboration.md)、
[`../tencent_gr/W1W2_MMD_PO_LOCALIZE_FSDS.md`](../tencent_gr/W1W2_MMD_PO_LOCALIZE_FSDS.md)。

---

## 1. 一句话分工

| 层 | 回答 | 方案 |
|---|---|---|
| **When** | 这一时刻图异常了吗？ | 连续监测：\(T_t^{\mathrm{MMD}}, T_t^{\mathrm{PO}}, T_t^{\mathrm{cmean}}\) → OnlineRFPerm → \(R_t\) |
| **Where** | 异常在哪些实体？ | Drill：`rank-average(MMD, PO, cmean)` → merchant→user→order → \(K^\star\) |
| **What tip** | 哪些特征 tip？ | \(K^\star\) 上 **一次** FSDS（可选 PO-VIMP prior） |
| **Graph algos** | 大规模图结构键从哪来？ | Leiden / 多跳 / 嵌入 —— **只作 upstream key**，不进 Drill 门 |

---

## 2. 连续监测框架里「能做什么」（用我们的方案）

### A. 长期 running monitor（必做主线）

1. **时间切片图流**  
   \(\Delta t\)（建议 1d / 6h / 1h）→ \(G_t=(V_t,E_t)\)，边来自 `[t, t+\Delta t)` 的 `(user,item,ts)`。  
   滚动参考窗 \(G_{t-h:t-1}\) 或 burn-in 固定 \(G_{\mathrm{ref}}\)。

2. **三路标量 OnlineRFPerm（并行，不互替）**
   \[
   T_t^{\mathrm{MMD}}=\widehat{\mathrm{MMD}}^2(X_{t-1},X_t)-e_{\mathrm{MMD}},\quad
   T_t^{\mathrm{PO}}=\mathrm{PO\text{-}gap}(t-1\!\to\!t)-e_{\mathrm{PO}},\quad
   T_t^{\mathrm{cmean}}=\|\mu_t-\mu_{\mathrm{ref}}\|_2-e_{\mathrm{cm}}.
   \]
   各自 EWMA-\(p\) + alpha-investing → \(R_t^{\bullet}\)。  
   **Fire** 可用并集或「MMD∧(PO∨cmean)」工程规则。

3. **Reject 才 Drill**  
   \(R_t=1\) → 在 \(S_t=E_{t-h:t}\)（或 \(E_t\) vs ref）上跑  
   `rank-average(cmean_l2, mmd2, po_tau2)` 实体排序 → 嵌套 merchant→user→order → 停层 \(K^\star\) → 一次 FSDS。  
   \(R_t=0\)：**不下钻、不跑 FSDS**（省算力）。

4. **Signed 方向卡片**  
   每次交付：\(\Delta\bar y\)、\(\mathrm{sign}(\delta_j)\)、实体 \(s_v\)（正向/负向）。

5. **Lead / AR 仪表**  
   Lead = \(t_{\mathrm{MMD}}-t_{\mathrm{PO}}\) 等；AR = 报警率（非 Type-I FAR）。

### B. 结合大规模图算法（可选 upstream，不进门）

| 图算法产出 | 怎么接我们的方案 | 不要做什么 |
|---|---|---|
| **Leiden / Louvain 社区** \(c(v)\) | 作为与 merchant 并行的 `drill_key`；比较 \(\mathrm{eff}/\pi\) | 用 modularity \(Q\) 当 Drill 门 |
| **多跳 / 共现图**（user–item–user） | 进 FE：`i_n_covisit_neighbors` 等；或作边权重进入 MMD 的 \(X\) | 用路径显著性替代 MMD/PO/cmean |
| **图嵌入**（Node2Vec / GNN / mm_emb） | 拼进边特征 \(x_e\)；cmean/MMD 在嵌入空间跑 | 单独 unsupervised AD 覆盖 OnlineRFPerm |
| **PageRank / 中心性** | 作实体 prior / 体积代理；或 rank blend 的第四列（可选） | 进 \(\mathrm{Drill}\) 产品式 |
| **动态图摘要**（sketch / reservoir） | 大规模下近似 \(X_t\)、近似 MMD | 改变 PO/cmean 语义定义 |

**组合拳：**  
`大规模图算法 → 结构键 / 富特征` **⊕** `MMD+PO+cmean Drill` **⊕** `OnlineRFPerm 连续监测`。

### C. 工程可落地的下一刀（优先级）

1. **Stream harness**：`Δt` 循环 + 三路 \(T_t\) 日志 + reject 触发 localize（接现有 `run_w1w2_mmd_po_localize_fsds`）。  
2. **Drill stop 埋点**：每层记 \((\|\delta\|,\mathrm{Tail},\mathrm{CV},R,\pi)\)，仅 `Drill=0` 或业务 cap 时 FSDS。  
3. **Leiden 旁路对比**：同一 reject 时刻，merchant-key vs community-key 的 \(\mathrm{eff}\)。  
4. **mm_emb 模态**：把 `emb_81_32` 等拼进 \(X\)，看 tip 是否从行为特征迁到视觉/文本。  
5. **Confirm timing（可选）**：tip-cmean 幅度 × LOCO confirm，挂在 reject 之后。

---

## 3. 数据接口（给你用的）

### 3.1 磁盘布局（TencentGR subset）

```
data/tencent_subset/
  seq/*.parquet          # 行为流（建图边）
  item_feat/*.parquet    # 物品侧（含 merchant≈col 122）
  user_feat/*.parquet    # 用户侧
  README.md
```

可选全量：`candidate/`、`mm_emb/emb_*/**/*.parquet`、`indexer.pkl`。

### 3.2 原始 schema

**`seq`**（每用户一条）：

| field | type | 含义 |
|---|---|---|
| `user_id` | int64 | RID |
| `seq` | List\[Dict\] | 每元素：`item_id`(int), `action_type`(0=曝光,1=click), `timestamp`(int unix) |

展开后的边：

```text
Edge = (user_id, item_id, action_type, timestamp)
```

**`item_feat`**：`item_id` + 加密列 `100`…`122`；**商户代理 = `122`**。  
**`user_feat`**：`user_id` + `103`…`110`（部分为 List）。

### 3.3 我们代码侧接口（已有）

```python
from pathlib import Path
from scripts.tencent_gr.time_window_feats import (
    scan_time_range,          # (t_min, t_max, n_users)
    propose_two_windows,      # 原型：两窗 + gap（对照用）
    feature_engineer,         # [t_start, t_end) → grid / edges / meta
    fsds_feature_columns,
)
from scripts.tencent_gr.run_three_step_subset_localize import (
    load_item_merchant_map,   # item_id → merchant_id (col 122)
    attach_merchant,
)

ROOT = Path("data/tencent_subset")

# 1) 时间轴
t_min, t_max, n = scan_time_range(ROOT, max_users=20_000)

# 2) 任一时间片上的「图面板」（连续监测的一帧）
#    G_t 的边特征 / 标签都在返回的 grid、edge_df 里
panel = feature_engineer(
    ROOT, t_start, t_end,
    max_users=20_000,
    window_name=f"slice_{t_start}",
    terminal_action=1,   # TencentGR: click = terminal success
)
grid = panel["grid"]          # 一行 ≈ (user,item) 聚合边 + u_/i_/e_ 特征
edge_df = panel["edge_df"]    # 更细的边表（含 last_ts）
meta = panel["meta"]          # n_days, counts, ...

# 3) 挂商户键（Drill L1）
imap = load_item_merchant_map(ROOT, merchant_col="122")
grid = attach_merchant(grid, imap)   # + merchant_id, merchant_mapped

# 4) FSDS 可用列（已去泄漏）
cols = fsds_feature_columns(grid)
```

`feature_engineer` 返回字典：

| key | 内容 |
|---|---|
| `grid` | `(user_id, item_id)` 特征行 + `y_convert` + `e_last_ts` |
| `edge_df` | 窗内边聚合 |
| `user_df` / `item_df` | 用户 / 物品汇总 |
| `convert_df` | 三段路径归因行 |
| `covisit_df` | 共现邻居 |
| `meta` | 窗名、天数、计数 |

### 3.4 连续监测建议接口（长期 running）

伪接口（实现时可新开 `scripts/tencent_gr/stream_graph_monitor.py`）：

```python
@dataclass
class GraphSlice:
    t: int                 # slice index or t_start
    t_start: int
    t_end: int
    grid: pd.DataFrame     # X rows for MMD/PO/cmean
    edges: pd.DataFrame    # raw (u,i,a,ts) or aggregated
    merchant_id: Optional[pd.Series]
    community_id: Optional[pd.Series]  # Leiden etc. optional

def iter_graph_slices(
    root: Path,
    *,
    delta_sec: int = 86400,   # Δt
    max_users: int = 20_000,
    t_min: Optional[int] = None,
    t_max: Optional[int] = None,
) -> Iterator[GraphSlice]:
    """Yield one graph panel per [t, t+Δt)."""
    ...

def monitor_step(
    prev: GraphSlice,
    cur: GraphSlice,
    state_mmd, state_po, state_cm,
) -> dict:
    """Return {T_mmd, T_po, T_cm, R_mmd, R_po, R_cm, fire}."""
    ...

def drill_on_reject(
    ref: GraphSlice,
    cur: GraphSlice,
    *,
    keys=("merchant_id", "user_id"),  # + optional community_id
) -> dict:
    """rank-average(MMD, PO, cmean) → K* → one FSDS."""
    ...
```

**边级最小契约（喂图算法 / 喂我们）：**

```python
# long format — 大规模图库（NetworkX / igraph / cuGraph / PyG）都能吃
edges_long = pd.DataFrame({
    "src": user_id,       # or bipartite left
    "dst": item_id,       # bipartite right
    "ts": timestamp,
    "action": action_type,
    "y": (action_type == 1).astype(int),  # click
})
# 可选节点属性
item_attr = item_feat[["item_id", "122"]].rename(columns={"122": "merchant_id"})
user_attr = user_feat  # as-is
```

### 3.5 CLI 入口（现状）

```bash
# 两窗原型（对照，非连续）
PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \
  --root data/tencent_subset

# 三层实体下钻原型
PYTHONPATH=. python3 scripts/tencent_gr/run_three_step_subset_localize.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30
```

连续 `iter_graph_slices` harness = 下一步工程（上表 C.1）。

---

## 4. 与 W1/W2 原型的关系

| | W1/W2 原型 | 连续监测 |
|---|---|---|
| 时间 | 两窗 + ≥30d gap | 每个 \(\Delta t\) 一张图 |
| 触发 | 人工跑 | OnlineRFPerm \(R_t\) |
| Drill | MMD+PO+cmean | **同一套**，仅在 reject 时 |
| 图算法 | 无 / 旁路 | 可选 community / emb upstream |

原型证明 **where**；连续框架解决 **when + 长期 running**。

---

## 5. 一页 takeaway

1. 监测三路：**MMD / PO / cmean** OnlineRFPerm；Drill 仍是三者 rank-average。  
2. 大规模图算法 → **键与特征**，不替代门。  
3. 数据从 `seq` 展成 `(u,i,a,t)`，经 `feature_engineer(t0,t1)` 得每片 `grid`；商户键 `item_feat.122`。  
4. 下一刀：把 `feature_engineer` 打成 `iter_graph_slices` + reject→drill harness。
