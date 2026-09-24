# 按边切窗：颗粒度与 columns

> 广告看板说的「每 1000/2000 条切窗」= **按 `feature_grid` 的边行切**，不是按日历小时、不是按 user、不是按 impression 原子事件。

表：`results/tencent_gr_tabular_ui/feature_grid.parquet`  
构建：`scripts/tencent_gr/run_tabular_user_item_feats.py` → `build_feature_grid`

---

## 1. 一行是什么？

| 主张 | 事实 |
|---|---|
| **一行** | 一条 **(user_id, item_id)** UI 边的聚合行 |
| **主键** | `(user_id, item_id)` — 表内 **无重复** |
| **行数** | 286,949 = `n_edges_ui`（≠ 原始 seq 里每一次曝光原子事件） |
| **时间序** | `e_last_ts`（该边上最后一次交互 unix 秒） |
| **标签 Y** | `y_convert` = 该边是否出现过转化（`has_convert`） |

所以「按边切窗」=：

```text
sort by e_last_ts
chunk_id = row_index // N     # N ∈ {1000, 2000}
```

相邻 chunk \(t\to t{+}1\) 各约 N 条 **边**，日历跨度 **不保证相等**（本机前几窗可差到 ~1000h vs ~200h）——这正是用样本量而不用「每小时」的理由：投放日志突发。

---

## 2. 列族（52 列 → 切窗用哪些）

### 2.1 键 / 序 / 标签（不进 X，只管切窗与 Y）

| column | 角色 |
|---|---|
| `user_id` | 用户键 |
| `item_id` | 创意/商品键（上游 merchant≈广告主在 `item_feat.122`） |
| `e_last_ts` | **排序键**（边粒度时间） |
| `y_convert` | **Y**：边是否转化 |

### 2.2 边层 `e_*`（这条 UI 边上的计数）

| column | 含义 |
|---|---|
| `e_n_exp` / `e_log1p_exp` | 边上曝光次数 |
| `e_n_clk` / `e_log1p_clk` | 边上点击次数 |
| `e_n_cnv` | 边上转化次数（**leak → 两面板都丢**） |
| `e_ctr` | 边 CTR |

### 2.3 用户层 `u_*`（挂到边上的用户图谱聚合）

活跃/规模：`u_n_events`, `u_n_exp`, `u_n_clk`, `u_n_uniq_items`, `u_span_sec`, `u_log1p_*`, `u_user_activity_rank`  
转化相关（**leak**）：`u_n_cnv`, `u_cvr`, `u_ctcvr`, `u_has_convert`, `u_log1p_n_cnv`, `u_user_ctcvr_rank`  
比率：`u_ctr`

### 2.4 物品/创意层 `i_*`

规模：`i_n_exp`, `i_n_clk`, `i_n_users`, `i_ctr`, `i_item_pop_rank`, `i_log1p_n_*`  
转化 leak：`i_n_cnv`, `i_cvr`, `i_ctcvr`, `i_n_as_convert_terminal`, `i_item_cnv_rank`, `i_log1p_n_cnv`  
**路径归因（content 核心）**：`i_credit_{first,last,linear}`, `i_share_{first,last,linear}`, `i_item_credit_rank`, `i_log1p_credit_linear`  
共现：`i_n_covisit_neighbors`, `i_log1p_n_covisit_neighbors`  
错配：`ui_pop_mismatch`

---

## 3. 切窗时两套 X（同一套边行）

同一 `sort(e_last_ts)` + `// N`，只换进模特征：

| 面板 | 保留 | 丢掉 |
|---|---|---|
| **ops** | 边/用户/物品 **量与 CTR** + 归因/共现等非 convert 列（约 35 维） | 全部 convert leak（`*cnv*`, `*cvr*`, `y_convert`） |
| **content** | 归因/共现/错配等（约 14 维） | leak **再加** volume（`e_n_exp`, `e_log1p_exp`, `i_n_exp`, `u_n_exp`, `*_ctr`, pop rank…） |

漏斗读法：

- ops top 常是 `e_log1p_exp` → **买量强度在传**
- content top 常是 `i_credit_last` / `i_share_last` → **末跳候选**

---

## 4. 和「原子曝光日志」差在哪

| | 原子 seq 事件 | 本表边行 |
|---|---|---|
| 一行 | 一次曝光/点击 | 一个 (u,i) 聚合 |
| 切 N=1000 | 1000 次曝光 | 1000 条 UI 边 |
| 适合 | 精细时序 CTR | 审出/归因/图谱 tip |

若业务要「按曝光 PV 切窗」，需要回 `seq` 原子表另建；**当前广告看板明确是边颗粒度**。

---

## 5. 机器可读清单

跑看板或：

```bash
PYTHONPATH=. python3 scripts/export_tencent_edge_grain_schema.py
```

写出 `results/sample_chunk_adjacent_board/edge_grain_schema.json`。

相关：[`Ad_Funnel_Scenario_Board.md`](./Ad_Funnel_Scenario_Board.md)
