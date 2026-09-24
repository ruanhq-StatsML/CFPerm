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

### 3.1 Pairwise 切窗（就是你说的这个）

边计数 `e_n_exp` 等是 **特征 X**，不是切窗单位。切窗单位是 **边行**：

```python
df = df.sort_values("e_last_ts").reset_index(drop=True)
# chunk ids: 0,0,...,0, 1,1,...,1, ...
from itertools import pairwise
N = 1000
n_chunks = len(df) // N   # (+ maybe a tail)
for t0, t1 in pairwise(range(n_chunks)):
    df1 = df.iloc[t0 * N : (t0 + 1) * N]   # train
    df2 = df.iloc[t1 * N : (t1 + 1) * N]   # test
    # fit SelectKBest+HGB on (X1, y1); score on (X2, y2)
```

代码里 `adjacent_chunk_pairs` = `list(pairwise(occupied_chunk_ids))`，语义同上。

### 3.2 ops ≈35 维 vs content ≈14 维 —— 进模逻辑

**共同前提**：两边都先丢掉 convert leak（`e_n_cnv`, `u_cvr`, `i_n_cnv`…），否则 Y 直接漏进 X。

| | ops | content |
|---|---|---|
| 问题 | 「量 + 结构」合在一起，下一段还能不能排转化？ | 去掉量之后，**路径/归因结构**还能否传？ |
| 典型 top | `e_log1p_exp`, `e_n_exp` | `i_credit_last`, `i_share_last` |
| 业务读 | 买量强度基线 | 末跳/路径候选 shortlist |
| 误读 | AUC≈1 ⇒ 可上线 | credit tip ⇒ 定罪末跳 |

content 十四维大致三类：

1. **路径归因**（核心）：`i_credit_{first,last,linear}`, `i_share_{first,last,linear}`, `i_log1p_credit_linear`, `i_item_credit_rank`
2. **共现/结构**：`i_n_covisit_neighbors`, `i_log1p_n_covisit_neighbors`, `ui_pop_mismatch`
3. **用户非量形状**：`u_span_sec`, `u_n_uniq_items`, `u_log1p_n_uniq_items`（会话跨度/多样性，不是曝光次数）

### 3.3 `i_credit_*` / `i_share_*` 是什么

来自转化路径上的 **触点归因**（first / last / linear），再归一成份额：

- `credit_last`：作为路径**末跳**拿到的归因质量  
- `share_last`：`credit_last / Σ credit_last`（该 item 在末跳归因池里的份额）  
- `*_linear`：路径线性分摊；`*_first`：首跳  

广告读法：content 面板里它们进 top = 「去量后，末跳/路径结构仍能帮着排转化」→ **L1 末跳候选**，对齐审出 tip 桶；**不是**自动限投证明。

---

## 4. 还有其他的吗？

| 还有 | 和本边 pairwise 的关系 |
|---|---|
| **日历宽窗**（5/10/60 min） | 另一套 grain（`run_diffusiondb_adjacent_board`）；投放突发时不等价于边数窗 |
| **其它 pack** | DiffDB=prompt 行；Metro/Beijing=小时行；Waymo=proxy 行——同一 pairwise 协议，不同一行含义 |
| **重叠/stride** | 当前 **不重叠** 相邻块；要更密可 `stride=N/2`（未做） |
| **原子 PV 切窗** | 需回 `seq` 曝光事件表；本 `feature_grid` 做不到 |

本广告场景默认：**边行 + pairwise + ops/content 双面板**。

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
