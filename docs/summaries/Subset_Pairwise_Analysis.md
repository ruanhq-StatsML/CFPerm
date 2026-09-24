# Subset pairwise + 日历宽窗 + pack 指标

曝光强度 = **买量基线**。有了它，subset-analysis 很直接：  
同一套 `sort(e_last_ts) → pairwise(N)`，只换 **边的子集**。

---

## 1. Subset-analysis 怎么切

| cohort | 切法 | 在问什么 |
|---|---|---|
| `all` | 全边流 | 参照 |
| `intensity_high/low` | `i_log1p_n_exp` 三分位（**物品曝光规模**；边 `e_n_exp` 多数为 1 切不开） | 买量基线是否吃「大曝光创意」 |
| `credit_high/low` | `i_item_credit_rank` 头/尾 25%（秩越小归因份额越高） | 末跳归因重的边，content 是否更干净 |

每个 cohort × `{ops, content}` 各跑一遍相邻 pairwise（默认 N=1000）。

读法：

1. **ops @ intensity_high ≈ ops @ intensity_low ≈ 0.99** → 买量强度基线在大/小曝光创意上都稳传（cohort-robust）  
2. **content 两端都 ~0.8** → 去量后仍有可传结构，不只是「有量才行」  
3. **credit_high content** → 应看到 credit/share 类 top；AUC 可略低于全量 content  
4. **credit_low y_rate→0 / 无 pairs** → 归因秩最差的边几乎无转化，pairwise FS 自然空——这是子集结论，不是脚本坏了  
5. Pack scorecard / 日历 scorecard 分表看，不混 grain

```bash
PYTHONPATH=. python3 scripts/run_subset_pairwise_analysis.py \
  --chunk-size 1000 --out results/subset_pairwise_analysis
```

---

## 2. Pack pairwise 指标（样本数 grain）

来自 `sample_chunk_adjacent_board` cross_pack — **同一协议、不同「一行」**：

| pack | 一行 | 典型 pairwise 读数 |
|---|---|---|
| tencent_gr (ops) | UI 边 | AUC≈0.99，Jac 高 → 强度基线 |
| tencent_gr_content | 同边、去量 | AUC 掉一截，top→credit |
| diffusiondb | prompt 行 | AUC≈0.60，弱传 |
| metro | 小时行 | 高 AUC + 低 Jac → 驱动在换 |
| beijing / waymo | 小时 / proxy | 高传 + 较稳 |

指标固定看：`mean_auc`, `mean_logreg_auc`, `mean_jaccard`, `mean_delta_Y`, `FSDS∩cmean`。  
**不要**和日历宽窗表混在一个格子里比绝对值。

---

## 3. 日历宽窗 pairwise（另一 grain）

DiffusionDB `adjacent_board`：`width_min` 切日历箱 → 相邻箱 train→test。

| | 样本数 pairwise | 日历宽窗 pairwise |
|---|---|---|
| 平衡什么 | 每窗 **n** | 每窗 **时间宽度** |
| 副作用 | 日历跨度狂跳（边流突发） | 每箱 n 失衡、空箱 |
| 广告默认 | **用这个**（边 N=1000/2000） | 作对照，不替代 |
| DiffDB | sample-chunk 板 | width 30/60/180 min 板 |

Equal-count 的日历副作用本机可见：同 N=1000 边，span 可从 ~1100h → ~200h。  
那不是 bug，是「为什么广告用边数窗」的证据。

---

## 4. 还有其他的吗？

| 还可做 | 状态 |
|---|---|
| intensity × credit 交叉四格 | 易加；先看边缘三分位 |
| stride=N/2 重叠 pairwise | 未做 |
| 按广告主 merchant 子集 | 需 `item_feat.122` 映射 |
| PV 原子事件 pairwise | 要回 seq，不是 feature_grid |

**切片列注意**：边 `e_n_exp` 绝大多数为 1，不能当买量三分位；用 **`i_log1p_n_exp`**。末跳用 **`i_item_credit_rank`** 头尾，比 raw `credit_last`（大量 0）可分。

当前交付：强度/credit 子集 pairwise + pack scorecard + 日历对照。

产物：`results/subset_pairwise_analysis/`  
相关：[`Tencent_Edge_Grain_Columns.md`](./Tencent_Edge_Grain_Columns.md) · [`Ad_Funnel_Scenario_Board.md`](./Ad_Funnel_Scenario_Board.md)
