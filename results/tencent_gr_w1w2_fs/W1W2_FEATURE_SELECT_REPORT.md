# W1 vs W2 → feature selection (落地)

比较两个时间窗的图谱特征，**用差异做 feature selection**。
选出的是时段漂移驱动特征，可直接进监控 / 重训 / 归因。

## Protocol
1. `feature_engineer(t_start,t_end)` ×2，gap = **30** days (≥30)
2. Localization subset（可选）: W1 `share_linear` top-300
3. **Compare W1 vs W2**: |SMD| + domain-AUC + KS + period-HGB perm-importance
4. Rank-average → select top-**15** (drop within-window `*_rank`)

## Windows
- span **231.3**d | W1 n=2051 | W2 n=1788 | feats=21
- period classifier AUC = **1.000** (how separable are the two windows)
- click rate snapshot: W1=0.0000 → W2=0.0000

## Selected shift drivers
| rank | feature | \|SMD\| | domain-AUC | KS | perm-imp |
|---:|---|---:|---:|---:|---:|
| 1 | `i_share_linear` | 2.136 | 1.000 | 1.000 | 0.4530 |
| 2 | `u_n_events` | 2.858 | 0.956 | 0.810 | 0.0000 |
| 3 | `u_n_uniq_items` | 2.855 | 0.956 | 0.811 | 0.0000 |
| 4 | `u_log1p_n_events` | 2.195 | 0.956 | 0.810 | 0.0000 |
| 5 | `u_log1p_n_uniq` | 2.193 | 0.956 | 0.811 | 0.0000 |
| 6 | `u_n_exp` | 2.759 | 0.953 | 0.810 | 0.0000 |
| 7 | `i_log1p_n_covisit` | 1.721 | 0.909 | 0.725 | 0.0000 |
| 8 | `i_n_covisit_neighbors` | 1.599 | 0.909 | 0.725 | 0.0000 |
| 9 | `i_log1p_n_users` | 1.639 | 0.883 | 0.688 | 0.0000 |
| 10 | `i_n_users` | 1.454 | 0.883 | 0.688 | 0.0000 |
| 11 | `i_log1p_n_exp` | 1.617 | 0.881 | 0.655 | 0.0000 |
| 12 | `i_n_exp` | 1.463 | 0.881 | 0.655 | 0.0000 |
| 13 | `ui_pop_mismatch` | 0.231 | 0.835 | 0.756 | 0.0000 |
| 14 | `i_share_last` | 0.864 | 0.675 | 0.618 | 0.0000 |
| 15 | `i_share_first` | 0.969 | 0.648 | 0.616 | 0.0000 |

## Why this lands
- 输入是业务已有的 `(user,item)` 图谱特征 + 时间窗
- 输出是 **W1→W2 漂移特征清单**，可挂告警 / 重训触发 / 运营解释
- 无图算法；防泄漏：两窗独立 FE，选择信号是 period W 而非 peek 未来标签训练
