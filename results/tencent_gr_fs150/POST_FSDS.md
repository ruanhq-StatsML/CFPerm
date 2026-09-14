# 转化粒 FSDS（User / Ctx / 交叉）

Y = 满窗 `y_post_clk_1d`。X 停在 `cnv_ts`。SKU 闸关，同品列不进。
不是预报榜。RF-domain = P(X)；PO/LOGO/LOCO = 哪包/哪列扛早/晚对 Y 的差。

daily mix：n_cnv=12866  sku_clk=0.0011  empty=0.484  sess0=0.984  GATE=0

中位切：塔和列的 LOGO/LOCO **全负**。Q1 vs Q4 家族 LOGO 仍负；LOCO ΔR>0 的是 `n_clk_before_7d`、`dt_any_min`、`n_clk_7d × sess`，量级 1e-5，R 仍是 3e-4。
概念没跳（和用户粒 rfperm fire=0 同一句）。RF-domain 分的是人（`n_prior_cnv`）。PO-VIMP 爱交叉，LOCO 不认——impurity ≠ 漂移源。SKU 两批都 ~0。

## 时钟 `seq_t_end_median`

n=12564 pos=0.113 W1=0.500  RF-domain AUC **0.680**  PO-risk R **0.000330**  τ(PO, RF)=0.554

构成（early vs late）：

| | n | pos | empty_any | sess0 | sku_clk |
|---|---:|---:|---:|---:|---:|
| early | 6286 | 0.077 | 0.493 | 0.987 | 0.0005 |
| late | 6278 | 0.150 | 0.489 | 0.981 | 0.0018 |

| tower | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| cross | 6 | 0.167 | 0.416 | -0.00004 | 0.000 |
| ctx | 6 | 0.236 | 0.210 | -0.00002 | 0.000 |
| user | 6 | 0.597 | 0.374 | -0.00006 | 0.000 |

PO-VIMP 序：

| rank | feat | tower | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `uc_n_clk_7d_log__x__dt_any_log` | cross | 0.1941 | 0.0929 |
| 2 | `uc_lag_post_clk_1d_rate__x__dt_any_log` | cross | 0.1808 | 0.0441 |
| 3 | `dt_any_min` | ctx | 0.1571 | 0.1420 |
| 4 | `n_prior_cnv` | user | 0.1211 | 0.3494 |
| 5 | `lag_empty_any` | user | 0.0840 | 0.0830 |
| 6 | `lag_post_clk_1d_rate` | user | 0.0723 | 0.0393 |
| 7 | `n_clk_before_7d` | user | 0.0523 | 0.0841 |
| 8 | `n_clk_before_1d` | user | 0.0398 | 0.0410 |
| 9 | `sess_pos` | ctx | 0.0292 | 0.0111 |
| 10 | `log1p_price` | ctx | 0.0227 | 0.0638 |
| 11 | `uc_n_clk_7d_log__x__sess_clk_before` | cross | 0.0210 | 0.0026 |
| 12 | `uc_lag_post_clk_1d_rate__x__sess_clk_before` | cross | 0.0198 | 0.0010 |
| 13 | `n_clk_before_1h` | user | 0.0045 | 0.0006 |
| 14 | `sess_clk_before` | ctx | 0.0009 | 0.0004 |
| 15 | `empty_any` | ctx | 0.0003 | 0.0083 |
| 16 | `dt_any_min_miss` | ctx | 0.0000 | 0.0102 |
| 17 | `uc_lag_post_clk_1d_rate__x__empty_any` | cross | 0.0000 | 0.0055 |
| 18 | `uc_lag_empty_any__x__empty_any` | cross | 0.0000 | 0.0208 |

LOCO ΔR（>0 = 这列在扛 PO-risk；impurity 头名可以是负的）：

| feat | tower | impurity | LOCO ΔR |
|---|---|---:|---:|
| `sess_clk_before` | ctx | 0.0009 | -4.68e-06 |
| `uc_lag_empty_any__x__empty_any` | cross | 0.0000 | -5.54e-06 |
| `uc_lag_post_clk_1d_rate__x__sess_clk_before` | cross | 0.0198 | -6.00e-06 |
| `lag_empty_any` | user | 0.0840 | -1.01e-05 |
| `n_clk_before_1d` | user | 0.0398 | -1.27e-05 |
| `sess_pos` | ctx | 0.0292 | -1.79e-05 |
| `uc_lag_post_clk_1d_rate__x__dt_any_log` | cross | 0.1808 | -1.85e-05 |
| `lag_post_clk_1d_rate` | user | 0.0723 | -2.12e-05 |
| `empty_any` | ctx | 0.0003 | -2.43e-05 |
| `n_clk_before_7d` | user | 0.0523 | -2.49e-05 |
| `n_prior_cnv` | user | 0.1211 | -2.65e-05 |
| `dt_any_min_miss` | ctx | 0.0000 | -3.20e-05 |
| `uc_lag_post_clk_1d_rate__x__empty_any` | cross | 0.0000 | -3.61e-05 |
| `log1p_price` | ctx | 0.0227 | -3.67e-05 |
| `uc_n_clk_7d_log__x__dt_any_log` | cross | 0.1941 | -4.11e-05 |
| `n_clk_before_1h` | user | 0.0045 | -4.49e-05 |
| `uc_n_clk_7d_log__x__sess_clk_before` | cross | 0.0210 | -5.37e-05 |
| `dt_any_min` | ctx | 0.1571 | -5.52e-05 |

## 时钟 `seq_t_end_q1q4`

n=6286 pos=0.095 W1=0.500  RF-domain AUC **0.732**  PO-risk R **0.000362**  τ(PO, RF)=0.466

构成（early vs late）：

| | n | pos | empty_any | sess0 | sku_clk |
|---|---:|---:|---:|---:|---:|
| early | 3142 | 0.056 | 0.491 | 0.989 | 0.0003 |
| late | 3144 | 0.135 | 0.524 | 0.981 | 0.0013 |

| tower | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| cross | 6 | 0.100 | 0.551 | -0.00004 | 0.000 |
| ctx | 6 | 0.382 | 0.174 | -0.00011 | 0.000 |
| user | 6 | 0.518 | 0.275 | -0.00010 | 0.000 |

PO-VIMP 序：

| rank | feat | tower | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `uc_lag_post_clk_1d_rate__x__dt_any_log` | cross | 0.2995 | 0.0288 |
| 2 | `uc_n_clk_7d_log__x__dt_any_log` | cross | 0.1836 | 0.0599 |
| 3 | `dt_any_min` | ctx | 0.1221 | 0.2222 |
| 4 | `n_prior_cnv` | user | 0.0988 | 0.3313 |
| 5 | `n_clk_before_1d` | user | 0.0606 | 0.0190 |
| 6 | `lag_empty_any` | user | 0.0597 | 0.0695 |
| 7 | `uc_n_clk_7d_log__x__sess_clk_before` | cross | 0.0542 | 0.0012 |
| 8 | `log1p_price` | ctx | 0.0384 | 0.1032 |
| 9 | `lag_post_clk_1d_rate` | user | 0.0271 | 0.0382 |
| 10 | `n_clk_before_7d` | user | 0.0207 | 0.0579 |
| 11 | `sess_pos` | ctx | 0.0139 | 0.0223 |
| 12 | `uc_lag_post_clk_1d_rate__x__sess_clk_before` | cross | 0.0137 | 0.0005 |
| 13 | `n_clk_before_1h` | user | 0.0076 | 0.0019 |
| 14 | `dt_any_min_miss` | ctx | 0.0000 | 0.0180 |
| 15 | `sess_clk_before` | ctx | 0.0000 | 0.0008 |
| 16 | `empty_any` | ctx | 0.0000 | 0.0156 |
| 17 | `uc_lag_post_clk_1d_rate__x__empty_any` | cross | 0.0000 | 0.0013 |
| 18 | `uc_lag_empty_any__x__empty_any` | cross | 0.0000 | 0.0083 |

LOCO ΔR（>0 = 这列在扛 PO-risk；impurity 头名可以是负的）：

| feat | tower | impurity | LOCO ΔR |
|---|---|---:|---:|
| `n_clk_before_7d` | user | 0.0207 | +3.75e-05 |
| `uc_n_clk_7d_log__x__sess_clk_before` | cross | 0.0542 | +3.20e-05 |
| `dt_any_min` | ctx | 0.1221 | +2.49e-05 |
| `n_clk_before_1h` | user | 0.0076 | +1.59e-05 |
| `uc_lag_post_clk_1d_rate__x__dt_any_log` | cross | 0.2995 | +1.20e-05 |
| `lag_post_clk_1d_rate` | user | 0.0271 | +7.11e-06 |
| `empty_any` | ctx | 0.0000 | -1.44e-05 |
| `uc_lag_empty_any__x__empty_any` | cross | 0.0000 | -1.61e-05 |
| `lag_empty_any` | user | 0.0597 | -2.38e-05 |
| `uc_lag_post_clk_1d_rate__x__empty_any` | cross | 0.0000 | -2.72e-05 |
| `uc_n_clk_7d_log__x__dt_any_log` | cross | 0.1836 | -2.74e-05 |
| `uc_lag_post_clk_1d_rate__x__sess_clk_before` | cross | 0.0137 | -3.32e-05 |
| `n_prior_cnv` | user | 0.0988 | -4.02e-05 |
| `n_clk_before_1d` | user | 0.0606 | -4.05e-05 |
| `log1p_price` | ctx | 0.0384 | -5.49e-05 |
| `dt_any_min_miss` | ctx | 0.0000 | -6.68e-05 |
| `sess_clk_before` | ctx | 0.0000 | -7.04e-05 |
| `sess_pos` | ctx | 0.0139 | -1.10e-04 |

## 时钟 `cnv_ts_median`

n=12564 pos=0.113 W1=0.500  RF-domain AUC **0.648**  PO-risk R **0.000365**  τ(PO, RF)=0.487

构成（early vs late）：

| | n | pos | empty_any | sess0 | sku_clk |
|---|---:|---:|---:|---:|---:|
| early | 6282 | 0.063 | 0.581 | 0.990 | 0.0003 |
| late | 6282 | 0.163 | 0.401 | 0.979 | 0.0019 |

| tower | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| cross | 6 | 0.244 | 0.428 | -0.00003 | 0.000 |
| ctx | 6 | 0.457 | 0.284 | -0.00004 | 0.000 |
| user | 6 | 0.299 | 0.288 | -0.00000 | 0.000 |

PO-VIMP 序：

| rank | feat | tower | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `uc_n_clk_7d_log__x__dt_any_log` | cross | 0.2263 | 0.1823 |
| 2 | `dt_any_min` | ctx | 0.2172 | 0.1991 |
| 3 | `uc_lag_post_clk_1d_rate__x__dt_any_log` | cross | 0.1783 | 0.0365 |
| 4 | `n_prior_cnv` | user | 0.1027 | 0.0648 |
| 5 | `lag_empty_any` | user | 0.0739 | 0.0279 |
| 6 | `n_clk_before_7d` | user | 0.0403 | 0.1302 |
| 7 | `log1p_price` | ctx | 0.0352 | 0.1081 |
| 8 | `lag_post_clk_1d_rate` | user | 0.0328 | 0.0533 |
| 9 | `sess_pos` | ctx | 0.0276 | 0.0213 |
| 10 | `uc_n_clk_7d_log__x__sess_clk_before` | cross | 0.0199 | 0.0046 |
| 11 | `n_clk_before_1d` | user | 0.0198 | 0.0210 |
| 12 | `n_clk_before_1h` | user | 0.0186 | 0.0016 |
| 13 | `sess_clk_before` | ctx | 0.0037 | 0.0002 |
| 14 | `uc_lag_post_clk_1d_rate__x__sess_clk_before` | cross | 0.0037 | 0.0017 |
| 15 | `dt_any_min_miss` | ctx | 0.0000 | 0.0356 |
| 16 | `empty_any` | ctx | 0.0000 | 0.0932 |
| 17 | `uc_lag_post_clk_1d_rate__x__empty_any` | cross | 0.0000 | 0.0044 |
| 18 | `uc_lag_empty_any__x__empty_any` | cross | 0.0000 | 0.0141 |

用户时钟 `_seq_t_end` 和 150 板同一把刀（后来的人）。`cnv_ts` 是单的早晚，跟满窗 follow 缠在一起，只作对照。
SKU 漏斗构成两批都接近 0，不进 X。
