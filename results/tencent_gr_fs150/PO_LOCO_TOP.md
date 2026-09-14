# PO-risk 哪维显著（LOCO，不是 impurity）

时钟 `_seq_t_end`。原始事件：`data/tencent_subset/seq/part-00000-*.parquet`。
用户矩阵：`results/tencent_gr_fs150/user_feats_selected.parquet`。

Impurity VIMP ≠ LOCO。中位切上 top impurity `life_cnv_share` 的 LOCO ΔR 是 **负的**。
Q1 vs Q4 才是 mix 变狠的时钟（RF-domain AUC 0.88）。

## 中位切 LOCO ΔR>0 的前几名

| feat | impurity | LOCO ΔR |
|---|---:|---:|
| dec_hl7d_dec_clk | 0.026 | **+9.4e-5** |
| x_life_cvr__7d_ctr | 0.026 | +6.1e-5 |
| attr_anyclk2cnv_min_std | 0.022 | +6.0e-5 |
| life_cnv_share | 0.045 | −4.4e-5 |

## Q1 vs Q4

| feat | impurity | LOCO ΔR |
|---|---:|---:|
| **dec_hl7d_dec_cnv** | **0.080** | **+3.6e-5** |
| dec_hl3d_dec_clk | 0.029 | +3.4e-5 |
| life_cnv_share | 0.003 | +2.9e-5 |

最显著：**7 天半衰期成交热度** `dec_hl7d_dec_cnv = Σ_{cnv} 2^{-age/7d}`。
Q1 均值 0.57 → Q4 0.04。后来的人最近没成交；lifetime `pay_cnt` 的 PO-mass≈0，扛不住这个差。
