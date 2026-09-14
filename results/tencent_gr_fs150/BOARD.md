# TencentGR subset: 150 vs 128, then PO-risk

Subset: `seq/part-00000*` + `item_feat/` + `user_feat/` + `indexer.pkl`
(~1.9 GB). `mm_emb/` skipped (tens of GB; FS does not read it).

Users **6000**, label `future_cnv` (prefix feats → suffix conversion,
pos rate 0.157). Same combinatorial 1000-d generator as the 128 board;
F-score top **150** is a **strict superset** of the 128 (overlap 128,
22 new, 0 dropped).

## Holdout (HGB / LogReg)

| model | k | AUC | AP | Acc |
|---|---:|---:|---:|---:|
| HGB | 128 | 0.8214 | 0.4805 | 0.8533 |
| HGB | **150** | **0.8290** | **0.4914** | **0.8547** |
| LogReg | 128 | 0.8140 | 0.4887 | 0.8573 |
| LogReg | 150 | 0.8141 | 0.4906 | 0.8553 |

+22 extra F-score dims help the tree (~+0.008 AUC) and barely move linear.

## 22 new features (ranks 128–150, not in the 128 board)

`life_cnv_price_sum`, `x_life_ctcvr__7d_ctr`, `x_arpu_sum_proxy__sess_bounce_rate`,
`sess_depth_clk_mean`, `3d_n_clk`, `30d_clk_share`, `trend_cvr_3d_minus_14d`,
`x_arpu_mean_proxy__sess_n`, `x_hist_len__item_entropy_clk`, `trans_exp_to_clk`,
`trend_cvr_7d_minus_30d`, `1d_n_cnv`, `dec_hl7d_dec_arpu`, `x_7d_cvr__hist_len`,
`sq_item_entropy_clk`, `x_7d_cvr__sess_n`, `x_life_ctcvr__pay_cnt`,
`x_sess_n__item_entropy_clk`, `trans_clk_to_cnv`, `x_arpu_mean_proxy__active_days`,
`3d_clk_share`, `x_arpu_sum_proxy__arpu_mean_proxy`

Full numbered list: `FEATURES_150.txt`.

## Gated PO-risk (current method)

Clock = last prefix event time. 24 × 200 users, `future_cnv` Acc ↑, γ=1.5.

| method | Acc | fire |
|---|---:|---:|
| uniform last-two | 0.8298 | — |
| DRE last-two | 0.8223 | always |
| rfperm √PO | 0.8298 | **0.00** |
| resid drop-old | 0.8298 | 0.00 |

This clock is **quiet**: consecutive OOS does not jump, so rfperm stays
uniform (correct). Always-on DRE costs ~0.8 Acc points — same pattern as
covariate-stable real clocks.

HGB holdout Acc 0.855 vs online RF Acc 0.830 is the learner gap
(shallow RF probe vs HistGB), not a hop.
