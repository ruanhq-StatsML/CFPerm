# FSDS multi-step overnight log

## OVERNIGHT_SUMMARY (iter01–10, wrap)

讲武德: **`W` = period**, PO-risk = `mean(τ̂²)` shift proxy — **not** an ATE.

### Locked DS recipe

```bash
PYTHONPATH=. python3 scripts/tencent_gr/po_help_fsds.py \
  --w1-grid results/tencent_gr_localize_fsds_time/W1_localized_grid.parquet \
  --w2-grid results/tencent_gr_localize_fsds_time/W2_localized_grid_sample.parquet \
  --select-k 15 --seed 0 --out-dir results/tencent_gr_fsds_iterate/po_help_cli
```

Or in Python: `po_help_select(...)` → hand `pool` to official FSDS.

### What worked

| recommendation | evidence |
|---|---|
| **PO-VIMP → FSDS** (`P_po_vimp` / CLI) | best mean W2 **~0.72** (@k=15); **0.724** @k=18 |
| **k=15** (try **18** for plain PO-VIMP) | k∈{10,12} hurts |
| **pool slack k+3** | tight pool=k collapses W2 |
| **rare-pos-capped π** for lower σ | σ **0.137** vs ~0.15–0.17 |
| **LOO-pos majority** to freeze one list | σ **0.114** @ mean 0.708 |
| always **seed mean±std** | seed gap ≫ method gap |

### What did not

| idea | result |
|---|---|
| hard corr@0.92 / J*-only | hurts W2 |
| α∈{0.3,0.5,0.7} tune | **identical** on d≈21 — keep 0.5 |
| bootstrap π alone | not free lunch |
| τ̂² row filter | ≈ plain P_po |
| seed-maj2 / avg-VIMP topk alone | hurts mean |

### Scope

图谱特征 only. No community / ego / GNN / ATE claims.

Cookbook: `docs/tencent_gr/PO_RISK_FOR_DS.md` · CLI: `scripts/tencent_gr/po_help_fsds.py`

---

## iter01 (2026-09-19 07:45 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_logreg_auc | W2_hgb_auc | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7679 | 0.4669 | 0.4064 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| C_cmean_then_F|W1hold | cmean|δ| prefilter→30 then F→15 | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7664 | 0.4669 | 0.4300 | u_span_sec,ui_pop_mismatch,i_n_exp,u_n_events,u_n_uniq_items |
| D_F_corr_prune|W1hold | F wide=21 → corr-prune@0.92 →10 | 10 | 0.4706 | 0.4924 | 0.8721 | 0.6437 | 0.7117 | 0.4502 | e_n_exp,u_n_exp,u_span_sec,u_log1p_n_uniq,i_n_covisit_neighbors |
| E_stable_pi_F|W1hold | 5-fold π-stable SelectKBest(F) | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7735 | 0.4669 | 0.4978 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| F_cmean_stable|W1hold | cmean pre→30 + π-stable F→15 | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7735 | 0.4669 | 0.4844 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| G_F_then_HGB_perm|W1hold | F screen=21 → HGB perm-imp →15 | 15 | 0.5789 | 0.6434 | 0.9468 | 0.7730 | 0.7042 | 2.1339 | i_n_exp,i_share_linear,ui_pop_mismatch,u_n_exp,u_span_sec |

See `iter01/FINDINGS.md`.

## iter02 (2026-09-19 07:54 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4300 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| C_cmean_then_F|W1hold | cmean|δ| prefilter→21 then F→15 | fused→official_FSDS_on_pool | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7664 | 0.0107 | 0.4669 | 0.3983 | u_span_sec,ui_pop_mismatch,i_n_exp,u_n_events,u_n_uniq_items |
| H_soft_corr|W1hold | F wide=21 → soft-corr@0.98 →13 | fused→official_FSDS_on_pool | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4184 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,u_log1p_n_uniq |
| I_stable_pi_F3|W1hold | 3-fold π-stable SelectKBest(F) (n_pos=3) | fused→official_FSDS_on_pool | 15 | 0.6667 | 0.4415 | 0.0006 | 0.9073 | 0.4871 | 0.0021 | 0.4609 | 0.3370 | u_n_uniq_items,i_n_covisit_neighbors,i_n_users,ui_pop_mismatch,u_log1p_n_events |
| J_delta_share_FSDS|W1hold | J* cumshare≥0.8 →15 cols then FSDS-F→15 | fused→official_FSDS_on_pool | 15 | 0.5789 | 0.6015 | 0.0009 | 0.9119 | 0.3960 | 0.0011 | 0.4527 | 0.3842 | u_span_sec,ui_pop_mismatch,i_n_covisit_neighbors,i_n_users,i_n_exp |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS_on_pool | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4430 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| G_F_then_HGB_perm|W1hold | F screen=21 → HGB perm-imp →15 | fused→official_FSDS_on_pool | 15 | 0.6667 | 0.6434 | 0.0010 | 0.9641 | 0.7730 | 0.0064 | 0.6912 | 1.8905 | i_n_exp,i_share_linear,ui_pop_mismatch,u_n_exp,u_span_sec |

See `iter02/FINDINGS.md`.

## iter03 (2026-09-19 07:58 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4146 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| C_cmean_then_F|W1hold | cmean|δ| prefilter→21 then F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7664 | 0.0107 | 0.4669 | 0.3815 | u_span_sec,ui_pop_mismatch,i_n_exp,u_n_events,u_n_uniq_items |
| H_soft_corr|W1hold | F wide=21 → soft-corr@0.98 →13 | fused→official_FSDS | 13 | 0.6471 | 0.7463 | 0.0010 | 0.8580 | 0.6821 | 0.0061 | 0.4168 | 0.3918 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,u_log1p_n_uniq |
| I_stable_pi_F3|W1hold | 3-fold π-stable SelectKBest(F) (n_pos=3) | fused→official_FSDS | 15 | 0.6667 | 0.4415 | 0.0006 | 0.9073 | 0.4871 | 0.0021 | 0.4609 | 0.3296 | u_n_uniq_items,i_n_covisit_neighbors,i_n_users,ui_pop_mismatch,u_log1p_n_events |
| J_delta_share_FSDS|W1hold | J* cumshare≥0.8 →15 cols then FSDS-F→15 | fused→official_FSDS | 15 | 0.5789 | 0.6015 | 0.0009 | 0.9119 | 0.3960 | 0.0011 | 0.4527 | 0.4300 | u_span_sec,ui_pop_mismatch,i_n_covisit_neighbors,i_n_users,i_n_exp |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4116 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| G_F_then_HGB_perm|W1hold | F screen=21 → HGB perm-imp →15 | fused→official_FSDS | 15 | 0.5789 | 0.6434 | 0.0010 | 0.9468 | 0.7730 | 0.0064 | 0.7042 | 1.8236 | i_n_exp,i_share_linear,ui_pop_mismatch,u_n_exp,u_span_sec |

See `iter03/FINDINGS.md`.

## iter03 (2026-09-19 07:58 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4245 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| H_soft_corr|W1hold | F wide=21 → soft-corr@0.98 →13 | fused→official_FSDS | 13 | 0.6471 | 0.7463 | 0.0010 | 0.8580 | 0.6821 | 0.0061 | 0.4168 | 0.4169 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,u_log1p_n_uniq |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4013 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| G_F_then_HGB_perm|W1hold | F screen=21 → HGB perm-imp →15 | fused→official_FSDS | 15 | 0.5789 | 0.6434 | 0.0010 | 0.9468 | 0.7730 | 0.0064 | 0.7042 | 1.8476 | i_n_exp,i_share_linear,ui_pop_mismatch,u_n_exp,u_span_sec |
| Z_combined|W1hold | COMBINED: cmean-soft→21 | F-wide→21 | π3 | soft-corr@0.98→13 | k=15 | fused→official_FSDS | 13 | 0.6471 | 0.7463 | 0.0010 | 0.8566 | 0.6968 | 0.0074 | 0.4161 | 0.4037 | e_log1p_exp,u_span_sec,i_share_first,i_log1p_n_exp,u_n_uniq_items |

See `iter03/FINDINGS.md`.

## iter03 (2026-09-19 07:58 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4135 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.3908 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| Z_combined|W1hold | COMBINED: cmean-soft→21 | F-wide→21 | π3 | soft-corr@0.98→13 | k=15 | fused→official_FSDS | 13 | 0.6471 | 0.7463 | 0.0010 | 0.8566 | 0.6968 | 0.0074 | 0.4161 | 0.4139 | e_log1p_exp,u_span_sec,i_share_first,i_log1p_n_exp,u_n_uniq_items |

See `iter03/FINDINGS.md`.

## iter03 (2026-09-19 07:59 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4373 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4225 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| Z_combined|W1hold | COMBINED: cmean-soft→21 | π3→13 | twin-dedupe@0.999 | k=15 | fused→official_FSDS | 13 | 0.5556 | 0.4415 | 0.0006 | 0.8682 | 0.4905 | 0.0021 | 0.4463 | 0.3509 | e_log1p_exp,e_n_exp,i_share_first,i_log1p_n_exp,u_n_uniq_items |

See `iter03/FINDINGS.md`.

## iter03 (2026-09-19 07:59 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4216 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4033 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| Z_combined|W1hold | COMBINED: cmean-soft→21 | π3→13 | twin-dedupe@0.999 | k=15 | fused→official_FSDS | 13 | 0.5556 | 0.4415 | 0.0006 | 0.8682 | 0.4905 | 0.0021 | 0.4463 | 0.3305 | e_log1p_exp,e_n_exp,i_share_first,i_log1p_n_exp,u_n_uniq_items |
| H_soft_corr|W1hold | F wide=21 → soft-corr@0.98 →13 | fused→official_FSDS | 13 | 0.6471 | 0.7463 | 0.0010 | 0.8580 | 0.6821 | 0.0061 | 0.4168 | 0.4074 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,u_log1p_n_uniq |

See `iter03/FINDINGS.md`.

## iter03 (2026-09-19 07:59 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4446 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4048 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.3930 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |

See `iter03/FINDINGS.md`.

## iter04/_tmp_k8_s0_A_baseline_F (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 8 | 1.0000 | 0.4748 | 0.0006 | 0.2533 | 0.3268 | 0.0018 | 0.2558 | 0.2096 | e_n_exp,u_n_uniq_items,u_span_sec,i_n_exp,i_log1p_n_exp |

See `iter04/_tmp_k8_s0_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k8_s0_Z_combined (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→16 + π-stable F→8 | fused→official_FSDS | 8 |  | 0.5775 | 0.0008 | 0.6050 | 0.2894 | 0.0008 | 0.2917 | 0.2648 | i_log1p_n_exp,u_n_uniq_items,i_n_exp,i_log1p_n_covisit,u_n_exp |

See `iter04/_tmp_k8_s0_Z_combined/FINDINGS.md`.

## iter04/_tmp_k8_s0_Z_combined_MI (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→16 + π-stable MI→8 | fused→official_FSDS | 8 |  | 0.4190 | 0.0006 | 0.6022 | 0.5969 | 0.0015 | 0.4512 | 1.4762 | u_log1p_n_uniq,i_credit_last,i_log1p_n_exp,ui_pop_mismatch,u_n_events |

See `iter04/_tmp_k8_s0_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k8_s1_A_baseline_F (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 8 | 1.0000 | 0.5815 | 0.0008 | 0.1438 | 0.5720 | 0.0012 | 0.3253 | 0.2565 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,i_n_users |

See `iter04/_tmp_k8_s1_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k8_s1_Z_combined (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→16 + π-stable F→8 | fused→official_FSDS | 8 |  | 0.5921 | 0.0007 | 0.1152 | 0.7394 | 0.0039 | 0.3235 | 0.3276 | u_span_sec,i_n_users,i_credit_last,i_n_covisit_neighbors,u_n_exp |

See `iter04/_tmp_k8_s1_Z_combined/FINDINGS.md`.

## iter04/_tmp_k8_s1_Z_combined_MI (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→16 + π-stable MI→8 | fused→official_FSDS | 8 |  | 0.6878 | 0.0010 | 0.3472 | 0.4159 | 0.0010 | 0.3091 | 1.3433 | u_n_exp,i_log1p_n_users,u_log1p_n_uniq,i_log1p_n_exp,i_credit_last |

See `iter04/_tmp_k8_s1_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k8_s2_A_baseline_F (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 8 | 1.0000 | 0.5708 | 0.0007 | 0.0754 | 0.8224 | 0.0086 | 0.3656 | 0.1953 | e_n_exp,u_log1p_n_events,i_n_exp,i_log1p_n_exp,e_log1p_exp |

See `iter04/_tmp_k8_s2_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k8_s2_Z_combined (2026-09-19 08:02 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→16 + π-stable F→8 | fused→official_FSDS | 8 |  | 0.7245 | 0.0009 | 0.6738 | 0.8191 | 0.0096 | 0.4102 | 0.3081 | i_n_exp,i_credit_first,i_log1p_n_exp,u_log1p_n_events,ui_pop_mismatch |

See `iter04/_tmp_k8_s2_Z_combined/FINDINGS.md`.

## iter04/_tmp_k8_s2_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→16 + π-stable MI→8 | fused→official_FSDS | 8 |  | 0.6400 | 0.0008 | 0.3413 | 0.7929 | 0.0050 | 0.4833 | 1.3420 | u_log1p_n_events,u_n_events,ui_pop_mismatch,u_n_uniq_items,i_credit_first |

See `iter04/_tmp_k8_s2_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k10_s0_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 10 | 1.0000 | 0.5638 | 0.0008 | 0.3650 | 0.3687 | 0.0010 | 0.2614 | 0.3230 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,i_n_exp |

See `iter04/_tmp_k10_s0_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k10_s0_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→20 + π-stable F→10 | fused→official_FSDS | 10 |  | 0.5632 | 0.0008 | 0.3650 | 0.3687 | 0.0010 | 0.2614 | 0.2761 | e_log1p_exp,e_n_exp,i_log1p_n_exp,i_share_first,u_n_uniq_items |

See `iter04/_tmp_k10_s0_Z_combined/FINDINGS.md`.

## iter04/_tmp_k10_s0_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→20 + π-stable MI→10 | fused→official_FSDS | 10 |  | 0.6542 | 0.0010 | 0.9521 | 0.2804 | 0.0016 | 0.7524 | 1.7060 | i_share_last,i_credit_linear,i_share_first,i_credit_last,u_n_exp |

See `iter04/_tmp_k10_s0_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k10_s1_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 10 | 1.0000 | 0.5572 | 0.0007 | 0.1149 | 0.6133 | 0.0013 | 0.3616 | 0.2656 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,i_n_users |

See `iter04/_tmp_k10_s1_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k10_s1_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→20 + π-stable F→10 | fused→official_FSDS | 10 |  | 0.5572 | 0.0007 | 0.1149 | 0.6060 | 0.0013 | 0.3616 | 0.2787 | u_span_sec,i_n_users,i_credit_last,i_share_last,e_log1p_exp |

See `iter04/_tmp_k10_s1_Z_combined/FINDINGS.md`.

## iter04/_tmp_k10_s1_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→20 + π-stable MI→10 | fused→official_FSDS | 10 |  | 0.7625 | 0.0012 | 0.1308 | 0.5598 | 0.0011 | 0.4350 | 1.6685 | i_share_last,i_share_first,e_log1p_exp,i_n_users,u_n_events |

See `iter04/_tmp_k10_s1_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k10_s2_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 10 | 1.0000 | 0.5879 | 0.0007 | 0.3893 | 0.6958 | 0.0032 | 0.3617 | 0.1957 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,i_n_exp |

See `iter04/_tmp_k10_s2_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k10_s2_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→20 + π-stable F→10 | fused→official_FSDS | 10 |  | 0.4529 | 0.0006 | 0.3258 | 0.6956 | 0.0025 | 0.3664 | 0.2805 | i_n_exp,e_log1p_exp,e_n_exp,i_share_first,i_credit_first |

See `iter04/_tmp_k10_s2_Z_combined/FINDINGS.md`.

## iter04/_tmp_k10_s2_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→20 + π-stable MI→10 | fused→official_FSDS | 10 |  | 0.6400 | 0.0008 | 0.3159 | 0.7929 | 0.0050 | 0.4823 | 1.6240 | u_log1p_n_events,u_n_events,ui_pop_mismatch,i_share_first,u_log1p_n_uniq |

See `iter04/_tmp_k10_s2_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k12_s0_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 12 | 1.0000 | 0.5909 | 0.0009 | 0.8252 | 0.3936 | 0.0016 | 0.3730 | 0.2871 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |

See `iter04/_tmp_k12_s0_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k12_s0_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→12 | fused→official_FSDS | 12 |  | 0.5913 | 0.0009 | 0.8252 | 0.3938 | 0.0016 | 0.3730 | 0.2889 | e_log1p_exp,e_n_exp,i_log1p_n_exp,i_share_first,i_credit_first |

See `iter04/_tmp_k12_s0_Z_combined/FINDINGS.md`.

## iter04/_tmp_k12_s0_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→12 | fused→official_FSDS | 12 |  | 0.5055 | 0.0007 | 0.6293 | 0.3465 | 0.0019 | 0.7972 | 1.8594 | e_log1p_exp,e_n_exp,i_share_last,i_share_linear,i_share_first |

See `iter04/_tmp_k12_s0_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k12_s1_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 12 | 1.0000 | 0.5921 | 0.0007 | 0.1397 | 0.7429 | 0.0039 | 0.4089 | 0.3991 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,i_n_users |

See `iter04/_tmp_k12_s1_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k12_s1_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→12 | fused→official_FSDS | 12 |  | 0.5921 | 0.0007 | 0.1397 | 0.7389 | 0.0039 | 0.4089 | 0.4094 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |

See `iter04/_tmp_k12_s1_Z_combined/FINDINGS.md`.

## iter04/_tmp_k12_s1_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→12 | fused→official_FSDS | 12 |  | 0.4854 | 0.0007 | 0.1401 | 0.5919 | 0.0013 | 0.3628 | 1.8095 | i_share_first,i_credit_last,i_log1p_n_users,u_log1p_n_uniq,u_n_uniq_items |

See `iter04/_tmp_k12_s1_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k12_s2_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 12 | 1.0000 | 0.7236 | 0.0009 | 0.5482 | 0.8068 | 0.0100 | 0.6758 | 0.3072 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,i_n_exp |

See `iter04/_tmp_k12_s2_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k12_s2_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→12 | fused→official_FSDS | 12 |  | 0.6490 | 0.0008 | 0.5457 | 0.7786 | 0.0081 | 0.6872 | 0.4064 | i_n_exp,e_log1p_exp,e_n_exp,i_share_first,i_credit_first |

See `iter04/_tmp_k12_s2_Z_combined/FINDINGS.md`.

## iter04/_tmp_k12_s2_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→12 | fused→official_FSDS | 12 |  | 0.6092 | 0.0007 | 0.3593 | 0.8410 | 0.0085 | 0.5505 | 1.7550 | u_n_events,u_log1p_n_events,ui_pop_mismatch,u_n_uniq_items,i_share_first |

See `iter04/_tmp_k12_s2_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k15_s0_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4263 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |

See `iter04/_tmp_k15_s0_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k15_s0_Z_combined (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 |  | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4150 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |

See `iter04/_tmp_k15_s0_Z_combined/FINDINGS.md`.

## iter04/_tmp_k15_s0_Z_combined_MI (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→15 | fused→official_FSDS | 15 |  | 0.5802 | 0.0008 | 0.8020 | 0.5402 | 0.0015 | 0.7980 | 1.8975 | e_log1p_exp,e_n_exp,i_share_last,i_credit_linear,i_share_first |

See `iter04/_tmp_k15_s0_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k15_s1_A_baseline_F (2026-09-19 08:03 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4540 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |

See `iter04/_tmp_k15_s1_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k15_s1_Z_combined (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 |  | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4690 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |

See `iter04/_tmp_k15_s1_Z_combined/FINDINGS.md`.

## iter04/_tmp_k15_s1_Z_combined_MI (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→15 | fused→official_FSDS | 15 |  | 0.6250 | 0.0010 | 0.1552 | 0.4158 | 0.0010 | 0.4164 | 1.7833 | e_log1p_exp,i_share_first,i_log1p_n_users,u_log1p_n_uniq,i_share_linear |

See `iter04/_tmp_k15_s1_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k15_s2_A_baseline_F (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3121 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |

See `iter04/_tmp_k15_s2_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k15_s2_Z_combined (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 |  | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 0.4461 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |

See `iter04/_tmp_k15_s2_Z_combined/FINDINGS.md`.

## iter04/_tmp_k15_s2_Z_combined_MI (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→15 | fused→official_FSDS | 15 |  | 0.4682 | 0.0005 | 0.4229 | 0.8192 | 0.0097 | 0.4850 | 1.8228 | u_n_events,u_log1p_n_events,ui_pop_mismatch,e_log1p_exp,e_n_exp |

See `iter04/_tmp_k15_s2_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k18_s0_A_baseline_F (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 18 | 1.0000 | 0.6441 | 0.0010 | 0.9204 | 0.7454 | 0.0058 | 0.4596 | 0.5559 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |

See `iter04/_tmp_k18_s0_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k18_s0_Z_combined (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→18 | fused→official_FSDS | 18 |  | 0.6441 | 0.0010 | 0.9204 | 0.7454 | 0.0058 | 0.4596 | 0.5408 | e_log1p_exp,e_n_exp,u_n_uniq_items,i_n_exp,i_log1p_n_covisit |

See `iter04/_tmp_k18_s0_Z_combined/FINDINGS.md`.

## iter04/_tmp_k18_s0_Z_combined_MI (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→18 | fused→official_FSDS | 18 |  | 0.5981 | 0.0009 | 0.9144 | 0.7465 | 0.0031 | 0.8223 | 1.8947 | e_log1p_exp,e_n_exp,i_credit_last,i_share_last,i_log1p_n_covisit |

See `iter04/_tmp_k18_s0_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k18_s1_A_baseline_F (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 18 | 1.0000 | 0.7666 | 0.0011 | 0.1738 | 0.5285 | 0.0012 | 0.3323 | 0.4950 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |

See `iter04/_tmp_k18_s1_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k18_s1_Z_combined (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→18 | fused→official_FSDS | 18 |  | 0.7666 | 0.0011 | 0.1738 | 0.5286 | 0.0012 | 0.3323 | 0.5287 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |

See `iter04/_tmp_k18_s1_Z_combined/FINDINGS.md`.

## iter04/_tmp_k18_s1_Z_combined_MI (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→18 | fused→official_FSDS | 18 |  | 0.6380 | 0.0010 | 0.1790 | 0.5819 | 0.0020 | 0.3982 | 1.9823 | i_credit_last,e_log1p_exp,i_share_last,e_n_exp,i_share_first |

See `iter04/_tmp_k18_s1_Z_combined_MI/FINDINGS.md`.

## iter04/_tmp_k18_s2_A_baseline_F (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 18 | 1.0000 | 0.6927 | 0.0010 | 0.5012 | 0.8672 | 0.0131 | 0.6612 | 0.4376 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,u_log1p_n_events |

See `iter04/_tmp_k18_s2_A_baseline_F/FINDINGS.md`.

## iter04/_tmp_k18_s2_Z_combined (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→18 | fused→official_FSDS | 18 |  | 0.6929 | 0.0010 | 0.5556 | 0.8671 | 0.0131 | 0.7212 | 0.4480 | e_log1p_exp,e_n_exp,i_share_first,i_credit_first,i_n_exp |

See `iter04/_tmp_k18_s2_Z_combined/FINDINGS.md`.

## iter04/_tmp_k18_s2_Z_combined_MI (2026-09-19 08:04 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→18 | fused→official_FSDS | 18 |  | 0.6160 | 0.0007 | 0.5372 | 0.8304 | 0.0176 | 0.7541 | 1.9056 | u_n_events,u_log1p_n_events,ui_pop_mismatch,e_log1p_exp,e_n_exp |

See `iter04/_tmp_k18_s2_Z_combined_MI/FINDINGS.md`.

## iter04 (2026-09-19 08:05 UTC)

- k×seed×{A,Z,Z-MI} sweep (45 runs)
- Seed variance dwarfs method gap; Z≈A at k=15 mean; MI hurts; k≈15–18 sweet.
- See iter04/FINDINGS.md


## iter05/s0 (2026-09-19 08:22 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4278 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4272 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→21 (risk=1.1919e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7664 | 0.0107 | 0.4669 | 2.7163 | ui_pop_mismatch,i_log1p_n_users,i_log1p_n_covisit,u_n_exp,u_n_uniq_items |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→21 | 5-fold π-stable SelectKBest(F) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 2.7288 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |

See `iter05/s0/FINDINGS.md`.

## iter05/s1 (2026-09-19 08:22 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4525 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4938 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→21 (risk=1.25052e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5451 | 0.0012 | 0.5088 | 2.6935 | i_n_users,i_log1p_n_users,i_n_covisit_neighbors,i_log1p_n_covisit,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→21 | 5-fold π-stable SelectKBest(F) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 2.7631 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |

See `iter05/s1/FINDINGS.md`.

## iter05/s2 (2026-09-19 08:22 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3329 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 0.4256 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→21 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 1.0000 | 0.5464 | 0.0007 | 0.4846 | 0.8412 | 0.0171 | 0.6457 | 2.5782 | ui_pop_mismatch,i_log1p_n_exp,u_n_exp,i_n_exp,u_n_uniq_items |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→21 | 5-fold π-stable SelectKBest(F) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 2.7340 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |

See `iter05/s2/FINDINGS.md`.

## iter05/tight_s0 (2026-09-19 08:23 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4209 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4184 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.1919e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.5819 | 0.0008 | 0.9274 | 0.7593 | 0.0061 | 0.4460 | 2.7891 | ui_pop_mismatch,i_log1p_n_users,i_log1p_n_covisit,i_n_covisit_neighbors,u_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.5819 | 0.0008 | 0.9274 | 0.7388 | 0.0057 | 0.4460 | 2.8002 | u_n_uniq_items,i_n_exp,i_log1p_n_users,i_n_covisit_neighbors,ui_pop_mismatch |

See `iter05/tight_s0/FINDINGS.md`.

## iter05/tight_s1 (2026-09-19 08:23 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4516 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4437 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.25052e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.7666 | 0.0011 | 0.2027 | 0.5362 | 0.0012 | 0.3296 | 2.7806 | ui_pop_mismatch,i_n_users,i_log1p_n_users,i_log1p_n_exp,i_n_covisit_neighbors |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7666 | 0.0011 | 0.2027 | 0.5286 | 0.0012 | 0.3296 | 2.8004 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |

See `iter05/tight_s1/FINDINGS.md`.

## iter05/tight_s2 (2026-09-19 08:23 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3182 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 0.4236 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6922 | 0.0010 | 0.5831 | 0.8710 | 0.0136 | 0.5153 | 2.7186 | ui_pop_mismatch,i_n_users,i_log1p_n_exp,u_n_exp,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.6929 | 0.0010 | 0.6322 | 0.8671 | 0.0131 | 0.6211 | 2.6992 | i_share_first,i_credit_first,i_n_exp,u_span_sec,i_log1p_n_covisit |

See `iter05/tight_s2/FINDINGS.md`.

## iter05 (2026-09-19 08:23 UTC)

- PO-risk helper + P_po_vimp_FSDS / Z_combined_PO
- PO-VIMP→FSDS mean W2 0.722 (best); ui_pop_mismatch dominates VIMP
- See iter05/FINDINGS.md


## iter06 (2026-09-19 08:25 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4262 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4323 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.1919e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.5819 | 0.0008 | 0.9274 | 0.7593 | 0.0061 | 0.4460 | 2.7986 | ui_pop_mismatch,i_log1p_n_users,i_log1p_n_covisit,i_n_covisit_neighbors,u_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.5819 | 0.0008 | 0.9274 | 0.7388 | 0.0057 | 0.4460 | 2.8451 | u_n_uniq_items,i_n_exp,i_log1p_n_users,i_n_covisit_neighbors,ui_pop_mismatch |
| P_po_boot_pi|W1hold | PO-boot-π: PO-help pool→18 | bootπ×40→15 | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.7253 | 0.0009 | 0.8337 | 0.6682 | 0.0050 | 0.4118 | 2.9179 | u_span_sec,i_log1p_n_users,ui_pop_mismatch,i_n_covisit_neighbors,i_credit_last |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.5932 | 0.0008 | 0.9320 | 0.7371 | 0.0056 | 0.5905 | 2.8116 | i_log1p_n_exp,u_n_uniq_items,i_log1p_n_users,i_n_covisit_neighbors,i_n_users |

See `iter06/FINDINGS.md`.

## iter06_s1 (2026-09-19 08:26 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4555 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4415 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.25052e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.7666 | 0.0011 | 0.2027 | 0.5362 | 0.0012 | 0.3296 | 2.7464 | ui_pop_mismatch,i_n_users,i_log1p_n_users,i_log1p_n_exp,i_n_covisit_neighbors |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7666 | 0.0011 | 0.2027 | 0.5286 | 0.0012 | 0.3296 | 2.7434 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |
| P_po_boot_pi|W1hold | PO-boot-π: PO-help pool→18 | bootπ×40→15 | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.5228 | 0.0008 | 0.1923 | 0.5649 | 0.0017 | 0.3368 | 2.8012 | i_n_users,i_n_covisit_neighbors,i_share_last,i_log1p_n_users,i_share_linear |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7310 | 0.0011 | 0.1938 | 0.5713 | 0.0014 | 0.3310 | 2.7395 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |

See `iter06_s1/FINDINGS.md`.

## iter06_s2 (2026-09-19 08:26 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3257 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 0.4337 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6922 | 0.0010 | 0.5831 | 0.8710 | 0.0136 | 0.5153 | 2.6927 | ui_pop_mismatch,i_n_users,i_log1p_n_exp,u_n_exp,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.6929 | 0.0010 | 0.6322 | 0.8671 | 0.0131 | 0.6211 | 2.7049 | i_share_first,i_credit_first,i_n_exp,u_span_sec,i_log1p_n_covisit |
| P_po_boot_pi|W1hold | PO-boot-π: PO-help pool→18 | bootπ×40→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6922 | 0.0010 | 0.7010 | 0.8707 | 0.0132 | 0.8291 | 2.8500 | i_share_first,i_credit_first,i_credit_last,ui_pop_mismatch,u_log1p_n_events |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7843 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |

See `iter06_s2/FINDINGS.md`.

## iter06 (PO-help + rare-pos / bootstrap π)

| variant | mean_W2 | std_W2 | mean_AP |
|---|---:|---:|---:|
| P_po_vimp_FSDS | 0.7221 | 0.1705 | 0.0070 |
| Z_combined | 0.7197 | 0.1572 | 0.0072 |
| Z_combined_PO_rare | 0.7171 | 0.1369 | 0.0054 |
| A_baseline_F | 0.7156 | 0.1537 | 0.0099 |
| Z_combined_PO | 0.7115 | 0.1709 | 0.0067 |
| P_po_boot_pi | 0.7013 | 0.1555 | 0.0066 |

See `iter06/FINDINGS.md`.

## iter07 (2026-09-19 08:53 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4122 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.1919e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.5819 | 0.0008 | 0.9274 | 0.7593 | 0.0061 | 0.4460 | 2.7870 | ui_pop_mismatch,i_log1p_n_users,i_log1p_n_covisit,i_n_covisit_neighbors,u_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO(α=0.5): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.5932 | 0.0008 | 0.9320 | 0.7371 | 0.0056 | 0.5905 | 2.8289 | i_log1p_n_exp,u_n_uniq_items,i_log1p_n_users,i_n_covisit_neighbors,i_n_users |
| Z_combined_PO_a03|W1hold | COMBINED-PO(α=0.3): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.5932 | 0.0008 | 0.9320 | 0.7371 | 0.0056 | 0.5905 | 2.8217 | i_log1p_n_exp,u_n_uniq_items,i_log1p_n_users,i_n_covisit_neighbors,i_n_users |
| Z_combined_PO_a07|W1hold | COMBINED-PO(α=0.7): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.5932 | 0.0008 | 0.9320 | 0.7371 | 0.0056 | 0.5905 | 2.7915 | i_log1p_n_exp,u_n_uniq_items,i_log1p_n_users,i_n_covisit_neighbors,i_n_users |
| P_po_vimp_rare_pi|W1hold | PO-VIMP+rare-π: pre→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.5932 | 0.0008 | 0.9320 | 0.7371 | 0.0056 | 0.5905 | 2.7951 | i_log1p_n_exp,u_n_uniq_items,i_log1p_n_users,i_n_covisit_neighbors,i_n_users |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.5932 | 0.0008 | 0.9320 | 0.7371 | 0.0056 | 0.5905 | 2.8164 | i_log1p_n_exp,u_n_uniq_items,i_log1p_n_users,i_n_covisit_neighbors,i_n_users |

See `iter07/FINDINGS.md`.

## iter07_s1 (2026-09-19 08:53 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4516 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.25052e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.7666 | 0.0011 | 0.2027 | 0.5362 | 0.0012 | 0.3296 | 2.7286 | ui_pop_mismatch,i_n_users,i_log1p_n_users,i_log1p_n_exp,i_n_covisit_neighbors |
| Z_combined_PO|W1hold | COMBINED-PO(α=0.5): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7310 | 0.0011 | 0.1938 | 0.5713 | 0.0014 | 0.3310 | 2.7426 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |
| Z_combined_PO_a03|W1hold | COMBINED-PO(α=0.3): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7310 | 0.0011 | 0.1938 | 0.5713 | 0.0014 | 0.3310 | 2.7647 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |
| Z_combined_PO_a07|W1hold | COMBINED-PO(α=0.7): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7310 | 0.0011 | 0.1938 | 0.5713 | 0.0014 | 0.3310 | 2.7708 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |
| P_po_vimp_rare_pi|W1hold | PO-VIMP+rare-π: pre→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7310 | 0.0011 | 0.1938 | 0.5713 | 0.0014 | 0.3310 | 2.7465 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.7310 | 0.0011 | 0.1938 | 0.5713 | 0.0014 | 0.3310 | 2.7376 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |

See `iter07_s1/FINDINGS.md`.

## iter07_s2 (2026-09-19 08:54 UTC)

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3235 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6922 | 0.0010 | 0.5831 | 0.8710 | 0.0136 | 0.5153 | 2.6972 | ui_pop_mismatch,i_n_users,i_log1p_n_exp,u_n_exp,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO(α=0.5): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7524 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| Z_combined_PO_a03|W1hold | COMBINED-PO(α=0.3): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7461 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| Z_combined_PO_a07|W1hold | COMBINED-PO(α=0.7): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7538 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| P_po_vimp_rare_pi|W1hold | PO-VIMP+rare-π: pre→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7437 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7643 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |

See `iter07_s2/FINDINGS.md`.

## iter07 (α ablate + PO×rare-π)

| variant | mean_W2 | std_W2 | mean_AP |
|---|---:|---:|---:|
| P_po_vimp_FSDS | 0.7221 | 0.1705 | 0.0070 |
| P_po_vimp_rare_pi | 0.7171 | 0.1369 | 0.0054 |
| Z_combined_PO | 0.7171 | 0.1369 | 0.0054 |
| Z_combined_PO_a07 | 0.7171 | 0.1369 | 0.0054 |
| Z_combined_PO_a03 | 0.7171 | 0.1369 | 0.0054 |
| Z_combined_PO_rare | 0.7171 | 0.1369 | 0.0054 |
| A_baseline_F | 0.7156 | 0.1537 | 0.0099 |

α∈{0.3,0.5,0.7} identical; tight pool k hurts (see FINDINGS).

See `iter07/FINDINGS.md`.

## iter08 (k-sweep + τ̂² row filter)

Best mean W2: P_po_vimp @k=18 → 0.724; Z_PO_rare @k=15 σ=0.137; k<15 hurts; τ̂²-rows ≈ P_po.

See `iter08/FINDINGS.md`.

Summary table:

```
variant,k,mean_W2,std_W2,mean_AP
A_baseline_F,10,0.5592465530546026,0.1701061746924145,0.001821965186553154
A_baseline_F,12,0.647781682611404,0.22239353791549465,0.005165173780251741
A_baseline_F,15,0.7155921401277437,0.1536922211137169,0.00987705041974968
A_baseline_F,18,0.7136978309578929,0.17153121604245025,0.006701876668429305
P_po_tau2_rows,15,0.7139906934798576,0.16005766132077662,0.006230654223196792
P_po_vimp_FSDS,10,0.4719107652544185,0.03135230886729512,0.0010582721420647043
P_po_vimp_FSDS,12,0.5527803345140807,0.18659196473900266,0.0016493015100453713
P_po_vimp_FSDS,15,0.7221455201331363,0.17046694918816666,0.006961923821458033
P_po_vimp_FSDS,18,0.7241978820926188,0.17117319703240805,0.006989842155813718
Z_combined_PO_rare,10,0.4237790422000948,0.10449305478338501,0.0010826452941216356
Z_combined_PO_rare,12,0.5623227716880969,0.11810028096087313,0.0018119454548845499
Z_combined_PO_rare,15,0.7170831822534609,0.13692106131117482,0.00535008085152226
Z_combined_PO_rare,18,0.7131283760540725,0.16416980871708853,0.007185296946281637
```

## iter09 (seed maj / LOO-pos stability)

Single-seed P_po still best mean (0.722). LOO-pos maj lowest σ (0.114) @0.708. Seed-maj2 hurts. Avg-VIMP topk alone hurts.

See `iter09/FINDINGS.md`.

```
variant,mean_W2,std_W2,mean_AP,n_feat
P_po_vimp_FSDS__single,0.7221455201331363,0.17046694918816666,0.006961923821458033,15
Z_combined_PO_rare__single,0.7170831822534609,0.13692106131117482,0.00535008085152226,15
A_baseline_F__single,0.7155921401277437,0.1536922211137169,0.00987705041974968,15
P_po__LOO_pos_maj,0.7077720134995676,0.11435653260939743,0.003253813678269336,15
Z_combined_PO_rare__intersect_fill,0.6875598508725442,0.19196005572788702,0.0057125334313500135,15
Z_combined_PO_rare__maj2,0.683739854405489,0.18714320629559722,0.005659996990687779,15
P_po_vimp_FSDS__intersect_fill,0.6776803893676959,0.15799646692809496,0.005475077246394265,15
P_po_vimp_FSDS__maj2,0.6740986342379531,0.1522593208009497,0.006098654265749095,15
A_baseline_F__maj2,0.5564608726373432,0.19932830102617202,0.003894810696854073,15
A_baseline_F__intersect_fill,0.5563353601279297,0.19911327075408525,0.003894810696854073,15
P_po_avgVIMP_topk,0.5508395392296321,0.18712717364432277,0.001665226715585677,15
P_po_vimp__seed_avg_topk,0.5508395392296321,0.18712717364432277,0.001665226715585677,15
P_po_avgVIMP_pool__maj2,0.42064355377048873,0.15709624189125035,0.0015291868994964056,15
```

## iter10 (DS handoff wrap)

CLI smoke mean W2≈0.720 @k=15. See `iter10/FINDINGS.md` + OVERNIGHT_SUMMARY above.
