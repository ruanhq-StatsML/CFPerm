# FSDS multi-step iteration `iter02`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4300 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| C_cmean_then_F|W1hold | cmean|δ| prefilter→21 then F→15 | fused→official_FSDS_on_pool | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7664 | 0.0107 | 0.4669 | 0.3983 | u_span_sec,ui_pop_mismatch,i_n_exp,u_n_events,u_n_uniq_items |
| H_soft_corr|W1hold | F wide=21 → soft-corr@0.98 →13 | fused→official_FSDS_on_pool | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4184 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,u_log1p_n_uniq |
| I_stable_pi_F3|W1hold | 3-fold π-stable SelectKBest(F) (n_pos=3) | fused→official_FSDS_on_pool | 15 | 0.6667 | 0.4415 | 0.0006 | 0.9073 | 0.4871 | 0.0021 | 0.4609 | 0.3370 | u_n_uniq_items,i_n_covisit_neighbors,i_n_users,ui_pop_mismatch,u_log1p_n_events |
| J_delta_share_FSDS|W1hold | J* cumshare≥0.8 →15 cols then FSDS-F→15 | fused→official_FSDS_on_pool | 15 | 0.5789 | 0.6015 | 0.0009 | 0.9119 | 0.3960 | 0.0011 | 0.4527 | 0.3842 | u_span_sec,ui_pop_mismatch,i_n_covisit_neighbors,i_n_users,i_n_exp |
| F_cmean_stable|W1hold | cmean pre→21 + π-stable F→15 | fused→official_FSDS_on_pool | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4430 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| G_F_then_HGB_perm|W1hold | F screen=21 → HGB perm-imp →15 | fused→official_FSDS_on_pool | 15 | 0.6667 | 0.6434 | 0.0010 | 0.9641 | 0.7730 | 0.0064 | 0.6912 | 1.8905 | i_n_exp,i_share_linear,ui_pop_mismatch,u_n_exp,u_span_sec |

## Findings

- Best **W2** HGB AUC: `H_soft_corr|W1hold` = 0.7735 (baseline A=0.7679)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.6621
- cmean-|δ| prefilter aligns FS with period-shift guidance (same δ story as drill doc)
- Scope locked: graph features only; selection steps only; no graph algorithms.
