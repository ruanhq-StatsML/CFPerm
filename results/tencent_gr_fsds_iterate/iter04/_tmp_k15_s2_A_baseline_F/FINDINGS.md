# FSDS multi-step iteration `iter04/_tmp_k15_s2_A_baseline_F`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3121 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |

## Findings

- Best **W2** HGB AUC: `A_baseline_F|W1hold` = 0.8363 (baseline A=0.8363)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.5459
- Scope locked: graph features only; selection steps only; no graph algorithms.
