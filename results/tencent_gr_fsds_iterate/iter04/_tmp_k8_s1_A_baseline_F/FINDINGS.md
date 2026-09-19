# FSDS multi-step iteration `iter04/_tmp_k8_s1_A_baseline_F`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=8
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 8 | 1.0000 | 0.5815 | 0.0008 | 0.1438 | 0.5720 | 0.0012 | 0.3253 | 0.2565 | e_n_exp,u_n_exp,u_n_uniq_items,u_span_sec,i_n_users |

## Findings

- Best **W2** HGB AUC: `A_baseline_F|W1hold` = 0.5720 (baseline A=0.5720)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.5815
- Scope locked: graph features only; selection steps only; no graph algorithms.
