# FSDS multi-step iteration `iter04/_tmp_k18_s0_A_baseline_F`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=18
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 18 | 1.0000 | 0.6441 | 0.0010 | 0.9204 | 0.7454 | 0.0058 | 0.4596 | 0.5559 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |

## Findings

- Best **W2** HGB AUC: `A_baseline_F|W1hold` = 0.7454 (baseline A=0.7454)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.6441
- Scope locked: graph features only; selection steps only; no graph algorithms.
