# FSDS multi-step iteration `iter04/_tmp_k8_s0_A_baseline_F`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=8
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 8 | 1.0000 | 0.4748 | 0.0006 | 0.2533 | 0.3268 | 0.0018 | 0.2558 | 0.2096 | e_n_exp,u_n_uniq_items,u_span_sec,i_n_exp,i_log1p_n_exp |

## Findings

- Best **W2** HGB AUC: `A_baseline_F|W1hold` = 0.3268 (baseline A=0.3268)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.4748
- Scope locked: graph features only; selection steps only; no graph algorithms.
