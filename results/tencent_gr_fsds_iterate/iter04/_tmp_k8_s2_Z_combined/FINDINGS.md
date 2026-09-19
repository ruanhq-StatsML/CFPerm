# FSDS multi-step iteration `iter04/_tmp_k8_s2_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=8
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→16 + π-stable F→8 | fused→official_FSDS | 8 |  | 0.7245 | 0.0009 | 0.6738 | 0.8191 | 0.0096 | 0.4102 | 0.3081 | i_n_exp,i_credit_first,i_log1p_n_exp,u_log1p_n_events,ui_pop_mismatch |

## Findings

