# FSDS multi-step iteration `iter04/_tmp_k15_s2_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→15 | fused→official_FSDS | 15 |  | 0.4682 | 0.0005 | 0.4229 | 0.8192 | 0.0097 | 0.4850 | 1.8228 | u_n_events,u_log1p_n_events,ui_pop_mismatch,e_log1p_exp,e_n_exp |

## Findings

