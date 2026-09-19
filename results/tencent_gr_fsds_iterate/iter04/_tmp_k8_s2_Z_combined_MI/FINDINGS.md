# FSDS multi-step iteration `iter04/_tmp_k8_s2_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=8
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→16 + π-stable MI→8 | fused→official_FSDS | 8 |  | 0.6400 | 0.0008 | 0.3413 | 0.7929 | 0.0050 | 0.4833 | 1.3420 | u_log1p_n_events,u_n_events,ui_pop_mismatch,u_n_uniq_items,i_credit_first |

## Findings

