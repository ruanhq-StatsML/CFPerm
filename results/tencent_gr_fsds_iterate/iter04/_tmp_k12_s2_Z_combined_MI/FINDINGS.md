# FSDS multi-step iteration `iter04/_tmp_k12_s2_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=12
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→12 | fused→official_FSDS | 12 |  | 0.6092 | 0.0007 | 0.3593 | 0.8410 | 0.0085 | 0.5505 | 1.7550 | u_n_events,u_log1p_n_events,ui_pop_mismatch,u_n_uniq_items,i_share_first |

## Findings

