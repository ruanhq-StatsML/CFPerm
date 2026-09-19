# FSDS multi-step iteration `iter04/_tmp_k10_s1_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=10
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→20 + π-stable MI→10 | fused→official_FSDS | 10 |  | 0.7625 | 0.0012 | 0.1308 | 0.5598 | 0.0011 | 0.4350 | 1.6685 | i_share_last,i_share_first,e_log1p_exp,i_n_users,u_n_events |

## Findings

