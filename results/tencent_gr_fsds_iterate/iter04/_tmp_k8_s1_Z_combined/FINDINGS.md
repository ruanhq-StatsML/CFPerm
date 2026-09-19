# FSDS multi-step iteration `iter04/_tmp_k8_s1_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=8
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→16 + π-stable F→8 | fused→official_FSDS | 8 |  | 0.5921 | 0.0007 | 0.1152 | 0.7394 | 0.0039 | 0.3235 | 0.3276 | u_span_sec,i_n_users,i_credit_last,i_n_covisit_neighbors,u_n_exp |

## Findings

