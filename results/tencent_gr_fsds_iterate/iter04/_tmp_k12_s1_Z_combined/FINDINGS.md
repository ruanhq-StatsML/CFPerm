# FSDS multi-step iteration `iter04/_tmp_k12_s1_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=12
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→12 | fused→official_FSDS | 12 |  | 0.5921 | 0.0007 | 0.1397 | 0.7389 | 0.0039 | 0.4089 | 0.4094 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |

## Findings

