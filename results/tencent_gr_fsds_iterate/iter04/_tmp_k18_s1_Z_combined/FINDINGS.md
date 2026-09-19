# FSDS multi-step iteration `iter04/_tmp_k18_s1_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=18
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→18 | fused→official_FSDS | 18 |  | 0.7666 | 0.0011 | 0.1738 | 0.5286 | 0.0012 | 0.3323 | 0.5287 | u_span_sec,i_n_users,i_n_covisit_neighbors,u_n_exp,u_n_uniq_items |

## Findings

