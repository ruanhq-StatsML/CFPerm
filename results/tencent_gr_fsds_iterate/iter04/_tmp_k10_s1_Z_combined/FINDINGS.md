# FSDS multi-step iteration `iter04/_tmp_k10_s1_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=10
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→20 + π-stable F→10 | fused→official_FSDS | 10 |  | 0.5572 | 0.0007 | 0.1149 | 0.6060 | 0.0013 | 0.3616 | 0.2787 | u_span_sec,i_n_users,i_credit_last,i_share_last,e_log1p_exp |

## Findings

