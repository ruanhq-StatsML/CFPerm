# FSDS multi-step iteration `iter04/_tmp_k15_s1_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→15 | fused→official_FSDS | 15 |  | 0.6250 | 0.0010 | 0.1552 | 0.4158 | 0.0010 | 0.4164 | 1.7833 | e_log1p_exp,i_share_first,i_log1p_n_users,u_log1p_n_uniq,i_share_linear |

## Findings

