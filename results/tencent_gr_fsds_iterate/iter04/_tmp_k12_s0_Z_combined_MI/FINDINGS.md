# FSDS multi-step iteration `iter04/_tmp_k12_s0_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=12
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→21 + π-stable MI→12 | fused→official_FSDS | 12 |  | 0.5055 | 0.0007 | 0.6293 | 0.3465 | 0.0019 | 0.7972 | 1.8594 | e_log1p_exp,e_n_exp,i_share_last,i_share_linear,i_share_first |

## Findings

