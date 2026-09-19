# FSDS multi-step iteration `iter04/_tmp_k12_s0_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=12
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→12 | fused→official_FSDS | 12 |  | 0.5913 | 0.0009 | 0.8252 | 0.3938 | 0.0016 | 0.3730 | 0.2889 | e_log1p_exp,e_n_exp,i_log1p_n_exp,i_share_first,i_credit_first |

## Findings

