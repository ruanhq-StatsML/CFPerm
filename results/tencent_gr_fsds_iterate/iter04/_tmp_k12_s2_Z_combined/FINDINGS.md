# FSDS multi-step iteration `iter04/_tmp_k12_s2_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=12
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→12 | fused→official_FSDS | 12 |  | 0.6490 | 0.0008 | 0.5457 | 0.7786 | 0.0081 | 0.6872 | 0.4064 | i_n_exp,e_log1p_exp,e_n_exp,i_share_first,i_credit_first |

## Findings

