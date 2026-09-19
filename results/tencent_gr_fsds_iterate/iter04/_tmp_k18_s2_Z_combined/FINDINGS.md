# FSDS multi-step iteration `iter04/_tmp_k18_s2_Z_combined`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=18
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→18 | fused→official_FSDS | 18 |  | 0.6929 | 0.0010 | 0.5556 | 0.8671 | 0.0131 | 0.7212 | 0.4480 | e_log1p_exp,e_n_exp,i_share_first,i_credit_first,i_n_exp |

## Findings

