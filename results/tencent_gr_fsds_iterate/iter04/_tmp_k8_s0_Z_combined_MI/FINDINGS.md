# FSDS multi-step iteration `iter04/_tmp_k8_s0_Z_combined_MI`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=8
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Z_combined_MI|W1hold | COMBINED-MI(=cmean+π-MI→FSDS): cmean pre→16 + π-stable MI→8 | fused→official_FSDS | 8 |  | 0.4190 | 0.0006 | 0.6022 | 0.5969 | 0.0015 | 0.4512 | 1.4762 | u_log1p_n_uniq,i_credit_last,i_log1p_n_exp,ui_pop_mismatch,u_n_events |

## Findings

