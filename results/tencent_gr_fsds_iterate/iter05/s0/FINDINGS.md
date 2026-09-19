# FSDS multi-step iteration `iter05/s0`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7679 | 0.0114 | 0.4669 | 0.4278 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 0.4272 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→21 (risk=1.1919e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7664 | 0.0107 | 0.4669 | 2.7163 | ui_pop_mismatch,i_log1p_n_users,i_log1p_n_covisit,u_n_exp,u_n_uniq_items |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→21 | 5-fold π-stable SelectKBest(F) | PO-risk=1.1919e-06 | →FSDS | fused→official_FSDS | 15 | 1.0000 | 0.6621 | 0.0010 | 0.9154 | 0.7735 | 0.0115 | 0.4669 | 2.7288 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |

## Findings

- Best **W2** HGB AUC: `Z_combined|W1hold` = 0.7735 (baseline A=0.7679)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.6621
- Scope locked: graph features only; selection steps only; no graph algorithms.
