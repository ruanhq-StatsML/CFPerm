# FSDS multi-step iteration `iter05/s2`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3329 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 0.4256 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→21 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 1.0000 | 0.5464 | 0.0007 | 0.4846 | 0.8412 | 0.0171 | 0.6457 | 2.5782 | ui_pop_mismatch,i_log1p_n_exp,u_n_exp,i_n_exp,u_n_uniq_items |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→21 | 5-fold π-stable SelectKBest(F) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 2.7340 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |

## Findings

- Best **W2** HGB AUC: `Z_combined|W1hold` = 0.8429 (baseline A=0.8363)
- Best **W1-hold** HGB AUC: `Z_combined|W1hold` = 0.6221
- Scope locked: graph features only; selection steps only; no graph algorithms.
