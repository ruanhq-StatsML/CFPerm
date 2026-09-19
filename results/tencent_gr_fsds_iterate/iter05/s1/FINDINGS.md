# FSDS multi-step iteration `iter05/s1`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8449 n_va=2700 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000355 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4525 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 0.4938 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→21 (risk=1.25052e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5451 | 0.0012 | 0.5088 | 2.6935 | i_n_users,i_log1p_n_users,i_n_covisit_neighbors,i_log1p_n_covisit,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→21 | 5-fold π-stable SelectKBest(F) | PO-risk=1.25052e-06 | →FSDS | fused→official_FSDS | 15 | 1.0000 | 0.7544 | 0.0010 | 0.1497 | 0.5426 | 0.0012 | 0.5088 | 2.7631 | u_span_sec,i_n_users,i_n_covisit_neighbors,i_credit_last,i_share_last |

## Findings

- Best **W2** HGB AUC: `P_po_vimp_FSDS|W1hold` = 0.5451 (baseline A=0.5426)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.7544
- Scope locked: graph features only; selection steps only; no graph algorithms.
