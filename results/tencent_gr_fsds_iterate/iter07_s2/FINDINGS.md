# FSDS multi-step iteration `iter07_s2`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3235 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6922 | 0.0010 | 0.5831 | 0.8710 | 0.0136 | 0.5153 | 2.6972 | ui_pop_mismatch,i_n_users,i_log1p_n_exp,u_n_exp,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO(α=0.5): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7524 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| Z_combined_PO_a03|W1hold | COMBINED-PO(α=0.3): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7461 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| Z_combined_PO_a07|W1hold | COMBINED-PO(α=0.7): cmean⋈PO-VIMP→18 | 3-fold π-stable SelectKBest(F) (n_pos=3) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7538 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| P_po_vimp_rare_pi|W1hold | PO-VIMP+rare-π: pre→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7437 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7643 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |

## Findings

- Best **W2** HGB AUC: `P_po_vimp_FSDS|W1hold` = 0.8710 (baseline A=0.8363)
- Best **W1-hold** HGB AUC: `P_po_vimp_FSDS|W1hold` = 0.6922
- Scope locked: graph features only; selection steps only; no graph algorithms.
