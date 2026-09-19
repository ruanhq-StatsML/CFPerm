# FSDS multi-step iteration `iter06_s2`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8315 n_va=2834 n_w2=11963 feats=21 k=15
- fuse_official_fsds=True (Scaler→Var→SelectKBest→model; W2 never selects)
- pos_train=0.000361 (rare-positive regime)

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_hgb_ap | W1_logreg_auc | W2_hgb_auc | W2_hgb_ap | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | metrics=official_run_fsds | 15 | 1.0000 | 0.5459 | 0.0007 | 0.4846 | 0.8363 | 0.0171 | 0.6457 | 0.3257 | e_n_exp,u_n_exp,u_n_uniq_items,u_log1p_n_events,u_log1p_n_uniq |
| Z_combined|W1hold | COMBINED(=cmean+π→FSDS): cmean pre→21 + π-stable F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6221 | 0.0005 | 0.5429 | 0.8429 | 0.0090 | 0.7445 | 0.4337 | i_share_first,i_credit_first,i_n_exp,e_log1p_exp,e_n_exp |
| P_po_vimp_FSDS|W1hold | PO-VIMP pre→18 (risk=1.15413e-06) then FSDS-F→15 | fused→official_FSDS | 15 | 0.6667 | 0.6922 | 0.0010 | 0.5831 | 0.8710 | 0.0136 | 0.5153 | 2.6927 | ui_pop_mismatch,i_n_users,i_log1p_n_exp,u_n_exp,i_n_exp |
| Z_combined_PO|W1hold | COMBINED-PO: cmean⋈PO-VIMP→18 | 5-fold π-stable SelectKBest(F) | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5789 | 0.6929 | 0.0010 | 0.6322 | 0.8671 | 0.0131 | 0.6211 | 2.7049 | i_share_first,i_credit_first,i_n_exp,u_span_sec,i_log1p_n_covisit |
| P_po_boot_pi|W1hold | PO-boot-π: PO-help pool→18 | bootπ×40→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6922 | 0.0010 | 0.7010 | 0.8707 | 0.0132 | 0.8291 | 2.8500 | i_share_first,i_credit_first,i_credit_last,ui_pop_mismatch,u_log1p_n_events |
| Z_combined_PO_rare|W1hold | COMBINED-PO-rare: cmean⋈PO→18 | 3-fold π (n_pos=3)→15 | PO-risk=1.15413e-06 | →FSDS | fused→official_FSDS | 15 | 0.5000 | 0.6221 | 0.0005 | 0.6841 | 0.8429 | 0.0090 | 0.6567 | 2.7843 | i_share_first,i_credit_first,i_log1p_n_exp,i_n_exp,u_span_sec |

## Findings

- Best **W2** HGB AUC: `P_po_vimp_FSDS|W1hold` = 0.8710 (baseline A=0.8363)
- Best **W1-hold** HGB AUC: `Z_combined_PO|W1hold` = 0.6929
- Scope locked: graph features only; selection steps only; no graph algorithms.
