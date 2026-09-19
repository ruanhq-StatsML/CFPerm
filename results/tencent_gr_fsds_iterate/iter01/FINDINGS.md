# FSDS multi-step iteration `iter01`

- grids: `W1_localized_grid.parquet` / `W2_localized_grid_sample.parquet`
- n_train=8310 n_va=2839 n_w2=11963 feats=21 k=15

## Variant table

| variant | notes | n_selected | jaccard_vs_baseline | W1_hgb_auc | W1_logreg_auc | W2_hgb_auc | W2_logreg_auc | sec | top5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_baseline_F|W1hold | Scaler→Var→SelectKBest(F) | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7679 | 0.4669 | 0.4064 | e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec |
| C_cmean_then_F|W1hold | cmean|δ| prefilter→30 then F→15 | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7664 | 0.4669 | 0.4300 | u_span_sec,ui_pop_mismatch,i_n_exp,u_n_events,u_n_uniq_items |
| D_F_corr_prune|W1hold | F wide=21 → corr-prune@0.92 →10 | 10 | 0.4706 | 0.4924 | 0.8721 | 0.6437 | 0.7117 | 0.4502 | e_n_exp,u_n_exp,u_span_sec,u_log1p_n_uniq,i_n_covisit_neighbors |
| E_stable_pi_F|W1hold | 5-fold π-stable SelectKBest(F) | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7735 | 0.4669 | 0.4978 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| F_cmean_stable|W1hold | cmean pre→30 + π-stable F→15 | 15 | 1.0000 | 0.6621 | 0.9154 | 0.7735 | 0.4669 | 0.4844 | u_n_uniq_items,ui_pop_mismatch,e_log1p_exp,e_n_exp,i_log1p_n_exp |
| G_F_then_HGB_perm|W1hold | F screen=21 → HGB perm-imp →15 | 15 | 0.5789 | 0.6434 | 0.9468 | 0.7730 | 0.7042 | 2.1339 | i_n_exp,i_share_linear,ui_pop_mismatch,u_n_exp,u_span_sec |

## Findings

- Best **W2** HGB AUC: `F_cmean_stable|W1hold` = 0.7735 (baseline A=0.7679)
- Best **W1-hold** HGB AUC: `A_baseline_F|W1hold` = 0.6621
- Corr-prune Jaccard vs baseline=0.471 (drops near-duplicate share/credit twins when present)
- cmean-|δ| prefilter aligns FS with period-shift guidance (same δ story as drill doc)
- Scope locked: graph features only; selection steps only; no graph algorithms.
