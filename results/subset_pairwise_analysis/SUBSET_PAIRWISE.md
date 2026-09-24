# Subset pairwise analysis

Subset pairwise = same edge sort + pairwise(N), but only edges in cohort. Intensity tertiles = buy-volume baseline slices; credit tertiles = last-touch slices. Pack scorecard = sample-count pairwise; calendar scorecard = DiffDB width-min board (different grain).

## Tencent subsets (edge pairwise)

| subset | panel | n | pairs | HGB | LogReg | ΔȲ | Jac | mean_span_h | reading | top0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| all | ops | 16000 | 2 | 0.997 | 0.997 | 0.011 | nan | 1455.7 | transfer holds; feature set partially stable | e_log1p_exp, e_n_exp, i_log1p_n_exp |
| all | content | 16000 | 2 | 0.771 | 0.760 | 0.011 | nan | 1455.7 | transfer holds; feature set partially stable | i_credit_last, i_share_last, i_item_credit_rank |
| intensity_high | ops | 4750 | 2 | 0.994 | 0.997 | -0.003 | nan | 1374.4 | transfer holds; feature set partially stable | e_log1p_exp, e_n_exp, u_ctr |
| intensity_high | content | 4750 | 2 | 0.725 | 0.843 | -0.003 | nan | 1374.4 | transfer holds; feature set partially stable | i_credit_last, i_share_last, i_log1p_n_covisit_neighbors |
| intensity_low | ops | 9034 | 2 | 0.992 | 0.993 | 0.0105 | nan | 1565.8 | transfer holds; feature set partially stable | e_n_exp, e_log1p_exp, i_n_exp |
| intensity_low | content | 9034 | 2 | 0.778 | 0.833 | 0.0105 | nan | 1565.8 | transfer holds; feature set partially stable | i_credit_last, i_share_last, i_item_credit_rank |
| credit_high | ops | 4015 | 3 | 0.998 | 0.997 | -0.014 | nan | 1710.3 | transfer holds; feature set partially stable | e_n_exp, e_log1p_exp, i_share_last |
| credit_high | content | 4015 | 3 | 0.662 | 0.713 | -0.014 | nan | 1710.3 | transfer holds; feature set partially stable | i_share_last, i_credit_last, u_span_sec |
| credit_low | ops | 8684 | 0 | nan | nan | 0 | nan | 1717.7 | no pairs |  |
| credit_low | content | 8684 | 0 | nan | nan | 0 | nan | 1717.7 | no pairs |  |
| I_high_C_high | ops | 2446 | 2 | 0.993 | 0.999 | -0.0065 | nan | 2067.2 | transfer holds; feature set partially stable | e_log1p_exp, e_n_exp, i_share_last |
| I_high_C_high | content | 2446 | 2 | 0.550 | 0.774 | -0.0065 | nan | 2067.2 | weak transfer (association does not travel) | i_share_last, i_credit_last, i_credit_linear |
| I_high_C_low | ops | 1686 | 0 | nan | nan | 0 | nan | 1926.2 | no pairs |  |
| I_high_C_low | content | 1686 | 0 | nan | nan | 0 | nan | 1926.2 | no pairs |  |
| I_low_C_high | ops | 1004 | — | — | — | — | — | — | too_small(need>=1500) | — |
| I_low_C_high | content | 1004 | — | — | — | — | — | — | too_small(need>=1500) | — |
| I_low_C_low | ops | 5770 | 0 | nan | nan | 0 | nan | 1792.8 | no pairs |  |
| I_low_C_low | content | 5770 | 0 | nan | nan | 0 | nan | 1792.8 | no pairs |  |
| item_9091091 | ops | 10 | — | — | — | — | — | — | too_small(need>=1500) | — |
| item_9091091 | content | 10 | — | — | — | — | — | — | too_small(need>=1500) | — |
| item_5082634 | ops | 8 | — | — | — | — | — | — | too_small(need>=1500) | — |
| item_5082634 | content | 8 | — | — | — | — | — | — | too_small(need>=1500) | — |
| item_8500991 | ops | 8 | — | — | — | — | — | — | too_small(need>=1500) | — |
| item_8500991 | content | 8 | — | — | — | — | — | — | too_small(need>=1500) | — |

## Pack pairwise scorecard (sample-count grain)

| pack | N | HGB | LogReg | Jac | ΔȲ | reading |
|---|---:|---:|---:|---:|---:|---|
| diffusiondb | 1000 | 0.605 | 0.616 | 0.24 | 0.01484 | weak transfer (association does not travel) |
| diffusiondb | 2000 | 0.609 | 0.619 | 0.34 | 0.02907 | weak transfer (association does not travel) |
| tencent_gr | 1000 | 0.996 | 0.832 | 0.70 | 0.007833 | strong transfer + stable top feats (persistent drivers) |
| tencent_gr | 2000 | 0.999 | 0.663 | 0.46 | 0.01733 | transfer holds; feature set partially stable |
| tencent_gr_content | 1000 | 0.835 | 0.845 | 0.57 | 0.007833 | strong transfer + stable top feats (persistent drivers) |
| tencent_gr_content | 2000 | 0.900 | 0.898 | 0.43 | 0.01733 | transfer holds; feature set partially stable |
| waymo_proxy | 1000 | 0.940 | 0.946 | 0.72 | 0.2307 | strong transfer + stable top feats (persistent drivers) |
| waymo_proxy | 2000 | 0.969 | 0.978 | 1.00 | 0.4618 | strong transfer + stable top feats (persistent drivers) |
| metro_interstate | 1000 | 0.972 | 0.609 | 0.23 | 20.37 | strong transfer but shifting drivers (regime / composition change) |
| metro_interstate | 2000 | 0.975 | 0.640 | 0.55 | 15.45 | strong transfer + stable top feats (persistent drivers) |
| beijing_pm25 | 1000 | 0.978 | 0.978 | 0.56 | 7.191 | strong transfer + stable top feats (persistent drivers) |
| beijing_pm25 | 2000 | 0.980 | 0.973 | 0.55 | 9.98 | strong transfer + stable top feats (persistent drivers) |

## Calendar-width scorecard (DiffDB — other grain)

| width_min | n_bins | n_reported | mean_auc | mean_ΔȲ |
|---:|---:|---:|---:|---:|
| 30.0 | 647 | 0 | nan | nan |
| 60.0 | 328 | 0 | nan | nan |
| 180.0 | 110 | 2 | 0.49901960784313726 | 0.02816794644741223 |

## Equal-count → calendar side effect (Tencent edges)

Same N edges ≠ same hours — burstiness:

| chunk | n | span_h | y_rate |
|---:|---:|---:|---:|
| 0 | 1000 | 3954.1 | 0.0500 |
| 1 | 1000 | 466.6 | 0.0660 |
| 2 | 1000 | 225.1 | 0.0650 |
| 3 | 1000 | 170.5 | 0.0350 |
| 4 | 1000 | 119.5 | 0.0490 |
| 5 | 1000 | 99.1 | 0.0270 |
| 6 | 1000 | 91.7 | 0.0220 |
| 7 | 1000 | 85.5 | 0.0220 |

## How to read

1. **intensity_high vs low (ops)**: if both stay ~1.0 AUC, intensity baseline is cohort-robust; if only high stays high, baseline is buy-volume concentrated.
2. **content on intensity_low**: last-touch signal without heavy spend — stronger creative/path story.
3. **credit_high content**: should surface `i_credit_*` / `i_share_*` tops.
4. **Pack scorecard**: compare DiffDB weak vs Tencent ops strong vs Metro shifting.
5. **Calendar scorecard**: width-min grain; do not mix cells with equal-count N.
6. **credit_low y_rate→0**: worst credit-rank edges often have no converts — pairwise FS has no pairs (expected); credit_high is the actionable slice.

