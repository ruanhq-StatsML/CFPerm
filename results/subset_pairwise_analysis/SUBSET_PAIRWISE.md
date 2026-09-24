# Subset pairwise analysis

Subset pairwise = same edge sort + pairwise(N), but only edges in cohort. Intensity tertiles = buy-volume baseline slices; credit tertiles = last-touch slices. Pack scorecard = sample-count pairwise; calendar scorecard = DiffDB width-min board (different grain).

## Tencent subsets (edge pairwise)

| subset | panel | n | pairs | HGB | LogReg | ΔȲ | Jac | mean_span_h | reading | top0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| all | ops | 12000 | 6 | 0.994 | 0.996 | -0.0035 | 0.60 | 585.8 | strong transfer + stable top feats (persistent drivers) | e_log1p_exp, e_n_exp, i_log1p_n_exp |
| all | content | 12000 | 6 | 0.800 | 0.798 | -0.0035 | 0.66 | 585.8 | strong transfer + stable top feats (persistent drivers) | i_credit_last, i_share_last, i_item_credit_rank |
| intensity_high | ops | 3539 | 2 | 0.994 | 0.994 | -0.009 | 0.67 | 1313.8 | strong transfer + stable top feats (persistent drivers) | e_log1p_exp, e_n_exp, i_credit_last |
| intensity_high | content | 3539 | 2 | 0.811 | 0.906 | -0.009 | 0.67 | 1313.8 | strong transfer + stable top feats (persistent drivers) | i_credit_last, i_share_last, i_log1p_credit_linear |
| intensity_low | ops | 6744 | 4 | 0.990 | 0.991 | -0.015 | 0.78 | 1004.3 | strong transfer + stable top feats (persistent drivers) | e_log1p_exp, e_n_exp, i_log1p_n_exp |
| intensity_low | content | 6744 | 4 | 0.811 | 0.840 | -0.015 | 0.43 | 1004.3 | transfer holds; feature set partially stable | i_credit_last, i_share_last, i_item_credit_rank |
| credit_high | ops | 3031 | 2 | 0.993 | 0.994 | -0.0205 | 0.67 | 1597.0 | strong transfer + stable top feats (persistent drivers) | e_log1p_exp, e_n_exp, i_share_last |
| credit_high | content | 3031 | 2 | 0.723 | 0.767 | -0.0205 | 0.25 | 1597.0 | strong transfer but shifting drivers (regime / composition change) | i_share_last, i_credit_last, i_credit_first |
| credit_low | ops | 6530 | 0 | nan | nan | 0 | nan | 1004.4 | no pairs |  |
| credit_low | content | 6530 | 0 | nan | nan | 0 | nan | 1004.4 | no pairs |  |

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
| 0 | 1000 | 4160.6 | 0.0700 |
| 1 | 1000 | 420.6 | 0.0610 |
| 2 | 1000 | 235.0 | 0.0460 |
| 3 | 1000 | 156.0 | 0.0440 |
| 4 | 1000 | 126.9 | 0.0180 |
| 5 | 1000 | 113.0 | 0.0220 |
| 6 | 1000 | 97.6 | 0.0210 |
| 7 | 1000 | 84.4 | 0.0140 |

## How to read

1. **intensity_high vs low (ops)**: if both stay ~1.0 AUC, intensity baseline is cohort-robust; if only high stays high, baseline is buy-volume concentrated.
2. **content on intensity_low**: last-touch signal without heavy spend — stronger creative/path story.
3. **credit_high content**: should surface `i_credit_*` / `i_share_*` tops.
4. **Pack scorecard**: compare DiffDB weak vs Tencent ops strong vs Metro shifting.
5. **Calendar scorecard**: width-min grain; do not mix cells with equal-count N.
6. **credit_low y_rate→0**: worst credit-rank edges often have no converts — pairwise FS has no pairs (expected); credit_high is the actionable slice.

