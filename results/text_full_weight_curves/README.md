# Feature vs sample weights — FULL holdout acc / AUC curves

Feature-level = X·(1+α·vimp); sample-level = membership w_i; both = product. Curves on FULL LIVE holdout overall acc + novel OVA-AUC. No EWMA.

KPI = **overall accuracy + mean novel OVA-AUC on full held-out LIVE** (not first-n). budgets=[20, 40, 60, 80, 100, 150, 200], seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, α_feat=3.0, **EWMA=off**.

## Endpoint overall accuracy (full eval)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.378 | 0.386 | 0.413 | 0.389 | 0.008 | 0.035 | 0.011 |
| news20_rec | 0.458 | 0.457 | 0.473 | 0.483 | -0.001 | 0.015 | 0.025 |
| news20_mix | 0.283 | 0.281 | 0.320 | 0.303 | -0.003 | 0.036 | 0.020 |
| cnae9_text | 0.863 | 0.889 | 0.878 | 0.888 | 0.027 | 0.015 | 0.025 |
| amazon_cat | 0.692 | 0.705 | 0.695 | 0.724 | 0.013 | 0.003 | 0.031 |

mean overall acc: uniform=0.535  feature_w=0.544  sample_w=0.556  both=0.557
mean Δacc vs U: feature_w=0.009  sample_w=0.021  both=0.023

## Endpoint mean novel OVA-AUC (full eval)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.683 | 0.674 | 0.682 | 0.675 | -0.010 | -0.002 | -0.009 |
| news20_rec | 0.693 | 0.694 | 0.688 | 0.685 | 0.001 | -0.005 | -0.008 |
| news20_mix | 0.696 | 0.669 | 0.702 | 0.681 | -0.028 | 0.005 | -0.015 |
| cnae9_text | 0.961 | 0.967 | 0.959 | 0.953 | 0.005 | -0.003 | -0.008 |
| amazon_cat | 0.808 | 0.819 | 0.810 | 0.814 | 0.011 | 0.002 | 0.006 |

mean OVA-AUC: uniform=0.768  feature_w=0.764  sample_w=0.768  both=0.762
mean ΔOVA-AUC vs U: feature_w=-0.004  sample_w=-0.000  both=-0.007

## Mean curves (acc / OVA-AUC vs novel budget)

| budget | acc U | acc feat | acc sample | acc both | AUC U | AUC feat | AUC sample | AUC both |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.312 | 0.312 | 0.348 | 0.341 | 0.624 | 0.616 | 0.637 | 0.619 |
| 40 | 0.377 | 0.384 | 0.425 | 0.429 | 0.674 | 0.672 | 0.689 | 0.688 |
| 60 | 0.425 | 0.437 | 0.468 | 0.469 | 0.703 | 0.705 | 0.713 | 0.708 |
| 80 | 0.450 | 0.465 | 0.486 | 0.493 | 0.719 | 0.722 | 0.725 | 0.725 |
| 100 | 0.472 | 0.484 | 0.493 | 0.511 | 0.730 | 0.733 | 0.734 | 0.733 |
| 150 | 0.514 | 0.532 | 0.540 | 0.547 | 0.763 | 0.763 | 0.764 | 0.762 |
| 200 | 0.535 | 0.544 | 0.556 | 0.557 | 0.768 | 0.764 | 0.768 | 0.762 |

```bash
python3 scripts/run_text_full_weight_curves.py --multiseed
```
