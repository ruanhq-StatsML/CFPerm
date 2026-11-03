# Feature vs sample weights — FULL holdout acc / AUC curves

Cold-start empowerment only: after a novel category reaches 200 samples, stop weighting it. Feature-level = X·(1+α·vimp); sample-level = membership w_i. FULL LIVE holdout. No EWMA.

KPI = **overall accuracy + mean novel OVA-AUC on full held-out LIVE** (not first-n). budget=novel/category [20, 40, 60, 80, 100, 150, 200, 250, 300, 400], cold_start_cap=**200** (stop weighting past cap), seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, α_feat=3.0, **EWMA=off**.

## At cold-start cap (budget=200/cat) — overall acc

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.493 | 0.456 | 0.471 | 0.485 | -0.037 | -0.022 | -0.008 |
| news20_rec | 0.534 | 0.542 | 0.528 | 0.533 | 0.007 | -0.007 | -0.001 |
| news20_mix | 0.476 | 0.470 | 0.470 | 0.473 | -0.006 | -0.006 | -0.004 |
| cnae9_text | 0.892 | 0.917 | 0.895 | 0.909 | 0.025 | 0.003 | 0.016 |
| amazon_cat | 0.760 | 0.751 | 0.748 | 0.754 | -0.009 | -0.011 | -0.006 |

## Endpoint overall accuracy (full eval, last budget)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.502 | 0.499 | 0.498 | 0.498 | -0.003 | -0.004 | -0.003 |
| news20_rec | 0.545 | 0.558 | 0.538 | 0.553 | 0.013 | -0.007 | 0.009 |
| news20_mix | 0.504 | 0.506 | 0.502 | 0.501 | 0.002 | -0.002 | -0.003 |
| cnae9_text | 0.894 | 0.913 | 0.906 | 0.909 | 0.019 | 0.012 | 0.015 |
| amazon_cat | 0.737 | 0.764 | 0.755 | 0.757 | 0.027 | 0.019 | 0.020 |

mean overall acc: uniform=0.636  feature_w=0.648  sample_w=0.640  both=0.644  both_hard=0.649  both_always=0.652
mean Δacc vs U: feature_w=0.012  sample_w=0.004  both=0.007  both_hard=0.013  both_always=0.016

## Endpoint mean novel OVA-AUC (full eval)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.783 | 0.773 | 0.771 | 0.775 | -0.010 | -0.012 | -0.008 |
| news20_rec | 0.749 | 0.748 | 0.739 | 0.742 | -0.000 | -0.010 | -0.007 |
| news20_mix | 0.841 | 0.838 | 0.843 | 0.840 | -0.002 | 0.002 | -0.001 |
| cnae9_text | 0.968 | 0.978 | 0.977 | 0.965 | 0.010 | 0.009 | -0.003 |
| amazon_cat | 0.859 | 0.870 | 0.872 | 0.863 | 0.010 | 0.013 | 0.004 |

mean OVA-AUC: uniform=0.840  feature_w=0.841  sample_w=0.841  both=0.837  both_hard=0.835  both_always=0.840
mean ΔOVA-AUC vs U: feature_w=0.001  sample_w=0.001  both=-0.003  both_hard=-0.005  both_always=-0.000

## Mean curves (acc / OVA-AUC vs novel budget)

| budget | acc U | acc feat | acc sample | acc both | AUC U | AUC feat | AUC sample | AUC both |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.456 | 0.471 | 0.464 | 0.477 | 0.729 | 0.724 | 0.725 | 0.730 |
| 40 | 0.534 | 0.524 | 0.539 | 0.542 | 0.769 | 0.758 | 0.765 | 0.760 |
| 60 | 0.567 | 0.578 | 0.572 | 0.580 | 0.781 | 0.784 | 0.781 | 0.783 |
| 80 | 0.568 | 0.572 | 0.579 | 0.576 | 0.779 | 0.777 | 0.784 | 0.780 |
| 100 | 0.595 | 0.598 | 0.593 | 0.597 | 0.806 | 0.803 | 0.800 | 0.799 |
| 150 | 0.624 | 0.625 | 0.623 | 0.627 | 0.820 | 0.820 | 0.820 | 0.818 |
| 200 | 0.631 | 0.627 | 0.623 | 0.631 | 0.829 | 0.825 | 0.823 | 0.824 |
| 250 | 0.631 | 0.649 | 0.636 | 0.645 | 0.833 | 0.839 | 0.835 | 0.833 |
| 300 | 0.643 | 0.653 | 0.644 | 0.660 | 0.837 | 0.842 | 0.835 | 0.842 |
| 400 | 0.636 | 0.648 | 0.640 | 0.644 | 0.840 | 0.841 | 0.841 | 0.837 |

```bash
python3 scripts/run_text_full_weight_curves.py --multiseed
```
