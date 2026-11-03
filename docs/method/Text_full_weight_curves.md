# Feature vs sample weights — FULL holdout acc / AUC curves

polish: OOD+align weights; recommended switches align off after cold-start frac drops (soft^γ=0.5). FULL holdout. No EWMA.

KPI = **overall accuracy + mean novel OVA-AUC on full held-out LIVE** (not first-n). budget=novel/category [20, 40, 60, 80, 100, 150, 200, 250, 300, 400], cold_start_cap=**200** (stop weighting past cap), seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, α_feat=3.0, **EWMA=off**.

## At cold-start cap (budget=200/cat) — overall acc

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.493 | 0.489 | 0.458 | 0.483 | -0.004 | -0.034 | -0.009 |
| news20_rec | 0.534 | 0.543 | 0.528 | 0.536 | 0.009 | -0.007 | 0.002 |
| news20_mix | 0.476 | 0.471 | 0.477 | 0.469 | -0.006 | 0.001 | -0.007 |
| cnae9_text | 0.892 | 0.916 | 0.892 | 0.904 | 0.024 | 0.000 | 0.012 |
| amazon_cat | 0.760 | 0.742 | 0.744 | 0.753 | -0.017 | -0.016 | -0.007 |

## Endpoint overall accuracy (full eval, last budget)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.502 | 0.497 | 0.512 | 0.512 | -0.005 | 0.010 | 0.011 |
| news20_rec | 0.545 | 0.555 | 0.545 | 0.561 | 0.010 | 0.000 | 0.016 |
| news20_mix | 0.504 | 0.498 | 0.509 | 0.505 | -0.006 | 0.005 | 0.001 |
| cnae9_text | 0.894 | 0.906 | 0.904 | 0.901 | 0.012 | 0.010 | 0.007 |
| amazon_cat | 0.737 | 0.753 | 0.748 | 0.757 | 0.016 | 0.011 | 0.020 |

mean overall acc: uniform=0.636  feature_w=0.642  sample_w=0.643  both=0.647  align_w=0.644  both_align=0.650  recommended=0.650  both_always=0.652
mean Δacc vs U: feature_w=0.005  sample_w=0.007  both=0.011  align_w=0.008  both_align=0.013  recommended=0.014  both_always=0.016

## Endpoint mean novel OVA-AUC (full eval)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.783 | 0.771 | 0.773 | 0.761 | -0.012 | -0.010 | -0.022 |
| news20_rec | 0.749 | 0.747 | 0.746 | 0.747 | -0.002 | -0.003 | -0.002 |
| news20_mix | 0.841 | 0.838 | 0.844 | 0.836 | -0.003 | 0.003 | -0.005 |
| cnae9_text | 0.968 | 0.967 | 0.974 | 0.972 | -0.001 | 0.006 | 0.004 |
| amazon_cat | 0.859 | 0.861 | 0.848 | 0.865 | 0.002 | -0.012 | 0.006 |

mean OVA-AUC: uniform=0.840  feature_w=0.837  sample_w=0.837  both=0.836  align_w=0.838  both_align=0.836  recommended=0.836  both_always=0.840
mean ΔOVA-AUC vs U: feature_w=-0.003  sample_w=-0.003  both=-0.004  align_w=-0.002  both_align=-0.004  recommended=-0.004  both_always=-0.000

## Mean curves (acc / OVA-AUC vs novel budget)

| budget | acc U | acc feat | acc sample | acc both | AUC U | AUC feat | AUC sample | AUC both |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.456 | 0.463 | 0.469 | 0.475 | 0.729 | 0.728 | 0.729 | 0.725 |
| 40 | 0.534 | 0.529 | 0.540 | 0.540 | 0.769 | 0.762 | 0.766 | 0.759 |
| 60 | 0.567 | 0.581 | 0.572 | 0.585 | 0.781 | 0.782 | 0.782 | 0.787 |
| 80 | 0.568 | 0.570 | 0.584 | 0.577 | 0.779 | 0.780 | 0.784 | 0.784 |
| 100 | 0.595 | 0.596 | 0.595 | 0.604 | 0.806 | 0.803 | 0.800 | 0.798 |
| 150 | 0.624 | 0.618 | 0.618 | 0.626 | 0.820 | 0.814 | 0.820 | 0.818 |
| 200 | 0.631 | 0.632 | 0.620 | 0.629 | 0.829 | 0.822 | 0.822 | 0.823 |
| 250 | 0.631 | 0.643 | 0.637 | 0.646 | 0.833 | 0.837 | 0.836 | 0.836 |
| 300 | 0.643 | 0.654 | 0.644 | 0.648 | 0.837 | 0.845 | 0.832 | 0.836 |
| 400 | 0.636 | 0.642 | 0.643 | 0.647 | 0.840 | 0.837 | 0.837 | 0.836 |

```bash
python3 scripts/run_text_full_weight_curves.py --multiseed
```
