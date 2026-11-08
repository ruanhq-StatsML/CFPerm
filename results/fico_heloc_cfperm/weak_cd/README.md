# FICO HELOC — weak-5 + highly nonlinear concept drift

Design justification: HELOC is driven by **several weak features**, not one
strong median axis. We randomly sample `K=5` features, keep `P(X)` identical
across batches, and redraw `Y` on the new batch from a highly nonlinear map
of those five (sin/tanh/interactions/cubic) → pure concept drift.

- Clean n=`9861`, draws=`5`, intensity=`2.5`, `B=150`, `n_perm=200`, `max_n_per_arm=1500`.

## Power (this design)

| scope | reject | any_driver@top5 | mean recall@5 | n |
|---|---:|---:|---:|---:|
| overall | 0.867 | 1.000 | 0.493 | 15 |
| loco | 0.600 | 1.000 | 0.520 | 5 |
| permucate | 1.000 | 1.000 | 0.480 | 5 |
| grf | 1.000 | 1.000 | 0.480 | 5 |

## vs median-split covariate shift (supplement)

| design | scope | reject | recover / driver-hit |
|---|---|---:|---:|
| median_split_covshift | overall | 0.222 | 0.111 |
| median_split_covshift | loco | 0.000 | 0.000 |
| median_split_covshift | permucate | 0.333 | 0.000 |
| median_split_covshift | grf | 0.333 | 0.333 |
| weak5_nonlinear_cd | overall | 0.867 | 1.000 |
| weak5_nonlinear_cd | loco | 0.600 | 1.000 |
| weak5_nonlinear_cd | permucate | 1.000 | 1.000 |
| weak5_nonlinear_cd | grf | 1.000 | 1.000 |

### Takeaway

- Median-split on one LR feature = **single strong covariate-shift** axis;
  power concentrates on GRF recovery of that axis.
- Weak-5 nonlinear CD = **several weak drivers** of `P(Y|X)`; rejection +
  recall@5 justify CFPerm under the multi-weak-feature regime that matches HELOC.

Files: `weak_cd_cells.csv`, `weak_cd_power.csv`, `power_comparison.csv`,
`weak_cd_posthoc_*.csv`.

