# Ensemble decorrelation under high modality correlation

## Justification

High pairwise grad cosine ⇒ modalities are **not** independent voters; the
cosine Gram is near-rank-1. Skip alone is blunt. Decorrelate the *ensemble*:

1. `effective_rank(G)` from pairwise `cos(g_m,g_{m'})` — near 1 ⇒ collinear
2. **leader** = argmax unique mass `α_m·(1−η·ρ_m)_+` (owns shared direction)
3. **diversifier** = residual uniqueness `r_m=(1−max(cos,0))_+` after leader
4. **redundant** = collinear with leader → damp adapt LR (FWD still on)

LR map (`soft_decorr`):

```
LR_m = α_to_lr(α_m · role_gain_m)
role_gain: leader↑ · diversifier·(1+boost·r_m) · redundant→floor
```

When `mean ρ` low / effective rank high → decorr inactive → reduces to soft LR.

## Synthetic isotropic sweep

| pair_cos | mean ρ | erank | decorr | leader | roles(v/t/a) | LR ratio |
|---|---:|---:|---|---|---|---:|
| +0.01 | 0.01 | 3.00 | False | video | D/D/D | 1.00 |
| +0.43 | 0.43 | 2.54 | False | video | D/D/D | 1.00 |
| +0.53 | 0.53 | 2.31 | False | video | D/D/D | 1.00 |
| +0.74 | 0.74 | 1.79 | True | video | L/D/D | 1.19 |
| +0.90 | 0.90 | 1.35 | True | video | L/R/R | 5.40 |

## Asymmetric cases (unequal α)

| case | mean ρ | erank | leader | roles | LR |
|---|---:|---:|---|---|---|
| `low_corr` | 0.15 | 2.94 | video | video=diversifier, text=diversifier, audio=diversifier | video=1.32, text=1.05, audio=0.64 |
| `high_iso` | 0.80 | 1.62 | video | video=leader, text=redundant, audio=redundant | video=1.78, text=0.26, audio=0.16 |
| `leader_plus_residual` | 0.63 | 1.77 | text | video=redundant, text=leader, audio=diversifier | video=0.33, text=1.41, audio=0.90 |
| `full_collinear` | 0.92 | 1.28 | video | video=leader, text=redundant, audio=redundant | video=1.78, text=0.26, audio=0.16 |

## From grad-cos trajectories (role replay)

| Dataset | source sched | mean ρ | erank | frac decorr | role frac L/D/R |
|---|---|---:|---:|---:|---|
| amazon | equal | 0.688 | 1.54 | 1.00 | 0.50/0.42/0.08 |
| amazon | soft | 0.669 | 1.56 | 1.00 | 0.50/0.42/0.08 |
| amazon | soft_gradcos | 0.710 | 1.51 | 1.00 | 0.50/0.33/0.17 |
| msrvtt | equal | 0.192 | 2.85 | 0.00 | 0.00/1.00/0.00 |
| msrvtt | soft | 0.217 | 2.82 | 0.00 | 0.00/1.00/0.00 |
| msrvtt | soft_gradcos | 0.215 | 2.83 | 0.00 | 0.00/1.00/0.00 |

## Online smoke (optional)

| Dataset | scheduler | mean pair_cos | frac decorr | mean Acc↑ | mean LR ratio |
|---|---|---:|---:|---:|---:|
| msrvtt | equal | 0.190 | 0.00 | +0.010 | 1.00 |
| msrvtt | soft | 0.216 | 0.00 | +0.006 | 2.83 |
| msrvtt | soft_decorr | 0.216 | 0.00 | +0.006 | 2.83 |

```bash
PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py
PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py --online msrvtt
```
