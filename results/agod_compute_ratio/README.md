# Computational ratio from modality correlation (green)

## Justification

Adapt efficiency claim is **adapt FLOPs**, not inference latency:

```
C_fwd = |M|·c_pf + c_sf          # always paid (FWD on)
C_bwd(A) = |A|·c_pb + c_sb       # A = adapt-active set
flops_rel(A) = (C_fwd + C_bwd(A)) / (C_fwd + C_bwd(M))
```

When modality gradients are **highly correlated** (`cos(g_m,g_{m'})↑`):

1. updates are nearly **collinear** → second tower's proj-BWD buys little new direction
2. define redundancy `ρ_m = mean_{m'≠m} max(cos_{mm'}, 0)` (conflict cos<0 does *not* justify skip)
3. unique mass `u_m = α_m·(1−η·ρ_m)_+` then keep top unique / drop collinear losers
4. **predict** `flops_rel` from the keep set *before* paying BWD

Closed-form savings proxy (no keep set needed):

```
savings_proxy ≈ ρ̄ · (|M|·c_pb) / C_full
flops_rel_proxy ≈ max(1 − savings_proxy, (C_fwd+c_sb)/C_full)
```

This is green by construction: FWD unchanged; only redundant adapt-BWD is skipped.

## Synthetic sweep (equal α)

| pair_cos | mean ρ | flops_proxy | flops_pred | n_keep |
|---|---:|---:|---:|---:|
| -0.20 | 0.00 | 1.000 | 1.000 | 3 |
| -0.03 | 0.00 | 1.000 | 1.000 | 3 |
| +0.49 | 0.49 | 0.720 | 1.000 | 3 |
| +0.72 | 0.72 | 0.589 | 0.619 | 1 |
| +0.89 | 0.89 | 0.490 | 0.619 | 1 |

## From grad-cos smoke trajectories

| Dataset | scheduler | mean ρ | flops_proxy | flops_pred | keep_frac |
|---|---|---:|---:|---:|---:|
| amazon | equal | 0.688 | 0.633 | 0.778 | 0.58 |
| amazon | soft | 0.669 | 0.643 | 0.822 | 0.67 |
| amazon | soft_gradcos | 0.710 | 0.622 | 0.778 | 0.58 |
| msrvtt | equal | 0.192 | 0.890 | 0.810 | 0.67 |
| msrvtt | soft | 0.217 | 0.876 | 0.810 | 0.67 |
| msrvtt | soft_gradcos | 0.215 | 0.877 | 0.810 | 0.67 |

```bash
PYTHONPATH=. python3 scripts/run_agod_compute_ratio_pred.py
```
