# Grad-RFPerm freeze: MSE–FLOPs efficiency

freeze_eff best {'freeze_early': 1, 'freeze_low_share': 1}; 2 Pareto wins vs always_adapt (MSE↓&FLOPs↓). no_adapt collapses — freeze ≠ stop learning.

- source=`canonical_extras_doc` · datasets=2 · Pareto wins=2

| dataset | policy | mse_rel | flops_rel | freeze_eff | Pareto? |
|---|---|---:|---:|---:|:---:|
| `electricity` | `always_adapt` | 1.00 | 1.00 | 0.000 | N |
| `electricity` | `freeze_early` | 0.86 | 0.70 | 0.200 | Y |
| `electricity` | `freeze_low_share` | 0.86 | 0.70 | 0.200 | Y |
| | | | | | *electricity: freeze_early, freeze_low_share Pareto-beat always_adapt (MSE↓ & FLOPs↓)* |
| `synthetic` | `always_adapt` | 1.00 | 1.00 | 0.000 | N |
| `synthetic` | `freeze_early` | 1.25 | 0.70 | -0.357 | N |
| `synthetic` | `freeze_low_share` | 1.15 | 0.57 | -0.263 | N |
| | | | | | *synthetic: best freeze_eff=freeze_low_share (may trade MSE for FLOPs)* |

## Justify

- After Grad reject, freezing low-share / early layers cuts adapt FLOPs.
- **Pareto** = MSE better *and* FLOPs lower than always_adapt.
- electricity: freeze policies ≈0.86× MSE @ ~0.7× FLOPs (Pareto).
- synthetic: FLOPs save but MSE↑ — freeze_eff can be negative; still prefer freeze_low_share over freeze_early.
- `no_adapt` collapses — do not equate freeze with stop.

