# Correlation buckets + R-run frequency (Amazon / MSR-VTT)

## Protocol

Per online window from grad-cos trajectories:

1. rebuild pair geometry from `mean_pair_cos` (or `pair_cos` if logged)
2. `soft_decorr` role split → `ρ̄`, `erank`, roles, LR
3. **bucket** by correlation:
   - **high**: `ρ̄ ≥ 0.55` OR `erank ≤ 1.55` (decorr trigger-aligned)
   - **low**: `ρ̄ < 0.35` AND `erank > 2.20` (independent voters)
   - **mid**: borderline
4. **R-run**: contiguous `role=redundant` while decorr on, length ≥ 2 (sticky);
   also count length-1 R events

Budget actuator stays **LR** (`γ_role`); FWD always on. This script characterizes
adjust-correlation geometry on logged streams (no retrain).

## Summary

| Dataset | source | mean ρ | erank | bucket L/M/H | frac decorr | frac any-R | sticky R-runs (rate/high) | R events (len≥1) | LR ratio (high) | Acc↑ source |
|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|
| amazon | equal | 0.688 | 1.538 | 0.000/0.000/1.000 | 1.000 | 0.167 | 0 (0.000/high) | 1 | 7.127 | -0.009 |
| amazon | soft | 0.669 | 1.563 | 0.000/0.000/1.000 | 1.000 | 0.167 | 0 (0.000/high) | 1 | 6.904 | -0.007 |
| amazon | soft_gradcos | 0.710 | 1.511 | 0.000/0.000/1.000 | 1.000 | 0.333 | 0 (0.000/high) | 2 | 8.332 | 0.000 |
| msrvtt | equal | 0.192 | 2.853 | 0.833/0.167/0.000 | 0.000 | 0.000 | 0 (0.000/high) | 0 | — | 0.010 |
| msrvtt | soft | 0.217 | 2.823 | 0.833/0.167/0.000 | 0.000 | 0.000 | 0 (0.000/high) | 0 | — | 0.006 |
| msrvtt | soft_gradcos | 0.215 | 2.829 | 0.833/0.167/0.000 | 0.000 | 0.000 | 0 (0.000/high) | 0 | — | -0.006 |

## High-bucket role mix

| Dataset | source | n_high | role frac L/D/R | leader-swap | Acc↑ in high |
|---|---|---:|---|---:|---:|
| amazon | equal | 6 | 0.500/0.417/0.083 | 0.000 | -0.009 |
| amazon | soft | 6 | 0.500/0.417/0.083 | 0.000 | -0.007 |
| amazon | soft_gradcos | 6 | 0.500/0.333/0.167 | 0.000 | 0.000 |

## Readout (honest)

- **Amazon**: windows sit in **high** bucket (ρ̄≈0.67–0.71, erank≈1.5); decorr always on.
  Length-1 R events appear; sticky R-runs (k≥2) are rare on this short smoke (n=6).
  High-bucket LR ratio is large → soft_decorr reallocates step budget.
- **MSR-VTT**: mostly **low** (+ occasional mid); decorr off; R-run rate = 0
  → negative control for adjust-correlation.

```bash
PYTHONPATH=. python3 scripts/run_agod_decorr_buckets.py
```
