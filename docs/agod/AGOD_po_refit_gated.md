# Rolling PO-learner refit (RF / XGBoost / MLP)

Not Ridge. Same family for the DR PO-learner, the residual gate,
and the next-batch predictor: `rf` (80 trees, depth 6), `xgb`
(80 rounds, depth 4), `mlp` (32-16, standardized).
Residual denominator is **K-fold CV MSE** on 上一批 (trees overfit
train residuals). Uniform stays the default; re-adjust only when
batches are clearly different (`gate_PO=1.25`, `gate_res=2`).

## Board (next-batch MSE)

### `rf`

| scene | uniform_pair | gated_pair | switch | resid | resid_po | oracle | uniform_hop | dre_hop |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `similar` | 0.582 | 0.582 | 0.555 | 0.540 | 0.540 | 0.540 | 0.725 | 0.725 |
| `covariate` | 0.762 | 0.758 | 0.784 | 0.888 | 0.958 | 0.728 | 1.206 | 2.489 |
| `concept` | 1.593 | 1.593 | 1.573 | 1.332 | 1.324 | 1.379 | 1.362 | 1.356 |
| `mixed` | 1.916 | 1.916 | 1.822 | 1.701 | 1.709 | 1.765 | 1.706 | 2.323 |

Fire rates (`rf`):

| scene | gated_pair | switch | resid | resid_po |
|---|---:|---:|---:|---:|
| `similar` | 0.00 | 0.12 | 0.00 | 0.00 |
| `covariate` | 0.08 | 0.21 | 0.17 | 0.17 |
| `concept` | 0.00 | 0.08 | 0.17 | 0.17 |
| `mixed` | 0.00 | 0.12 | 0.25 | 0.25 |

Concept hop path (`rf`, cut at `B_4`):

| t (train) | test | uniform_pair | switch | resid | resid_po | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.755 | 0.557 | 0.557 | 0.557 | 0.557 |
| 2 | `B_3` | 0.561 | 0.557 | 0.557 | 0.557 | 0.557 |
| 3 | `B_4` ← cut | 5.028 | 5.054 | 5.054 | 5.054 | 5.054 |
| 4 | `B_5` | 2.148 | 2.149 | 0.755 | 0.708 | 0.755 |
| 5 | `B_6` | 0.495 | 0.556 | 0.503 | 0.503 | 0.619 |
| 6 | `B_7` | 0.570 | 0.568 | 0.568 | 0.568 | 0.731 |

### `xgb`

| scene | uniform_pair | gated_pair | switch | resid | resid_po | oracle | uniform_hop | dre_hop |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `similar` | 0.426 | 0.443 | 0.415 | 0.396 | 0.396 | 0.396 | 0.556 | 0.562 |
| `covariate` | 0.545 | 0.564 | 0.597 | 0.614 | 0.622 | 0.511 | 0.842 | 1.117 |
| `concept` | 1.669 | 1.614 | 1.654 | 1.339 | 1.342 | 1.394 | 1.339 | 1.339 |
| `mixed` | 1.998 | 1.956 | 1.811 | 1.697 | 1.671 | 1.749 | 1.697 | 1.760 |

Fire rates (`xgb`):

| scene | gated_pair | switch | resid | resid_po |
|---|---:|---:|---:|---:|
| `similar` | 0.42 | 0.21 | 0.00 | 0.00 |
| `covariate` | 0.58 | 0.33 | 0.17 | 0.17 |
| `concept` | 0.29 | 0.12 | 0.17 | 0.17 |
| `mixed` | 0.46 | 0.25 | 0.25 | 0.25 |

Concept hop path (`xgb`, cut at `B_4`):

| t (train) | test | uniform_pair | switch | resid | resid_po | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.588 | 0.407 | 0.407 | 0.407 | 0.407 |
| 2 | `B_3` | 0.411 | 0.443 | 0.411 | 0.411 | 0.411 |
| 3 | `B_4` ← cut | 5.917 | 5.917 | 5.917 | 5.917 | 5.917 |
| 4 | `B_5` | 2.334 | 2.334 | 0.535 | 0.555 | 0.535 |
| 5 | `B_6` | 0.374 | 0.434 | 0.374 | 0.374 | 0.486 |
| 6 | `B_7` | 0.391 | 0.391 | 0.391 | 0.391 | 0.609 |

### `mlp`

| scene | uniform_pair | gated_pair | switch | resid | resid_po | oracle | uniform_hop | dre_hop |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `similar` | 0.448 | 0.457 | 0.437 | 0.426 | 0.426 | 0.426 | 0.600 | 0.605 |
| `covariate` | 0.773 | 0.837 | 1.065 | 1.096 | 1.142 | 0.674 | 1.512 | 1.857 |
| `concept` | 2.119 | 2.063 | 1.983 | 1.611 | 1.620 | 1.694 | 1.754 | 1.754 |
| `mixed` | 2.521 | 2.597 | 2.371 | 2.157 | 2.208 | 2.241 | 2.460 | 2.527 |

Fire rates (`mlp`):

| scene | gated_pair | switch | resid | resid_po |
|---|---:|---:|---:|---:|
| `similar` | 0.25 | 0.17 | 0.00 | 0.00 |
| `covariate` | 0.67 | 0.42 | 0.21 | 0.21 |
| `concept` | 0.29 | 0.25 | 0.17 | 0.17 |
| `mixed` | 0.54 | 0.42 | 0.33 | 0.33 |

Concept hop path (`mlp`, cut at `B_4`):

| t (train) | test | uniform_pair | switch | resid | resid_po | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.574 | 0.465 | 0.442 | 0.442 | 0.442 |
| 2 | `B_3` | 0.481 | 0.481 | 0.481 | 0.481 | 0.481 |
| 3 | `B_4` ← cut | 7.211 | 7.060 | 7.211 | 7.211 | 7.211 |
| 4 | `B_5` | 3.618 | 2.946 | 0.703 | 0.756 | 0.703 |
| 5 | `B_6` | 0.413 | 0.529 | 0.413 | 0.413 | 0.622 |
| 6 | `B_7` | 0.418 | 0.418 | 0.418 | 0.418 | 0.708 |

- **similar / covariate**: `P(Y|X)` stable → uniform slightly better.
- **resid**: residual hop-gate drops the old batch; should track oracle
  after the cut is observed.
- **resid_po**: same train-set switch plus √PO. Often a wash vs resid.
- **gated_pair / switch**: PO-ratio gate; misses a global map flip.
- **dre_hop**: X-only; overreacts to covariate hops.
