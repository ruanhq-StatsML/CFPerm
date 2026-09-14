# Gated online-RF √po_risk0 vs PO-tail localization

Probe is the shallow IPTW RF (`n_estimators=20`, `max_depth=4`)
on 上一批 as T=0. Instance `po_risk0` is `|Y−μ0(X)|` mixed with the
batch gap (`instance_po_risk`, mix=0.5). The next-batch model is
`rf` / `xgb` / `mlp` (not Ridge). Gate γ=1.5:

`e_now = err(μ0 fit B_{t-1} → B_t)`, `e_prev = err(μ0 fit B_{t-2} → B_{t-1})`

Fire only if `e_now / e_prev ≥ γ` (skip first hop). In-sample
`e1>e0` would fire every hop on trees — do not use it.

- **uniform_pair**: last two batches, w=1 (default).
- **rfperm**: on fire, train current batch with `w=√po_risk0`.
- **local**: same gate; train the high-`po_risk0` tail (`q=0.3`).
- **resid**: consecutive residual hop-gate (full learner MSE).
- **oracle**: knows the concept cut (train-set upper bound).

Always-on DRE / always-on √PO are off this board.

## Board (next-batch MSE)

### `rf`

| scene | uniform_pair | rfperm | local | resid | oracle |
|---|---:|---:|---:|---:|---:|
| `similar` | 0.540 | 0.540 | 0.540 | 0.540 | 0.540 |
| `covariate` | 0.728 | 0.845 | 1.087 | 0.852 | 0.728 |
| `concept` | 1.565 | 1.331 | 1.542 | 1.332 | 1.379 |
| `mixed` | 1.889 | 1.672 | 2.126 | 1.657 | 1.765 |

Fire rates (`rf`):

| scene | rfperm | local | resid |
|---|---:|---:|---:|
| `similar` | 0.00 | 0.00 | 0.00 |
| `covariate` | 0.04 | 0.04 | 0.12 |
| `concept` | 0.17 | 0.17 | 0.17 |
| `mixed` | 0.17 | 0.17 | 0.17 |

Concept hop path (`rf`, cut at `B_4`):

| t (train) | test | uniform_pair | rfperm | local | resid | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.557 | 0.557 | 0.557 | 0.557 | 0.557 |
| 2 | `B_3` | 0.557 | 0.557 | 0.557 | 0.557 | 0.557 |
| 3 | `B_4` ← cut | 5.054 | 5.054 | 5.054 | 5.054 | 5.054 |
| 4 | `B_5` | 2.149 | 0.748 | 2.017 | 0.755 | 0.755 |
| 5 | `B_6` | 0.503 | 0.503 | 0.503 | 0.503 | 0.619 |
| 6 | `B_7` | 0.568 | 0.568 | 0.568 | 0.568 | 0.731 |

### `xgb`

| scene | uniform_pair | rfperm | local | resid | oracle |
|---|---:|---:|---:|---:|---:|
| `similar` | 0.396 | 0.396 | 0.396 | 0.396 | 0.396 |
| `covariate` | 0.511 | 0.565 | 0.897 | 0.557 | 0.511 |
| `concept` | 1.639 | 1.339 | 1.551 | 1.339 | 1.394 |
| `mixed` | 1.957 | 1.648 | 1.849 | 1.661 | 1.749 |

Fire rates (`xgb`):

| scene | rfperm | local | resid |
|---|---:|---:|---:|
| `similar` | 0.00 | 0.00 | 0.00 |
| `covariate` | 0.04 | 0.04 | 0.04 |
| `concept` | 0.17 | 0.17 | 0.17 |
| `mixed` | 0.17 | 0.17 | 0.17 |

Concept hop path (`xgb`, cut at `B_4`):

| t (train) | test | uniform_pair | rfperm | local | resid | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.407 | 0.407 | 0.407 | 0.407 | 0.407 |
| 2 | `B_3` | 0.411 | 0.411 | 0.411 | 0.411 | 0.411 |
| 3 | `B_4` ← cut | 5.917 | 5.917 | 5.917 | 5.917 | 5.917 |
| 4 | `B_5` | 2.334 | 0.536 | 1.807 | 0.535 | 0.535 |
| 5 | `B_6` | 0.374 | 0.374 | 0.374 | 0.374 | 0.486 |
| 6 | `B_7` | 0.391 | 0.391 | 0.391 | 0.391 | 0.609 |

- **similar / covariate**: `P(Y|X)` stable → uniform; rfperm/local
  should stay quiet.
- **concept**: fire at the cut hop, then drop the old batch.
- **local** vs **rfperm**: same gate; localization keeps the high-PO
  tail instead of IPTW on the whole new batch.
