# Gated online-RF √po_risk0 reweighting

Overview: `docs/agod/AGOD_overview.md` (anneal to uniform;
reweight only on a consecutive OOS jump).
LaTeX: `docs/agod/AGOD_po_refit_gated_tables_only.tex`,
`docs/agod/AGOD_performance_tables.tex`.

Probe is the shallow IPTW RF (`n_estimators=20`, `max_depth=4`)
on 上一批 as T=0. Instance `po_risk0` is `|Y−μ0(X)|` mixed with the
batch gap (`instance_po_risk`, mix=0.5). The next-batch model is
`rf` / `xgb` / `mlp` (not Ridge). Gate γ=1.5:

`e_now = err(μ0 fit B_{t-1} → B_t)`, `e_prev = err(μ0 fit B_{t-2} → B_{t-1})`

Fire only if `e_now / e_prev ≥ γ` (skip first hop). In-sample
`e1>e0` would fire every hop on trees — do not use it.

- **uniform_pair**: last two batches, w=1 (default).
- **dre**: same last-two rows; always-on logistic density-ratio
  on T=1 vs T=0 (X only).
- **rfperm**: same rows; on fire, T=0 stays 1 and T=1 gets
  `w=√po_risk0` (mean 1).
- **resid**: drop the old batch when consecutive residual MSE jumps.
- **oracle**: knows the concept cut (train-set upper bound).

Always-on √PO / PO-tail subset are off this board. DRE is the
p(x) baseline on the *same* last-two rows.

## Board (next-batch MSE)

### `rf`

| scene | uniform | DRE | rfperm | drop-old | oracle |
|---|---:|---:|---:|---:|---:|
| `similar` | 0.540 | 0.541 | 0.540 | 0.540 | 0.540 |
| `covariate` | 0.728 | 0.874 | 0.720 | 0.852 | 0.728 |
| `concept` | 1.565 | 1.584 | 1.466 | 1.332 | 1.379 |
| `mixed` | 1.889 | 2.020 | 1.818 | 1.657 | 1.765 |

Fire rates (`rf`):

| scene | rfperm | resid |
|---|---:|---:|
| `similar` | 0.00 | 0.00 |
| `covariate` | 0.04 | 0.12 |
| `concept` | 0.17 | 0.17 |
| `mixed` | 0.17 | 0.17 |

Concept hop path (`rf`, cut at `B_4`):

| t (train) | test | uniform | DRE | rfperm | drop-old | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.557 | 0.569 | 0.557 | 0.557 | 0.557 |
| 2 | `B_3` | 0.557 | 0.547 | 0.557 | 0.557 | 0.557 |
| 3 | `B_4` ← cut | 5.054 | 5.051 | 5.054 | 5.054 | 5.054 |
| 4 | `B_5` | 2.149 | 2.254 | 1.556 | 0.755 | 0.755 |
| 5 | `B_6` | 0.503 | 0.506 | 0.503 | 0.503 | 0.619 |
| 6 | `B_7` | 0.568 | 0.575 | 0.568 | 0.568 | 0.731 |

### `xgb`

| scene | uniform | DRE | rfperm | drop-old | oracle |
|---|---:|---:|---:|---:|---:|
| `similar` | 0.396 | 0.393 | 0.396 | 0.396 | 0.396 |
| `covariate` | 0.511 | 0.586 | 0.509 | 0.557 | 0.511 |
| `concept` | 1.639 | 1.651 | 1.592 | 1.339 | 1.394 |
| `mixed` | 1.957 | 2.132 | 1.920 | 1.661 | 1.749 |

Fire rates (`xgb`):

| scene | rfperm | resid |
|---|---:|---:|
| `similar` | 0.00 | 0.00 |
| `covariate` | 0.04 | 0.04 |
| `concept` | 0.17 | 0.17 |
| `mixed` | 0.17 | 0.17 |

Concept hop path (`xgb`, cut at `B_4`):

| t (train) | test | uniform | DRE | rfperm | drop-old | oracle |
|---|---|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.407 | 0.382 | 0.407 | 0.407 | 0.407 |
| 2 | `B_3` | 0.411 | 0.412 | 0.411 | 0.411 | 0.411 |
| 3 | `B_4` ← cut | 5.917 | 5.956 | 5.917 | 5.917 | 5.917 |
| 4 | `B_5` | 2.334 | 2.380 | 2.051 | 0.535 | 0.535 |
| 5 | `B_6` | 0.374 | 0.371 | 0.374 | 0.374 | 0.486 |
| 6 | `B_7` | 0.391 | 0.405 | 0.391 | 0.391 | 0.609 |

- **similar / covariate**: `P(Y|X)` stable → uniform; rfperm
  should stay quiet. DRE reweights p(x) and can hurt covariate.
- **concept**: P(X) stable, P(Y|X) flips → DRE ≈ uniform;
  rfperm upweights the new map. Drop-old is the harder lever.
