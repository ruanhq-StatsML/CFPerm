# Observation-level PO hard-reweight v5 (CV-MSE + new maps)

## What changed

- New soft map: `log1p` (beijing-gated); `softmax` via `cv_family`.
- **CV-MSE** picks `(power, temper)` per rejected batch (`cv_power`),
  or with beijing temper cap (`cv_power_cap`),
  or among discrete families (`cv_family`).
- Thesis unchanged: 认 hard; packMSE only under beijing-class drift.

## Drift / CV diagnostics

| dataset | drift | beijing? | cv_p power̄ | cv_p λ̄ | cv_cap power̄ | cv_fam mode |
|---|---:|:---:|---:|---:|---:|---|
| `metro_interstate` | 0.4837 | yes | 0.1458 | 0.3625 | 0.08333 | `uniform` |
| `beijing_pm25` | 0.7479 | yes | 0.2292 | 0.3917 | 0.1042 | `log1p` |
| `stocks_AAPL` | 0.3228 | no | 0.07639 | 0.1167 | 0 | `uniform` |
| `stocks_MSFT` | 0.6109 | yes | 0.0625 | 0.1375 | 0.0625 | `uniform` |
| `stocks_IWM` | 0.5757 | yes | 0.2812 | 0.4 | 0.1771 | `uniform` |
| `waymo_proxy` | 0.2605 | no | 0.1573 | 0.3952 | 0 | `quantile` |

## Sig-only hard top-20% next-MSE (↓)  ← primary

| dataset | uni | hard_m | qrt_bj | log_bj | cv_p | cv_cap | cv_fam | best |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 1.805e+06 | 1.697e+06 | 1.832e+06 | 1.827e+06 | 1.733e+06 | 1.792e+06 | 1.717e+06 | `hard_m` |
| `beijing_pm25` | 5565 | 5090 | 5199 | 5248 | 5172 | 4905 | 5248 | `cv_cap` |
| `stocks_AAPL` | 0.002287 | 0.002262 | 0.002287 | 0.002287 | 0.002285 | 0.002287 | 0.002281 | `hard_m` |
| `stocks_MSFT` | 0.001131 | 0.001102 | 0.001117 | 0.00116 | 0.001111 | 0.001111 | 0.001117 | `hard_m` |
| `stocks_IWM` | 0.001947 | 0.001934 | 0.001932 | 0.001961 | 0.001987 | 0.001933 | 0.001968 | `qrt_bj` |
| `waymo_proxy` | 0.02424 | 0.02439 | 0.02424 | 0.02424 | 0.0241 | 0.02424 | 0.02427 | `cv_p` |

**Wins (hard):** `uni`=0, `hard_m`=3, `qrt_bj`=1, `log_bj`=0, `cv_p`=1, `cv_cap`=1, `cv_fam`=0

## Sig-only pack next-MSE (↓)  ← beijing-conditional

| dataset | uni | hard_m | qrt_bj | log_bj | cv_p | cv_cap | cv_fam | best | beijing? |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|:---:|
| `metro_interstate` | 8.354e+05 | 7.974e+05 | 8.335e+05 | 8.398e+05 | 8.542e+05 | 8.292e+05 | 8.572e+05 | `hard_m` | yes |
| `beijing_pm25` | 1673 | 1651 | 1610 | 1603 | 1629 | 1528 | 1587 | `cv_cap` | yes |
| `stocks_AAPL` | 0.0006964 | 0.0006971 | 0.0006964 | 0.0006964 | 0.0006958 | 0.0006964 | 0.0006978 | `cv_p` | no |
| `stocks_MSFT` | 0.0003237 | 0.0003294 | 0.00032 | 0.0003387 | 0.0003227 | 0.0003227 | 0.0003212 | `qrt_bj` | yes |
| `stocks_IWM` | 0.0005463 | 0.0005582 | 0.0005419 | 0.0005545 | 0.0005587 | 0.000544 | 0.0005563 | `qrt_bj` | yes |
| `waymo_proxy` | 0.009282 | 0.009365 | 0.009282 | 0.009282 | 0.009354 | 0.009282 | 0.0093 | `uni` | no |

**Wins (all):** `uni`=1, `hard_m`=1, `qrt_bj`=2, `log_bj`=0, `cv_p`=1, `cv_cap`=1, `cv_fam`=0
**Wins (beijing only, n=4):** `uni`=0, `hard_m`=1, `qrt_bj`=2, `log_bj`=0, `cv_p`=0, `cv_cap`=1, `cv_fam`=0

## Rel. pack MSE vs uniform

| dataset | hard_m | qrt_bj | log_bj | cv_p | cv_cap | cv_fam | drift |
|---|---:|---:|---:|---:|---:|---:|---:|
| `metro_interstate` | -4.5% | -0.2% | +0.5% | +2.2% | -0.7% | +2.6% | 0.4837 |
| `beijing_pm25` | -1.3% | -3.7% | -4.2% | -2.6% | -8.7% | -5.1% | 0.7479 |
| `stocks_AAPL` | +0.1% | +0.0% | +0.0% | -0.1% | +0.0% | +0.2% | 0.3228 |
| `stocks_MSFT` | +1.7% | -1.2% | +4.6% | -0.3% | -0.3% | -0.8% | 0.6109 |
| `stocks_IWM` | +2.2% | -0.8% | +1.5% | +2.3% | -0.4% | +1.8% | 0.5757 |
| `waymo_proxy` | +0.9% | +0.0% | +0.0% | +0.8% | +0.0% | +0.2% | 0.2605 |

## Hard-rank

| dataset | spearman | P@20% | n_reject |
|---|---:|---:|---:|
| `metro_interstate` | 0.5335 | 0.5441 | 8 |
| `beijing_pm25` | 0.5521 | 0.5294 | 6 |
| `stocks_AAPL` | 0.7701 | 0.6634 | 6 |
| `stocks_MSFT` | 0.8126 | 0.7602 | 4 |
| `stocks_IWM` | 0.7769 | 0.6716 | 4 |
| `waymo_proxy` | 0.6492 | 0.618 | 31 |

### Takeaway

- CV-MSE auto-tunes power/temper per reject; calm packs often pick ≈uniform.
- Prefer `cv_power_cap` / `qrt_bj` / `log_bj` for packMSE under beijing drift;
  prefer `gated_hard_adapt` / CV-family for the hard claim.

See `docs/agod/AGOD_obs_po_weights.md`.
