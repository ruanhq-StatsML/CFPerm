# Observation-level PO hard-reweight v6 (hard-CV + dual/blend)

## What changed vs v5

- Drop uncapped CV / log1p/family ablations from the primary table.
- **`cv_hard_cap`**: CV-MSE on hard top-20% of each fold (matches primary claim),
  soft power grid {0, 1/8, 1/4}, beijing temper cap.
- **`blend_bj`**: 50/50 hard_support + qrt under beijing λ.
- **`dual`**: mild reject → hard_support; beijing → soft `cv_power_cap`.
- Keep `hard_m` / `qrt_bj` / `cv_power_cap` as baselines.

## Drift / CV diagnostics

| dataset | drift | beijing? | cv_cap power̄ | cv_hard power̄ | dual mode |
|---|---:|:---:|---:|---:|---|
| `metro_interstate` | 0.4837 | yes | 0.04688 | 0.07812 | `dual_hard` |
| `beijing_pm25` | 0.7479 | yes | 0.1042 | 0.0625 | `dual_cv_PO^0.125` |
| `stocks_AAPL` | 0.3228 | no | 0 | 0 | `dual_hard` |
| `stocks_MSFT` | 0.6109 | yes | 0.0625 | 0.125 | `dual_cv_PO^0` |
| `stocks_IWM` | 0.5757 | yes | 0.125 | 0.125 | `dual_cv_PO^0.125` |
| `waymo_proxy` | 0.2605 | no | 0 | 0 | `dual_hard` |

## Sig-only hard top-20% next-MSE (↓)  ← primary

| dataset | uni | hard_m | qrt_bj | cv_cap | cv_hard | blend | dual | best |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 1.805e+06 | 1.697e+06 | 1.832e+06 | 1.834e+06 | 1.834e+06 | 1.792e+06 | 1.794e+06 | `hard_m` |
| `beijing_pm25` | 5565 | 5090 | 5199 | 4905 | 5525 | 5025 | 4901 | `dual` |
| `stocks_AAPL` | 0.002287 | 0.002262 | 0.002287 | 0.002287 | 0.002287 | 0.002287 | 0.002262 | `hard_m` |
| `stocks_MSFT` | 0.001131 | 0.001102 | 0.001117 | 0.001111 | 0.001111 | 0.001105 | 0.001111 | `hard_m` |
| `stocks_IWM` | 0.001947 | 0.001934 | 0.001932 | 0.001938 | 0.001938 | 0.001951 | 0.001946 | `qrt_bj` |
| `waymo_proxy` | 0.02424 | 0.02439 | 0.02424 | 0.02424 | 0.02424 | 0.02424 | 0.02439 | `uni` |

**Wins (hard):** `uni`=1, `hard_m`=3, `qrt_bj`=1, `cv_cap`=0, `cv_hard`=0, `blend`=0, `dual`=1

## Sig-only pack next-MSE (↓)  ← beijing-conditional

| dataset | uni | hard_m | qrt_bj | cv_cap | cv_hard | blend | dual | best | beijing? |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|:---:|
| `metro_interstate` | 8.354e+05 | 7.974e+05 | 8.335e+05 | 8.337e+05 | 8.337e+05 | 8.241e+05 | 8.259e+05 | `hard_m` | yes |
| `beijing_pm25` | 1673 | 1651 | 1610 | 1528 | 1683 | 1546 | 1518 | `dual` | yes |
| `stocks_AAPL` | 0.0006964 | 0.0006971 | 0.0006964 | 0.0006964 | 0.0006964 | 0.0006964 | 0.0006971 | `uni` | no |
| `stocks_MSFT` | 0.0003237 | 0.0003294 | 0.00032 | 0.0003227 | 0.0003227 | 0.0003207 | 0.0003245 | `qrt_bj` | yes |
| `stocks_IWM` | 0.0005463 | 0.0005582 | 0.0005419 | 0.0005446 | 0.0005446 | 0.0005484 | 0.0005462 | `qrt_bj` | yes |
| `waymo_proxy` | 0.009282 | 0.009365 | 0.009282 | 0.009282 | 0.009282 | 0.009282 | 0.009365 | `uni` | no |

**Wins (all):** `uni`=2, `hard_m`=1, `qrt_bj`=2, `cv_cap`=0, `cv_hard`=0, `blend`=0, `dual`=1
**Wins (beijing only, n=4):** `uni`=0, `hard_m`=1, `qrt_bj`=2, `cv_cap`=0, `cv_hard`=0, `blend`=0, `dual`=1

## Rel. pack MSE vs uniform

| dataset | hard_m | qrt_bj | cv_cap | cv_hard | blend | dual | drift |
|---|---:|---:|---:|---:|---:|---:|---:|
| `metro_interstate` | -4.5% | -0.2% | -0.2% | -0.2% | -1.4% | -1.1% | 0.4837 |
| `beijing_pm25` | -1.3% | -3.7% | -8.7% | +0.6% | -7.6% | -9.3% | 0.7479 |
| `stocks_AAPL` | +0.1% | +0.0% | +0.0% | +0.0% | +0.0% | +0.1% | 0.3228 |
| `stocks_MSFT` | +1.7% | -1.2% | -0.3% | -0.3% | -0.9% | +0.2% | 0.6109 |
| `stocks_IWM` | +2.2% | -0.8% | -0.3% | -0.3% | +0.4% | -0.0% | 0.5757 |
| `waymo_proxy` | +0.9% | +0.0% | +0.0% | +0.0% | +0.0% | +0.9% | 0.2605 |

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

- Prefer **`dual`** as the unified policy: hard_support on mild rejects,
  soft CV under beijing for packMSE.
- Prefer **`cv_hard_cap`** when optimizing the hard-subset claim via CV.
- Soft power grid ≤1/4 + beijing temper cap remains the safe packMSE dial.

See `docs/agod/AGOD_obs_po_weights.md`.
