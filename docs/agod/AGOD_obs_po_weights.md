# Observation-level PO hard-reweight v7 (blend mix + dual_r2)

## What changed vs v6

- Keep **`dual`** (v6 beijing packMSE winner) and **`hard_m`**.
- **`dual_r2`**: same dual policy with `n_recent=2` for PO/μ0 windows.
- **Blend mix scan**: `blend_25` / `blend_50` / `blend_75` (hard_support vs qrt).
- Drop cv_hard_cap from primary table (did not transfer in v6).

## Drift / dual diagnostics

| dataset | drift | beijing? | dual mode | dual_r2 mode |
|---|---:|:---:|---|---|
| `metro_interstate` | 0.4837 | yes | `dual_hard` | `dual_hard` |
| `beijing_pm25` | 0.7479 | yes | `dual_cv_PO^0.125` | `dual_cv_PO^0` |
| `stocks_AAPL` | 0.3228 | no | `dual_hard` | `dual_hard` |
| `stocks_MSFT` | 0.6109 | yes | `dual_cv_PO^0` | `dual_cv_PO^0` |
| `stocks_IWM` | 0.5757 | yes | `dual_cv_PO^0.125` | `dual_cv_PO^0.125` |
| `waymo_proxy` | 0.2605 | no | `dual_hard` | `dual_hard` |

## Sig-only hard top-20% next-MSE (↓)  ← primary

| dataset | uni | hard_m | dual | dual_r2 | b25 | b50 | b75 | best |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 1.805e+06 | 1.697e+06 | 1.794e+06 | 1.872e+06 | 1.809e+06 | 1.792e+06 | 1.788e+06 | `hard_m` |
| `beijing_pm25` | 5565 | 5090 | 4901 | 5380 | 5062 | 5025 | 5420 | `dual` |
| `stocks_AAPL` | 0.002287 | 0.002262 | 0.002262 | 0.002341 | 0.002287 | 0.002287 | 0.002287 | `hard_m` |
| `stocks_MSFT` | 0.001131 | 0.001102 | 0.001111 | 0.001161 | 0.001128 | 0.001105 | 0.001155 | `hard_m` |
| `stocks_IWM` | 0.001947 | 0.001934 | 0.001946 | 0.001943 | 0.001933 | 0.001951 | 0.001966 | `b25` |
| `waymo_proxy` | 0.02424 | 0.02439 | 0.02439 | 0.02435 | 0.02424 | 0.02424 | 0.02424 | `uni` |

**Wins (hard):** `uni`=1, `hard_m`=3, `dual`=1, `dual_r2`=0, `b25`=1, `b50`=0, `b75`=0

## Sig-only pack next-MSE (↓)  ← beijing-conditional

| dataset | uni | hard_m | dual | dual_r2 | b25 | b50 | b75 | best | beijing? |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|:---:|
| `metro_interstate` | 8.354e+05 | 7.974e+05 | 8.259e+05 | 8.413e+05 | 8.262e+05 | 8.241e+05 | 8.314e+05 | `hard_m` | yes |
| `beijing_pm25` | 1673 | 1651 | 1518 | 1629 | 1570 | 1546 | 1705 | `dual` | yes |
| `stocks_AAPL` | 0.0006964 | 0.0006971 | 0.0006971 | 0.0007053 | 0.0006964 | 0.0006964 | 0.0006964 | `uni` | no |
| `stocks_MSFT` | 0.0003237 | 0.0003294 | 0.0003245 | 0.0003295 | 0.0003276 | 0.0003207 | 0.0003356 | `b50` | yes |
| `stocks_IWM` | 0.0005463 | 0.0005582 | 0.0005462 | 0.0005623 | 0.0005429 | 0.0005484 | 0.0005546 | `b25` | yes |
| `waymo_proxy` | 0.009282 | 0.009365 | 0.009365 | 0.009385 | 0.009282 | 0.009282 | 0.009282 | `uni` | no |

**Wins (all):** `uni`=2, `hard_m`=1, `dual`=1, `dual_r2`=0, `b25`=1, `b50`=1, `b75`=0
**Wins (beijing only, n=4):** `uni`=0, `hard_m`=1, `dual`=1, `dual_r2`=0, `b25`=1, `b50`=1, `b75`=0

## Rel. pack MSE vs uniform

| dataset | hard_m | dual | dual_r2 | b25 | b50 | b75 | drift |
|---|---:|---:|---:|---:|---:|---:|---:|
| `metro_interstate` | -4.5% | -1.1% | +0.7% | -1.1% | -1.4% | -0.5% | 0.4837 |
| `beijing_pm25` | -1.3% | -9.3% | -2.6% | -6.1% | -7.6% | +1.9% | 0.7479 |
| `stocks_AAPL` | +0.1% | +0.1% | +1.3% | +0.0% | +0.0% | +0.0% | 0.3228 |
| `stocks_MSFT` | +1.7% | +0.2% | +1.8% | +1.2% | -0.9% | +3.7% | 0.6109 |
| `stocks_IWM` | +2.2% | -0.0% | +2.9% | -0.6% | +0.4% | +1.5% | 0.5757 |
| `waymo_proxy` | +0.9% | +0.9% | +1.1% | +0.0% | +0.0% | +0.0% | 0.2605 |

## Hard-rank

| dataset | spearman | P@20% | n_reject |
|---|---:|---:|---:|
| `metro_interstate` | 0.5335 | 0.5441 | 8 |
| `beijing_pm25` | 0.5521 | 0.5294 | 6 |
| `stocks_AAPL` | 0.7701 | 0.6634 | 6 |
| `stocks_MSFT` | 0.8126 | 0.7602 | 4 |
| `stocks_IWM` | 0.7769 | 0.6716 | 4 |
| `waymo_proxy` | 0.6492 | 0.618 | 31 |

### Takeaway / locked recipe (v7)

- **`n_recent=1`**: `dual_r2` loses to `dual` on beijing (−2.6% vs −9.3%) and hurts calm packs → reject r2.
- **Blend mix ≤ 0.5**: `b50` (−7.6% beijing) > `b25` (−6.1%) ≫ `b75` (+1.9%). Prefer more qrt than hard in the blend.
- **Default unified policy remains `dual`** (mild→hard_support, BJ→soft CV, n_recent=1).
- **Hard-only claim** still uses `hard_m` (wins hard-subset 3/6).

See `docs/agod/AGOD_obs_po_weights.md`.
