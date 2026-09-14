# Observation-level PO hard-reweight v3 (drift-adaptive)

## Logic (do we buy it?)

Yes, with a split:

1. **Hard-rank always** — obs PO identifies hard rows (Spearman ~0.5–0.8).
   That alone justifies PO as a *hardness score* for reweight targeting.
2. **Pack MSE only under shift** — when hard-tail = drift signal (beijing-like),
   soft temper helps next-MSE; when hard ≈ noise (calm stocks), uniform wins.
3. Therefore **λ should track drift intensity**, not a fixed temper.

v3: `λ = adaptive_temper(drift_intensity(PO, p, T))` with gate at drift≈0.2.
Also report **hard-subset next-MSE** (top-20% of next batch by PO).

## Drift intensity (mean on reject batches)

| dataset | drift_mean (qrt_adapt) | lam_mean |
|---|---:|---:|
| `metro_interstate` | 0.4837 | 0.266 |
| `beijing_pm25` | 0.7479 | 0.5136 |
| `stocks_AAPL` | 0.3228 | 0.1212 |
| `stocks_MSFT` | 0.6109 | 0.3853 |
| `stocks_IWM` | 0.5757 | 0.3522 |
| `waymo_proxy` | 0.2605 | 0.05668 |

## Sig-only next MSE — all rows (↓)

| dataset | uniform | qrtλ.5 | qrt_adapt | qrt_ad+r2 | cbrt_adapt | best |
|---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 8.354e+05 | 8.328e+05 | 8.383e+05 | 8.366e+05 | 8.344e+05 | `qrtλ.5` |
| `beijing_pm25` | 1673 | 1536 | 1653 | 1498 | 1625 | `qrt_ad+r2` |
| `stocks_AAPL` | 0.0006964 | 0.000702 | 0.0007022 | 0.0006953 | 0.0007022 | `qrt_ad+r2` |
| `stocks_MSFT` | 0.0003237 | 0.0003307 | 0.0003246 | 0.0003366 | 0.0003158 | `cbrt_adapt` |
| `stocks_IWM` | 0.0005463 | 0.0005551 | 0.0005558 | 0.0005462 | 0.0005621 | `qrt_ad+r2` |
| `waymo_proxy` | 0.009282 | 0.009368 | 0.00928 | 0.009269 | 0.009277 | `qrt_ad+r2` |

**Wins (all-row MSE):** `uniform`=0, `qrtλ.5`=1, `qrt_adapt`=0, `qrt_ad+r2`=4, `cbrt_adapt`=1

## Sig-only next MSE — hard top-20% of next batch (↓)  ← primary claim

| dataset | uniform | qrtλ.5 | qrt_adapt | qrt_ad+r2 | cbrt_adapt | best |
|---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 1.805e+06 | 1.755e+06 | 1.775e+06 | 1.871e+06 | 1.767e+06 | `qrtλ.5` |
| `beijing_pm25` | 5565 | 5017 | 5381 | 4959 | 5166 | `qrt_ad+r2` |
| `stocks_AAPL` | 0.002287 | 0.002268 | 0.002276 | 0.002321 | 0.00229 | `qrtλ.5` |
| `stocks_MSFT` | 0.001131 | 0.001132 | 0.001119 | 0.001172 | 0.001067 | `cbrt_adapt` |
| `stocks_IWM` | 0.001947 | 0.001951 | 0.001958 | 0.001954 | 0.002011 | `uniform` |
| `waymo_proxy` | 0.02424 | 0.02398 | 0.02429 | 0.0245 | 0.0242 | `qrtλ.5` |

**Wins (hard-subset MSE):** `uniform`=1, `qrtλ.5`=3, `qrt_adapt`=0, `qrt_ad+r2`=1, `cbrt_adapt`=1

## Rel. all-row MSE vs uniform (qrt_adapt / qrtλ.5)

| dataset | qrt_adapt | qrtλ.5 | drift |
|---|---:|---:|---:|
| `metro_interstate` | +0.3% | -0.3% | 0.4837 |
| `beijing_pm25` | -1.2% | -8.2% | 0.7479 |
| `stocks_AAPL` | +0.8% | +0.8% | 0.3228 |
| `stocks_MSFT` | +0.3% | +2.2% | 0.6109 |
| `stocks_IWM` | +1.7% | +1.6% | 0.5757 |
| `waymo_proxy` | -0.0% | +0.9% | 0.2605 |

## Hard-rank (unchanged claim)

| dataset | spearman | P@20% | n_reject |
|---|---:|---:|---:|
| `metro_interstate` | 0.5335 | 0.5441 | 8 |
| `beijing_pm25` | 0.5521 | 0.5294 | 6 |
| `stocks_AAPL` | 0.7701 | 0.6634 | 6 |
| `stocks_MSFT` | 0.8126 | 0.7602 | 4 |
| `stocks_IWM` | 0.7769 | 0.6716 | 4 |
| `waymo_proxy` | 0.6492 | 0.618 | 31 |

### Takeaway

- **Agree:** hard-rank is the robust justification; pack MSE is conditional on drift.
- **v3 test:** adaptive λ should fire on high-drift packs and stay near 0 on calm ones.
- Prefer hard-subset next-MSE when claiming hard-reweight benefit.

See `docs/agod/AGOD_obs_po_weights.md`.
