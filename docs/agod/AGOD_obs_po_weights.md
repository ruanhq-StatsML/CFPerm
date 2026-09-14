# Observation-level PO hard-reweight v4

## Do we buy hard / beijing-drift / packMSE?

**Yes, with a clean split (this is the locked thesis):**

1. **认 hard** — obs PO is a hardness score. Spearman ~0.5–0.8, P@20% ≫ random.
   The matching mechanism is **hard_support** (boost only the hard top-k),
   not diffuse soft IPTW. Primary metric = **hard-subset next-MSE**.
2. **条件认 beijing类漂移 → packMSE** — all-row pack MSE is only a fair
   claim when hard-tail ≈ shift signal (drift ≫ mild, gate≈0.45).
   On calm packs, hard ≈ noise → uniform wins packMSE; do not force lift.
3. **不认** chasing packMSE on every RFPerm reject, or treating PO as an
   image-OOD detector. Reject gate stays OnlineRFPerm; PO is post-hoc reweight.

v4: hard_support λ via mild gate (0.2) vs beijing gate (0.45);
soft qrt kept as the packMSE-oriented comparator under the same gates.

## Drift intensity (mean on reject batches)

| dataset | drift_mean | beijing_class? | lam qrt_mild | lam qrt_bj | lam hard_mild |
|---|---:|:---:|---:|---:|---:|
| `metro_interstate` | 0.4837 | yes | 0.266 | 0.08085 | 0.266 |
| `beijing_pm25` | 0.7479 | yes | 0.5136 | 0.4298 | 0.5136 |
| `stocks_AAPL` | 0.3228 | no | 0.1212 | 0 | 0.1212 |
| `stocks_MSFT` | 0.6109 | yes | 0.3853 | 0.2314 | 0.3853 |
| `stocks_IWM` | 0.5757 | yes | 0.3522 | 0.1916 | 0.3522 |
| `waymo_proxy` | 0.2605 | no | 0.05668 | 0 | 0.05668 |

## Sig-only next MSE — hard top-20% (↓)  ← primary claim (认 hard)

| dataset | uniform | qrt_mild | qrt_bj | hard_mild | hard_bj | best |
|---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 1.805e+06 | 1.775e+06 | 1.832e+06 | 1.697e+06 | 1.793e+06 | `hard_mild` |
| `beijing_pm25` | 5565 | 5381 | 5199 | 5090 | 4811 | `hard_bj` |
| `stocks_AAPL` | 0.002287 | 0.002276 | 0.002287 | 0.002262 | 0.002287 | `hard_mild` |
| `stocks_MSFT` | 0.001131 | 0.001119 | 0.001117 | 0.001102 | 0.001189 | `hard_mild` |
| `stocks_IWM` | 0.001947 | 0.001958 | 0.001932 | 0.001934 | 0.001964 | `qrt_bj` |
| `waymo_proxy` | 0.02424 | 0.02429 | 0.02424 | 0.02439 | 0.02424 | `uniform` |

**Wins (hard-subset):** `uniform`=1, `qrt_mild`=0, `qrt_bj`=1, `hard_mild`=3, `hard_bj`=1

## Sig-only next MSE — all rows / pack (↓)  ← only claim under beijing drift

| dataset | uniform | qrt_mild | qrt_bj | hard_mild | hard_bj | best | beijing? |
|---|---:---|---:---|---:---|---:---|---:|---|:---:|
| `metro_interstate` | 8.354e+05 | 8.383e+05 | 8.335e+05 | 7.974e+05 | 8.322e+05 | `hard_mild` | yes |
| `beijing_pm25` | 1673 | 1653 | 1610 | 1651 | 1533 | `hard_bj` | yes |
| `stocks_AAPL` | 0.0006964 | 0.0007022 | 0.0006964 | 0.0006971 | 0.0006964 | `uniform` | no |
| `stocks_MSFT` | 0.0003237 | 0.0003246 | 0.00032 | 0.0003294 | 0.0003416 | `qrt_bj` | yes |
| `stocks_IWM` | 0.0005463 | 0.0005558 | 0.0005419 | 0.0005582 | 0.0005516 | `qrt_bj` | yes |
| `waymo_proxy` | 0.009282 | 0.00928 | 0.009282 | 0.009365 | 0.009282 | `qrt_mild` | no |

**Wins (all packs):** `uniform`=1, `qrt_mild`=1, `qrt_bj`=2, `hard_mild`=1, `hard_bj`=1
**Wins among beijing-class packs only (n=4):** `uniform`=0, `qrt_mild`=0, `qrt_bj`=2, `hard_mild`=1, `hard_bj`=1

## Rel. pack MSE vs uniform (soft qrt paths)

| dataset | qrt_mild | qrt_bj | hard_mild | hard_bj | drift | beijing? |
|---|---:|---:|---:|---:|---:|:---:|
| `metro_interstate` | +0.3% | -0.2% | -4.5% | -0.4% | 0.4837 | yes |
| `beijing_pm25` | -1.2% | -3.7% | -1.3% | -8.4% | 0.7479 | yes |
| `stocks_AAPL` | +0.8% | +0.0% | +0.1% | +0.0% | 0.3228 | no |
| `stocks_MSFT` | +0.3% | -1.2% | +1.7% | +5.5% | 0.6109 | yes |
| `stocks_IWM` | +1.7% | -0.8% | +2.2% | +1.0% | 0.5757 | yes |
| `waymo_proxy` | -0.0% | +0.0% | +0.9% | +0.0% | 0.2605 | no |

## Hard-rank (认 hard — unchanged)

| dataset | spearman | P@20% | n_reject |
|---|---:|---:|---:|
| `metro_interstate` | 0.5335 | 0.5441 | 8 |
| `beijing_pm25` | 0.5521 | 0.5294 | 6 |
| `stocks_AAPL` | 0.7701 | 0.6634 | 6 |
| `stocks_MSFT` | 0.8126 | 0.7602 | 4 |
| `stocks_IWM` | 0.7769 | 0.6716 | 4 |
| `waymo_proxy` | 0.6492 | 0.618 | 31 |

### Takeaway

- **认 hard**: use PO to find / boost hard support; judge by hard-subset MSE + rank.
- **条件认 packMSE**: only advertise all-row lift on beijing-class drift packs;
  prefer `qrt_bj` / `hard_bj` (high gate) so calm rejects stay near uniform.
- Mild-gate soft qrt remains a useful ablation, not the default packMSE claim.

See `docs/agod/AGOD_obs_po_weights.md`.
