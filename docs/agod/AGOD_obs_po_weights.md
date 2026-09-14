# Observation-level PO-risk hard-reweight (iterated)

Default = **uniform**. OnlineRFPerm reject → obs PO on OOD batch →
`w ∝ PO / √PO / ∛PO / quantile / hybrid`. Sig-only next-MSE + hard-rank.

## Sig-only next MSE (↓ better)

| dataset | uniform | g_prop | g_sqrt | g_cbrt | g_quantile | g_hybrid | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 835439.5663 | 1185248.0535 | 917613.6149 | 878576.0213 | 913325.1359 | 1024913.3308 | `uniform` |
| `beijing_pm25` | 1672.8836 | 2039.2833 | 1752.2192 | 1621.3274 | 1701.3844 | 1963.7620 | `gated_cbrt` |
| `stocks_AAPL` | 0.0007 | 0.0008 | 0.0008 | 0.0007 | 0.0007 | 0.0008 | `uniform` |
| `stocks_MSFT` | 0.0003 | 0.0004 | 0.0003 | 0.0003 | 0.0003 | 0.0004 | `uniform` |
| `stocks_IWM` | 0.0005 | 0.0006 | 0.0006 | 0.0006 | 0.0006 | 0.0006 | `uniform` |
| `waymo_proxy` | 0.0093 | 0.0107 | 0.0099 | 0.0098 | 0.0101 | 0.0108 | `uniform` |

**Wins:** `uniform`=5, `gated_prop`=0, `gated_sqrt`=0, `gated_cbrt`=1, `gated_quantile`=0, `gated_hybrid`=0

## Hard-rank on reject batches (PO score vs oracle residual)

| dataset | spearman √PO / ∛PO / quantile / hybrid | P@20% √PO / ∛PO / quantile / hybrid |
|---|---|---|
| `metro_interstate` | 0.5335 / 0.5335 / 0.5335 / 0.5335 | 0.5441 / 0.5441 / 0.5441 / 0.5441 |
| `beijing_pm25` | 0.5521 / 0.5521 / 0.5521 / 0.5521 | 0.5294 / 0.5294 / 0.5294 / 0.5294 |
| `stocks_AAPL` | 0.7701 / 0.7701 / 0.7701 / 0.7701 | 0.6634 / 0.6634 / 0.6634 / 0.6634 |
| `stocks_MSFT` | 0.8126 / 0.8126 / 0.8126 / 0.8126 | 0.7602 / 0.7602 / 0.7602 / 0.7602 |
| `stocks_IWM` | 0.7769 / 0.7769 / 0.7769 / 0.7769 | 0.6716 / 0.6716 / 0.6716 / 0.6716 |
| `waymo_proxy` | 0.6492 / 0.6492 / 0.6492 / 0.6492 | 0.6180 / 0.6180 / 0.6180 / 0.6180 |

### Takeaway

1. **Obs-level PO ranks hard rows well** — Spearman vs oracle residual ≈ 0.53–0.81
   on reject batches (stocks highest). Monotone maps (√ / ∛ / quantile) share that ranking.
2. **IPTW into next-MSE is delicate** — `prop` overshoots; soft `cbrt` is safest and
   wins on beijing; many packs still prefer uniform for MSE.
3. Use PO weights as **hard-reweight after RFPerm**, not as an image-OOD detector.

Image-OOD stays with Mahalanobis / gradient / RF-binary.
