# Observation-level PO hard-reweight v2 (temper / soft / top-k)

Default = **uniform**. On OnlineRFPerm reject, obs PO → soft weights.
v2 adds PO^{1/4}, temper(λ=0.5), top-20% boost, soft p-gate.

## Sig-only next MSE (↓ better)

| dataset | uniform | cbrt | qrt | cbrtλ.5 | qrtλ.5 | topk20 | cbrt_soft | best |
|---|---:---|---:---|---:---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | 835439.5663 | 878576.0213 | 854390.7049 | 838005.4775 | 832847.9943 | 856440.5587 | 840610.7980 | `qrtλ.5` |
| `beijing_pm25` | 1672.8836 | 1621.3274 | 1643.6229 | 1623.2712 | 1536.1622 | 1626.4155 | 1632.1986 | `qrtλ.5` |
| `stocks_AAPL` | 0.0007 | 0.0007 | 0.0007 | 0.0007 | 0.0007 | 0.0007 | 0.0007 | `uniform` |
| `stocks_MSFT` | 0.0003 | 0.0003 | 0.0003 | 0.0003 | 0.0003 | 0.0003 | 0.0003 | `uniform` |
| `stocks_IWM` | 0.0005 | 0.0006 | 0.0006 | 0.0006 | 0.0006 | 0.0006 | 0.0006 | `uniform` |
| `waymo_proxy` | 0.0093 | 0.0098 | 0.0095 | 0.0095 | 0.0094 | 0.0098 | 0.0095 | `uniform` |

**Wins:** `uniform`=4, `cbrt`=0, `qrt`=0, `cbrtλ.5`=0, `qrtλ.5`=2, `topk20`=0, `cbrt_soft`=0

## Hard-rank (PO vs oracle residual) on reject batches

| dataset | spearman | P@20% | n_reject |
|---|---:|---:|---:|
| `metro_interstate` | 0.5335 | 0.5441 | 8 |
| `beijing_pm25` | 0.5521 | 0.5294 | 6 |
| `stocks_AAPL` | 0.7701 | 0.6634 | 6 |
| `stocks_MSFT` | 0.8126 | 0.7602 | 4 |
| `stocks_IWM` | 0.7769 | 0.6716 | 4 |
| `waymo_proxy` | 0.6492 | 0.6180 | 31 |

### Takeaway

1. Obs PO still **ranks hard rows** (use for hard-reweight targeting).
2. Prefer **tempered / soft** maps (∛·λ0.5, ¼, top-k) over raw prop/√ for MSE.
3. Not an image-OOD detector — Mahalanobis / gradient / RF-binary for that.

See `docs/agod/AGOD_obs_po_weights.md`.
