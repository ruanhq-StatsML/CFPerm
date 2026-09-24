# PO-risk adaptation efficiency

rank_eff wins {'refit': 6}; mse_eff wins {'refit': 1, 'probe': 5} — ranking skill ≠ MSE win; FLOPs decide which PO to keep.

- batches=40 · batch_size=100 · n_control=1 · datasets=6

## Per dataset

| dataset | n_rej | ref ρ | probe Δρ | refit Δρ | probe rank_eff | refit rank_eff | probe mse_eff | refit mse_eff | reading |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 7 | 0.3021 | 0.0134 | 0.0176 | 0.0187 | 0.1396 | 1.06e+05 | 8.192e+05 | refit: more Spearman per reject-FLOP than probe's always-on cost |
| `beijing_pm25` | 5 | 0.0844 | 0.3087 | 0.3149 | 0.4287 | 3.4986 | -382.4970 | -3488 | refit: more Spearman per reject-FLOP than probe's always-on cost |
| `stocks_AAPL` | 7 | 0.2865 | 0.3839 | 0.3789 | 0.5332 | 3.0073 | -0.0001382 | -0.0009838 | refit: more Spearman per reject-FLOP than probe's always-on cost |
| `waymo_proxy` | 22 | 0.1149 | 0.4345 | 0.4143 | 0.6035 | 1.0462 | -0.001571 | -0.001992 | refit: more Spearman per reject-FLOP than probe's always-on cost |
| `stocks_MSFT` | 9 | 0.3778 | 0.3755 | 0.3436 | 0.5215 | 2.1209 | -0.0001187 | -0.0005367 | refit: more Spearman per reject-FLOP than probe's always-on cost |
| `stocks_IWM` | 7 | 0.6084 | 0.1169 | 0.0802 | 0.1624 | 0.6369 | -7.459e-05 | -0.0002671 | refit: more Spearman per reject-FLOP than probe's always-on cost |

## How to read

- **rank_eff**: ΔSpearman vs frozen ref per million adaptation FLOPs.
- **mse_eff**: (uniform − mode) sig-only MSE per million FLOPs (positive = cheaper error drop).
- Probe pays *always-on* fits; refit pays only on reject — so a tiny Δρ for refit can still beat probe on rank_eff.
- Negative mse_eff with positive rank_eff ⇒ keep for triage, not for IPTW.

