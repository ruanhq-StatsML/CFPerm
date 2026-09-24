# PO-risk adaptation efficiency

mean duty=0.244 (breakeven=1.00) ⇒ E[refit]/E[probe]≈0.244; 6/6 packs prefer refit on budget; rank_eff wins {'refit': 6}; mse_eff wins {'refit': 1, 'probe': 5}.

- batches=40 · batch_size=100 · n_control=1 · datasets=6 · mean_duty=0.2435897435897436 · duty_breakeven=1.0 · prefer_refit_budget=6/6 · E[refit]/E[probe]≈0.2435897435897436

## Per dataset

| dataset | duty | E[refit]/E[probe] | ref ρ | probe rank_eff | refit rank_eff | refit rank_eff_E | probe mse_eff | refit mse_eff | reading |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 0.1795 | 0.1795 | 0.3021 | 0.0187 | 0.1396 | 0.1361 | 1.06e+05 | 8.192e+05 | duty=0.18: refit wins rank_eff (reject-only vs always-on probe) |
| `beijing_pm25` | 0.1282 | 0.1282 | 0.0844 | 0.4287 | 3.4986 | 3.4112 | -382.4970 | -3488 | duty=0.13: refit wins rank_eff (reject-only vs always-on probe) |
| `stocks_AAPL` | 0.1795 | 0.1795 | 0.2865 | 0.5332 | 3.0073 | 2.9321 | -0.0001382 | -0.0009838 | duty=0.18: refit wins rank_eff (reject-only vs always-on probe) |
| `waymo_proxy` | 0.5641 | 0.5641 | 0.1149 | 0.6035 | 1.0462 | 1.0201 | -0.001571 | -0.001992 | duty=0.56: refit wins rank_eff (reject-only vs always-on probe) |
| `stocks_MSFT` | 0.2308 | 0.2308 | 0.3778 | 0.5215 | 2.1209 | 2.0679 | -0.0001187 | -0.0005367 | duty=0.23: refit wins rank_eff (reject-only vs always-on probe) |
| `stocks_IWM` | 0.1795 | 0.1795 | 0.6084 | 0.1624 | 0.6369 | 0.6210 | -7.459e-05 | -0.0002671 | duty=0.18: refit wins rank_eff (reject-only vs always-on probe) |

## How to read

- **duty**: OnlineRFPerm gate reject rate (Bernoulli planning rate).
- **E[refit]/E[probe] ≈ duty · n_control** — ex-ante budget ratio.
- **rank_eff**: ΔSpearman vs frozen ref / realized FLOPs.
- **rank_eff_E**: same Δρ but / **expected** FLOPs (duty × n_batches).
- **mse_eff**: (uniform − mode) sig-only MSE / FLOPs (positive = cheaper error drop).
- Low duty ⇒ refit is the cheap PO path; probe cost is duty-invariant.

