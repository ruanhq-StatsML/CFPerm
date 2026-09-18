# Significant-batch-only next-MSE

Non-reject steps: gated modes = **uniform** → exclude from the mean.
Compare methods only where OnlineRFPerm opens the gate.

## RFPerm → T0/T1 √PO re-fit

| dataset | n_sig | duty | uniform | sqrt | sqrt_gated | sqrt_gated_refit | dre | best |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 6/38 | 0.16 | 8.602e+05 | 8.455e+05 | 7.839e+05 | **7.977e+05** | 1.642e+06 | `sqrt_gated` |
| `beijing_pm25` | 5/38 | 0.13 | 2364 | 2528 | 2627 | **2391** | 3577 | `uniform` |
| `stocks_AAPL` | 6/38 | 0.16 | 0.0007784 | 0.0008755 | 0.0008796 | **0.0008912** | 0.0008401 | `uniform` |
| `waymo_proxy` | 22/38 | 0.58 | 0.009992 | 0.01101 | 0.01118 | **0.01082** | 0.01154 | `uniform` |
| `stocks_MSFT` | 8/38 | 0.21 | 0.0006164 | 0.0006963 | 0.0006919 | **0.0006875** | 0.0005892 | `dre` |
| `stocks_IWM` | 7/38 | 0.18 | 0.0005025 | 0.0005516 | 0.0005588 | **0.0005479** | 0.0005117 | `uniform` |

**Wins (sig-only):** `uniform`=4, `sqrt`=0, `sqrt_gated`=1, `sqrt_gated_refit`=0, `dre`=1

Relative to uniform (sig-only, ↓ better):

| dataset | sqrt | sqrt_gated | sqrt_gated_refit | dre |
|---|---:|---:|---:|---:|
| `metro_interstate` | 0.983× | 0.911× | 0.927× | 1.909× |
| `beijing_pm25` | 1.069× | 1.111× | 1.012× | 1.513× |
| `stocks_AAPL` | 1.125× | 1.130× | 1.145× | 1.079× |
| `waymo_proxy` | 1.102× | 1.119× | 1.083× | 1.155× |
| `stocks_MSFT` | 1.130× | 1.122× | 1.115× | 0.956× |
| `stocks_IWM` | 1.098× | 1.112× | 1.090× | 1.018× |

## PO^{1/3} vs PO^{1/2} (same gate)

| dataset | n_sig | duty | uniform | sqrt | cbrt | gated_sqrt | gated_cbrt | dre | best |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 6/38 | 0.16 | 8.602e+05 | 8.455e+05 | 8.537e+05 | 7.977e+05 | **8.034e+05** | 1.642e+06 | `gated_sqrt` |
| `beijing_pm25` | 5/38 | 0.13 | 2364 | 2528 | 2682 | 2391 | **2523** | 3577 | `uniform` |
| `stocks_AAPL` | 6/38 | 0.16 | 0.0007784 | 0.0008755 | 0.0008205 | 0.0008912 | **0.0008584** | 0.0008401 | `uniform` |
| `waymo_proxy` | 22/38 | 0.58 | 0.009992 | 0.01101 | 0.01082 | 0.01082 | **0.01049** | 0.01154 | `uniform` |
| `stocks_MSFT` | 8/38 | 0.21 | 0.0006164 | 0.0006963 | 0.0006804 | 0.0006875 | **0.0006323** | 0.0005892 | `dre` |
| `stocks_IWM` | 7/38 | 0.18 | 0.0005025 | 0.0005516 | 0.0005267 | 0.0005479 | **0.000522** | 0.0005117 | `uniform` |

**Wins (sig-only):** `uniform`=4, `sqrt`=0, `cbrt`=0, `gated_sqrt`=1, `gated_cbrt`=0, `dre`=1

Relative to uniform (sig-only, ↓ better):

| dataset | sqrt | cbrt | gated_sqrt | gated_cbrt | dre |
|---|---:|---:|---:|---:|---:|
| `metro_interstate` | 0.983× | 0.992× | 0.927× | 0.934× | 1.909× |
| `beijing_pm25` | 1.069× | 1.135× | 1.012× | 1.068× | 1.513× |
| `stocks_AAPL` | 1.125× | 1.054× | 1.145× | 1.103× | 1.079× |
| `waymo_proxy` | 1.102× | 1.083× | 1.083× | 1.050× | 1.155× |
| `stocks_MSFT` | 1.130× | 1.104× | 1.115× | 1.026× | 0.956× |
| `stocks_IWM` | 1.098× | 1.048× | 1.090× | 1.039× | 1.018× |

```
metric = mean(mse_next[t] for t where reject_t)
# non-significant t: all gated modes ≡ uniform → drop
```
