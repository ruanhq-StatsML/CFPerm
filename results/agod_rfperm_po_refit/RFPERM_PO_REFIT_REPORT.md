# OnlineRFPerm → T=0/T=1 PO re-fit → √PO (post-hoc)

Default = **uniform**. On significant OnlineRFPerm reject only:

- `T=0`: recent control batch(es) before the pair
- `T=1`: previous ∪ current batch
- re-fit μ0 (optional μ1); `PO=|Y−μ0(X)|`; `w=√PO` on current batch

| dataset | unif | always-√ | gated | **refit** | dre | refit-duty | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 1.471e+06 | 1.621e+06 | 1.459e+06 | **1.461e+06** | 2.174e+06 | 0.18 | `sqrt_gated` |
| `beijing_pm25` | 1864 | 2020 | 1899 | **1868** | 2929 | 0.13 | `uniform` |
| `stocks_AAPL` | 0.0007237 | 0.0008474 | 0.0007397 | **0.0007415** | 0.0008296 | 0.18 | `uniform` |
| `waymo_proxy` | 0.0106 | 0.01158 | 0.01129 | **0.01108** | 0.01233 | 0.56 | `uniform` |
| `stocks_MSFT` | 0.0004298 | 0.0004715 | 0.0004456 | **0.0004447** | 0.0004434 | 0.23 | `uniform` |
| `stocks_IWM` | 0.0002679 | 0.000291 | 0.0002783 | **0.0002763** | 0.0002774 | 0.18 | `uniform` |

**refit vs DRE:** `5/6`
**refit vs always-√:** `6/6`
**refit vs probe-gated:** `4/6`

Expectation: uniform slightly best overall; refit fires sparsely on
clearly-shifted batches and should beat always-√ / DRE there.

## Significant batches only

| dataset | n_sig | duty | uniform | sqrt | sqrt_gated | sqrt_gated_refit | dre | best |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 6/38 | 0.16 | 8.602e+05 | 8.455e+05 | 7.839e+05 | **7.977e+05** | 1.642e+06 | `sqrt_gated` |
| `beijing_pm25` | 5/38 | 0.13 | 2364 | 2528 | 2627 | **2391** | 3577 | `uniform` |
| `stocks_AAPL` | 6/38 | 0.16 | 0.0007784 | 0.0008755 | 0.0008796 | **0.0008912** | 0.0008401 | `uniform` |
| `waymo_proxy` | 22/38 | 0.58 | 0.009992 | 0.01101 | 0.01118 | **0.01082** | 0.01154 | `uniform` |
| `stocks_MSFT` | 8/38 | 0.21 | 0.0006164 | 0.0006963 | 0.0006919 | **0.0006875** | 0.0005892 | `dre` |
| `stocks_IWM` | 7/38 | 0.18 | 0.0005025 | 0.0005516 | 0.0005588 | **0.0005479** | 0.0005117 | `uniform` |

**Wins (sig-only):** `uniform`=4, `sqrt`=0, `sqrt_gated`=1, `sqrt_gated_refit`=0, `dre`=1
