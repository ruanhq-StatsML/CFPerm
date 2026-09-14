# OnlineRFPerm-gated √PO — multi-pack benchmark

Gate = PDF OnlineRFPerm (rank/EWMA p + alpha-investing FDR).
Post-hoc: **only after significant p** → `w=√PO`; else uniform.

| dataset | unif | always-√ | **gated-√** | dre | duty | cumR(g−dre) | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `stocks_MSFT` | 0.0004298 | 0.0004715 | **0.0004456** | 0.0004434 | 0.23 | 8.603e-05 | `uniform` |
| `stocks_IWM` | 0.0002679 | 0.000291 | **0.0002783** | 0.0002774 | 0.18 | 3.486e-05 | `uniform` |
| `stocks_AAPL` | 0.0007237 | 0.0008474 | **0.0007397** | 0.0008296 | 0.18 | -0.003417 | `uniform` |
| `metro_interstate` | 1.471e+06 | 1.621e+06 | **1.459e+06** | 2.174e+06 | 0.18 | -2.718e+07 | `sqrt_gated` |
| `beijing_pm25` | 1864 | 2020 | **1899** | 2929 | 0.13 | -3.913e+04 | `uniform` |
| `waymo_proxy` | 0.0106 | 0.01158 | **0.01129** | 0.01233 | 0.56 | -0.0396 | `uniform` |

**gated-√ vs DRE:** `4/6`
**gated-√ vs always-√:** `6/6`

```python
T = MSE(f_ref, batch) - E_ref
p = rank/EWMA vs historical T; online FDR → reject  # small p = significance
w = sqrt(PO) if reject else 1                       # post-hoc weighting
```
