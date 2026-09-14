# OnlineRFPerm-gated √PO — held-out stream packs

Gate = PDF OnlineRFPerm (rank/EWMA p + alpha-investing FDR).
Reweight only when reject: `w=√PO`; else uniform.

| dataset | unif | always-√ | **gated-√** | dre | duty | cumR(g−dre) | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `stocks_AAPL` | 0.0007237 | 0.0008474 | **0.0007397** | 0.0008296 | 0.18 | -0.003417 | `uniform` |
| `metro_interstate` | 1.471e+06 | 1.621e+06 | **1.459e+06** | 2.174e+06 | 0.18 | -2.718e+07 | `sqrt_gated` |
| `beijing_pm25` | 1864 | 2020 | **1899** | 2929 | 0.13 | -3.913e+04 | `uniform` |
| `waymo_proxy` | 0.0106 | 0.01158 | **0.01129** | 0.01233 | 0.56 | -0.0396 | `uniform` |

```python
# OnlineRFPerm gate (PDF Alg.1)
T = MSE(f_ref, batch) - E_ref
p = rank/EWMA vs historical T; online FDR → reject
w = sqrt(PO) if reject else 1
```
