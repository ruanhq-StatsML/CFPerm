# OnlineRFPerm-gated √PO — held-out stream packs

Gate = PDF OnlineRFPerm (rank/EWMA p + alpha-investing FDR).
Reweight only when reject: `w=√PO`; else uniform.

| dataset | unif | always-√ | **gated-√** | dre | duty | cumR(g−dre) | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `stocks_MSFT` | 0.0004298 | 0.0004715 | **0.0004456** | 0.0004434 | 0.23 | 8.603e-05 | `uniform` |
| `stocks_IWM` | 0.0002679 | 0.000291 | **0.0002783** | 0.0002774 | 0.18 | 3.486e-05 | `uniform` |

```python
# OnlineRFPerm gate (PDF Alg.1)
T = MSE(f_ref, batch) - E_ref
p = rank/EWMA vs historical T; online FDR → reject
w = sqrt(PO) if reject else 1
```
