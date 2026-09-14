# PO^{1/3} vs PO^{1/2} adaptation (OnlineRFPerm-gated)

Default = **uniform**. On significant reject, T0/T1 re-fit then:
`gated_sqrt`: w∝PO^{1/2}; `gated_cbrt`: w∝PO^{1/3} (softer).

| dataset | unif | always-√ | always-∛ | gated-√ | **gated-∛** | dre | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 1.471e+06 | 1.621e+06 | 1.567e+06 | 1.461e+06 | **1.462e+06** | 2.174e+06 | `gated_sqrt` |
| `beijing_pm25` | 1864 | 2020 | 1995 | 1868 | **1885** | 2929 | `uniform` |
| `stocks_AAPL` | 0.0007237 | 0.0008474 | 0.0007879 | 0.0007415 | **0.0007364** | 0.0008296 | `uniform` |
| `waymo_proxy` | 0.0106 | 0.01158 | 0.01136 | 0.01108 | **0.01089** | 0.01233 | `uniform` |
| `stocks_MSFT` | 0.0004298 | 0.0004715 | 0.0004639 | 0.0004447 | **0.0004331** | 0.0004434 | `uniform` |
| `stocks_IWM` | 0.0002679 | 0.000291 | 0.0002831 | 0.0002763 | **0.0002715** | 0.0002774 | `uniform` |

**Wins:** `uniform`=5, `sqrt`=0, `cbrt`=0, `gated_sqrt`=1, `gated_cbrt`=0, `dre`=0
**gated-∛ ≤ gated-√:** `4/6`

```python
w = po ** (1/3)   # softer than sqrt; closer to uniform
```
