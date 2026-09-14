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

## Significant batches only

| dataset | n_sig | duty | uniform | sqrt | cbrt | gated_sqrt | gated_cbrt | dre | best |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 6/38 | 0.16 | 8.602e+05 | 8.455e+05 | 8.537e+05 | 7.977e+05 | **8.034e+05** | 1.642e+06 | `gated_sqrt` |
| `beijing_pm25` | 5/38 | 0.13 | 2364 | 2528 | 2682 | 2391 | **2523** | 3577 | `uniform` |
| `stocks_AAPL` | 6/38 | 0.16 | 0.0007784 | 0.0008755 | 0.0008205 | 0.0008912 | **0.0008584** | 0.0008401 | `uniform` |
| `waymo_proxy` | 22/38 | 0.58 | 0.009992 | 0.01101 | 0.01082 | 0.01082 | **0.01049** | 0.01154 | `uniform` |
| `stocks_MSFT` | 8/38 | 0.21 | 0.0006164 | 0.0006963 | 0.0006804 | 0.0006875 | **0.0006323** | 0.0005892 | `dre` |
| `stocks_IWM` | 7/38 | 0.18 | 0.0005025 | 0.0005516 | 0.0005267 | 0.0005479 | **0.000522** | 0.0005117 | `uniform` |

**Wins (sig-only):** `uniform`=4, `sqrt`=0, `cbrt`=0, `gated_sqrt`=1, `gated_cbrt`=0, `dre`=1
