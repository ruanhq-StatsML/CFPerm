# Stream packs: PO-OOD vs DRE — MSE & regret

Real streams + Waymo kinematics proxy. `batch_size=100`.

| dataset | n | d | shift | √PO MSE | DRE MSE | unif MSE | √/unif | DRE/unif | cumR(√−DRE) | win% vs DRE | best |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 4000 | 19 | 0.3162 | 1.621e+06 | 2.174e+06 | 1.471e+06 | 1.102 | 1.478 | -2.101e+07 | 0.68 | `uniform` |
| `beijing_pm25` | 4000 | 13 | 0.1715 | 2020 | 2929 | 1864 | 1.083 | 1.571 | -3.453e+04 | 0.76 | `uniform` |
| `stocks_SPY` | 4000 | 11 | 0.6367 | 0.0001944 | 0.0001974 | 0.0001809 | 1.075 | 1.092 | -0.0001143 | 0.32 | `inv` |
| `stocks_QQQ` | 4000 | 11 | 0.6395 | 0.0003784 | 0.0004439 | 0.0003502 | 1.080 | 1.268 | -0.002491 | 0.34 | `inv` |
| `affec` | 4000 | 32 | 0.7294 | 4.836 | 12.72 | 4.518 | 1.070 | 2.816 | -299.7 | 0.97 | `uniform` |
| `waymo_proxy` | 4000 | 9 | 0.3672 | 0.01158 | 0.01233 | 0.0106 | 1.093 | 1.163 | -0.02825 | 0.63 | `inv` |

**√PO vs DRE head-to-head:** `6/6`

```python
w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))
```
