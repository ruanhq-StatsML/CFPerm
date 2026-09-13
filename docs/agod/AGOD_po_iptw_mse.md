# PO-risk IPTW continuous-batch MSE

Soft high-PO upweight (no DRE/DGA):

```python
w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))
w = w / w.mean()   # mean-1 normalize for RF/LR scale
rf.fit(X, y, sample_weight=w)
```

## Affec RF (10×256) — `sqrt` wins

| mode | MSE_next mean | std | p50 | p90 | MSE_cur mean |
|---|---:|---:|---:|---:|---:|
| `uniform` | 0.1656 | 0.1385 | 0.1166 | 0.2803 | 0.0760 |
| `prop` | 0.1774 | 0.1218 | 0.1363 | 0.3426 | 0.0631 |
| **`sqrt`** | **0.1370** | **0.1039** | 0.1185 | **0.2373** | 0.0363 |
| `inv` | 1.8349 | 0.6336 | 1.7955 | 2.5365 | 1.8275 |

- `prop` (raw ∝PO): overfits current batch → next-MSE worse than uniform.
- `sqrt` (∝√PO): soft upweight → best next-MSE + lowest std/p90.
- `inv` (∝1/PO): collapses on Affec.

## Synth (hard rotate) — `inv` slightly ahead

| mode | MSE_next mean | std | p50 | p90 | MSE_cur mean |
|---|---:|---:|---:|---:|---:|
| `uniform` | 3.5242 | 1.2434 | 3.4737 | 5.0107 | 0.7561 |
| `prop` | 4.0115 | 1.7046 | 3.5730 | 6.1364 | 0.7647 |
| `sqrt` | 3.7390 | 1.4974 | 3.4982 | 5.6240 | 0.6034 |
| `inv` | 3.4237 | 1.1270 | 3.4679 | 4.8196 | 2.2811 |

Default mode in `agod.po_iptw.po_iptw_weights` is **`sqrt`**.
