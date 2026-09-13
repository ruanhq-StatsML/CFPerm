# PO-risk IPTW continuous-batch MSE

| mode | MSE_next mean | std | p50 | p90 | MSE_cur mean |
|---|---:|---:|---:|---:|---:|
| `uniform` | 3.5242 | 1.2434 | 3.4737 | 5.0107 | 0.7561 |
| `prop` | 4.0115 | 1.7046 | 3.5730 | 6.1364 | 0.7647 |
| `inv` | 3.4237 | 1.1270 | 3.4679 | 4.8196 | 2.2811 |

**Best next-batch MSE mean:** `inv`

```python
w = po / po.mean()           # prop
w = (1/po) / (1/po).mean()   # inv
rf.fit(X, y, sample_weight=w)
```
