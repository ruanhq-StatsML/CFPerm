# PO-risk IPTW continuous-batch MSE

| mode | MSE_next mean | std | p50 | p90 | MSE_cur mean |
|---|---:|---:|---:|---:|---:|
| `uniform` | 0.1656 | 0.1385 | 0.1166 | 0.2803 | 0.0760 |
| `prop` | 0.1774 | 0.1218 | 0.1363 | 0.3426 | 0.0631 |
| `inv` | 1.8349 | 0.6336 | 1.7955 | 2.5365 | 1.8275 |

**Best next-batch MSE mean:** `uniform`

```python
w = po / po.mean()           # prop
w = (1/po) / (1/po).mean()   # inv
rf.fit(X, y, sample_weight=w)
```
