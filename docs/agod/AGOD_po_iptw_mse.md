# PO-risk IPTW continuous-batch MSE

| mode | MSE_next mean | std | p50 | p90 | MSE_cur mean |
|---|---:|---:|---:|---:|---:|
| `uniform` | 1522.9567 | 1014.7202 | 1489.6241 | 2786.4544 | 574.0128 |
| `prop` | 7154.5548 | 11279.5158 | 2382.6280 | 18141.8068 | 1096.0886 |
| `inv` | 2807.8169 | 2251.6389 | 1696.8640 | 5747.6041 | 941.8731 |

**Best next-batch MSE mean:** `uniform`

```python
w = po / po.mean()           # prop
w = (1/po) / (1/po).mean()   # inv
rf.fit(X, y, sample_weight=w)
```
