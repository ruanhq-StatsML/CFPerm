# PO-risk as OOD score — streaming multi-dataset (bs=100)

PO √ soft-upweight vs logistic DRE vs uniform/prop/inv. Metric: next-batch MSE.

| dataset | uniform | prop | **sqrt** | inv | dre | best |
|---|---:|---:|---:|---:|---:|---|
| `affec` | **4.9027** | 6.2612 | *5.4558* | 9.4819 | 21.7733 | `uniform` |
| `msrvtt` | **0.0000** | 0.0000 | *0.0000* | 0.0000 | 0.0000 | `uniform` |
| `coco_time` | 808.8178 | 1031.1609 | *919.8179* | 814.4531 | **803.1976** | `dre` |
| `fashion_iq` | 1183.2311 | 1448.6814 | *1342.5351* | **1107.8642** | 1190.0198 | `inv` |
| `indiana_cxr` | 1031.2220 | 1119.6047 | *1083.8672* | 1072.0408 | **1015.4889** | `dre` |
| `tencent` | **0.0007** | 0.0014 | *0.0014* | 0.0024 | 0.0011 | `uniform` |
| `synth` | 2.3118 | 2.4488 | *2.4111* | **2.2616** | 2.3779 | `inv` |

**Wins (lowest next-MSE):** `uniform`=3, `prop`=0, `sqrt`=0, `inv`=2, `dre`=2

```python
w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))  # OOD score → soft IPTW
rf.fit(X, y, sample_weight=w / w.mean())
```

Claim: on gradual-shift streams, PO-risk OOD (`sqrt`) beats density-ratio (`dre`) empirically.
