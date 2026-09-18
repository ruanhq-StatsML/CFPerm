# PO-risk as OOD score — streaming real-data (bs=100)

No synth. Continuous → next **MSE** (↓); discrete → next **Acc** (↑).
PO-√ soft IPTW vs logistic density-ratio (DRE).

| dataset | task | uniform | prop | sqrt | inv | dre | best | √PO vs DRE |
|---|---|---:|---:|---:|---:|---:|---|---|
| `affec` | mse | **3.8743** | 4.8779 | *4.1471* | 7.0848 | 11.7758 | `uniform` | −7.6286 MSE |
| `tencent` | mse | **0.0007** | 0.0014 | *0.0014* | 0.0024 | 0.0010 | `uniform` | +0.0004 MSE |
| `msrvtt` | acc | 0.4872 | 0.4633 | *0.4944* | 0.4917 | **0.4967** | `dre` | -0.0022 Acc |
| `coco_time` | acc | 0.6361 | 0.6306 | *0.6311* | **0.6422** | 0.6283 | `inv` | +0.0028 Acc |
| `fashion_iq` | acc | **0.6033** | 0.5989 | *0.6017* | 0.5939 | 0.5811 | `uniform` | +0.0206 Acc |
| `indiana_cxr` | acc | **0.5256** | 0.5056 | *0.5244* | 0.5217 | 0.4439 | `uniform` | +0.0806 Acc |

**Wins:** `uniform`=4, `prop`=0, `sqrt`=0, `inv`=1, `dre`=1
**sqrt vs dre (head-to-head):** `4/6` favor PO-√ OOD over DRE.

```python
w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))  # OOD score → soft IPTW
model.fit(X, y, sample_weight=w / w.mean())
```

DRE baseline: logistic `w ∝ p(cur|x)/p(ref|x)` on X only — ignores label risk.
