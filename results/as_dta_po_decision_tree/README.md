# Product use case — Mag×Var PermuCATE decision tree

## The tree

```
Bootstrap Var(u), u=(φ−τ̂)² — variance stable?
├─ YES → mean↑?
│   ├─ YES → Benign covariate / Label delay (time traj)
│   └─ NO  → Quiet
└─ NO  → kurtosis > 10?
    ├─ YES → Malignant sparse (VAI token blocks)
    └─ NO  → Malignant dense / Schema (traj + VAI)
```

## Estimator

- **PermuCATE RF** text: $Y$=response, $W$=REF/LIVE, $X$=TF-IDF→SVD
- Unit $u_i=(\varphi_i-\hat\tau_i)^2$ — **not** $1-\hat p(Y\mid x)$
- VAI = token drop → re-featurize → frozen PermuCATE models

## Product use case

**Shift triage gate for integrity / moderation / ranking monitors**

One bootstrap Mag×Var panel on PermuCATE unit risk replaces 'OOD high ⇒ retrain'. Stable-var mean lifts get cheap conformal; kurtotic sparse attacks get VAI-guided freeze; dense/schema get heavier review.

## Effect

- mean family accuracy: **1.00**
- mean branch accuracy: **0.75**

| board | family_acc | branch_acc |
| --- | --- | --- |
| hc3_all | 1.00 | 0.75 |
| hc3_finance | 1.00 | 0.75 |
| news20_sci | 1.00 | 0.75 |

![](po_decision_tree.png)

```bash
python3 scripts/run_as_dta_po_decision_tree.py --quick
python3 scripts/run_as_dta_po_decision_tree.py
```

Code: `Python/src/as_dta_po_decision_tree.py`
