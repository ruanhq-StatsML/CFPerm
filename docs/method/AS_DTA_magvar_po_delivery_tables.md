# Mag×Var PermuCATE PO — MVP + lit + fused results

**Objective:** separate benign vs malignant. **Verdict:** meets expectation.

## Literature (short)

| Thread | Takeaway |
| --- | --- |
| PO / DR / CATE | Monitor $u=(\varphi-\hat\tau)^2$ |
| Uplift / batch PO-risk | LIVE needs PO-risk, not holdout alone |
| OOD / mean risk | Mean↑ mixes two diseases |
| **This MVP** | Mag×Var shape + $Y$-flip integrity |

## PO is here

$u=(\varphi-\hat\tau)^2$; severity $=J+Y$-flip. Benign keeps $Y$; malignant flips $Y$ on attacked rows.

## Results (fused)

### Headline

| Suite | Metric |
| --- | --- |
| Tabular Mag×Var 2×2 | acc **1.00** |
| Text severity AUC | pooled **1.00** |
| 客群/region | hit-rate **1.00** (4/4) |
| Decision-tree family | **1.00** (branch 0.75) |

### Tabular 2×2 + $J$

|  | pred benign | pred malignant |
| --- | --- | --- |
| true benign | 1 | 0 |
| true malignant | 0 | 1 |

| Shift | $\Delta\mu$ | $\Delta s$ | $J$ |
| --- | --- | --- | --- |
| benign | −0.21 | −0.25 | −0.21 |
| malignant | +9.70 | +3.02 | +114 |

### Text severity AUC

| Board | AUC(severity) | AUC($J$) | AUC($Y$-flip) | AUC($\Delta\mu$) |
| --- | --- | --- | --- | --- |
| hc3_all | **1.00** | 0.61 | 1.00 | 0.69 |
| hc3_finance | **1.00** | 0.42 | 1.00 | 0.56 |
| news20_sci | **1.00** | 0.25 | 1.00 | 0.22 |
| **pooled** | **1.00** | 0.43 | 1.00 | 0.48 |

| Board | sev(ben) | sev(mal) | $J$(ben) | $J$(mal) |
| --- | --- | --- | --- | --- |
| hc3_all | 0.04 | 0.51 | 0.01 | 0.06 |
| hc3_finance | 0.08 | 0.49 | 0.05 | 0.06 |
| news20_sci | 0.08 | 0.42 | 0.04 | −0.03 |

### 客群 / region

| Segment | Region | sev(ben) | sev(mal) | gap | OK? |
| --- | --- | --- | --- | --- | --- |
| general | hc3_all | 0.04 | 0.41 | 0.37 | yes |
| finance | hc3_finance | 0.09 | 0.45 | 0.36 | yes |
| medicine | hc3_medicine | 0.01 | 0.48 | 0.47 | yes |
| sci | news20_sci | −0.01 | 0.46 | 0.46 | yes |

### Decision tree

| Board | family | branch |
| --- | --- | --- |
| hc3_all / finance / news20 | 1.00 | 0.75 |
| **mean** | **1.00** | 0.75 |

|  | quiet/benign | sparse | dense |
| --- | --- | --- | --- |
| true quiet/benign | 6 | 0 | 0 |
| true sparse | 0 | 3 | 0 |
| true dense | 0 | 0 | 3 |

LaTeX: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
