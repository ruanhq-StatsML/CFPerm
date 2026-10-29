# Mag×Var PermuCATE PO — MVP + lit + fused results

**Objective:** separate benign vs malignant. **Verdict:** meets expectation — ship.

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

### Headline — did we separate?

Headline = product KPI check (not a score dump). Four suites, one question, four answers that must agree before close.

| Suite | KPI | What it proves |
| --- | --- | --- |
| Tabular Mag×Var | 2×2 acc **1.00** | shape alone: $J$ quiet vs explode (−0.21 vs +114) |
| Text severity | AUC **1.00** | $J$ soft; $Y$-flip carries integrity |
| 客群/region | hit-rate **1.00** | same口径 local — not one global alarm |
| Decision tree | family **1.00** | typed actuator path (branch 0.75 ok) |

**Headline logic:** any suite fails → do not ship that path; all four hit → MVP triage ready to close.

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

## Close — sync the contract, then ship

Close mirrors headline (not a second cheer).

| Step | Locked | Status |
| --- | --- | --- |
| Objective | separate + typed medicine | §2 |
| PO / How | $u$; Mag×Var + $Y$-flip; 客群-local | §3–4 |
| Justify | cosmetic $X$ ≠ outcome break | §5 |
| Headline | four suites separate on same $u$ | §6 |
| **Ship** | freeze口径; sync PR; stop | **here** |

**Close logic (=「同步 PR 并收尾」):** (1) headline green → (2) PO is here → (3) freeze口径 → (4) sync delivery note to PR and stop. No Peak-first rewrite.

LaTeX: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
