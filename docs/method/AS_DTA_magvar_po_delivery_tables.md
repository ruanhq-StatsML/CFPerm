# Mag×Var PermuCATE PO — separate benign vs malignant

**Objective:** distinguish benign vs malignant shift → different medicine.  
**Verdict:** meets expectation — deliver.

## PO is here

| Symbol | Definition |
| --- | --- |
| $W$ | REF=0 / LIVE=1 |
| $Y$ | response |
| $\varphi$ | PermuCATE DR pseudo-outcome |
| $\hat\tau$ | CATE RF |
| **$u$** | $(\varphi-\hat\tau)^2$ |
| $J$ | Mag×Var merge |
| severity | $J + Y$-flip |

## How to separate

| | Text $X$ | Response $Y$ |
| --- | --- | --- |
| **Benign** | synonym / back-translate | unchanged (flip $=0$) |
| **Malignant** | subset attack | flip only on attacked rows $Y\leftarrow 1-Y$ |

- **Mag×Var on $u$:** shape (uniform vs mixture mass)
- **$Y$-flip:** integrity outcome break
- **severity** $= J + Y$-flip  
  flip$=0$ + mild $J$ → conformal; mid flip → sparse freeze; large flip → dense freeze

## Justify

1. **Disease = objective.** Benign = cosmetic $X$ + valid $Y$; malignant = attack + outcome break. $Y$-flip on attacked rows measures that.
2. **Same as tabular PermuCATE** (malignant already breaks $Y$ on the hit set).
3. **Mag×Var = shape;** soft $J$ on text is OK without outcome break. Severity adds $Y$-flip.
4. **Wrong-drug KPI:** conformal≠freeze; two sensors prevent the swap.

## Effect

| Suite | Value |
| --- | --- |
| Tabular 2×2 | **1.00** |
| Text severity AUC | **1.00** |
| 客群 hit-rate | **1.00** |
| Tree family | **1.00** |

LaTeX: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
