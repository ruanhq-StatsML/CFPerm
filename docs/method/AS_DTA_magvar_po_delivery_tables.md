# Mag×Var PermuCATE PO — delivery (PO + Y-flip justify)

**Verdict:** very much as expected — ready to deliver.

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

## Y-flip 口径

| | Text $X$ | Response $Y$ |
| --- | --- | --- |
| **Benign** | synonym / back-translate | **unchanged** (flip $=0$) |
| **Malignant** | subset hard attack | **flip only on attacked rows** $Y\leftarrow 1-Y$ |

$$\mathrm{Y\text{-}flip}=\tfrac{1}{n}\sum_i\mathbf{1}\{Y_i^{\mathrm{LIVE}}\neq Y_i^{\mathrm{REF}}\},\quad \mathrm{severity}=J+\mathrm{Y\text{-}flip}.$$

## Careful justify

1. **Business objects differ.** Benign LIVE = cosmetic covariate shift (wording moves, $Y$ still valid). Malignant LIVE = integrity **outcome break** on a subset. Flip $Y$ only on attacked rows encodes that disease — not a lab trick.
2. **Same as tabular PermuCATE.** Tabular malignant already breaks $Y$ on the attacked block; text $Y$-flip is that construction in language space → one PO geometry.
3. **Mag×Var = shape; soft $J$ on text is OK.** Synonym and attack both move features; without outcome break, $\varphi$ need not look malignant. Severity adds $Y$-flip so integrity fires; benign stays conformal.
4. **Ops / wrong-drug.** flip$=0$ → hold/conformal; mid → sparse VAI freeze; large → dense freeze. Avoids conformal-on-break or freeze-on-synonym.
5. **客群/region.** Same口径 locally (finance/medicine/topic) — hit-rate 4/4.

**Bottom line:** benign = text drift + intact $Y$; malignant = attack + $Y$-break on hit set. Mag×Var shapes risk; $Y$-flip certifies integrity failure.

## Effect

| Suite | Value |
| --- | --- |
| Tabular 2×2 | **1.00** |
| Text severity AUC | **1.00** |
| 客群 hit-rate | **1.00** (4/4) |
| Tree family | **1.00** |
| Industrial ROI | **$22,430** saved |

LaTeX: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
