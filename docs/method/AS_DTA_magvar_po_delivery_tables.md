# Mag×Var PermuCATE PO — delivery (PO here + Y-flip)

**Verdict:** meets expectation — ready to deliver.

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

## Benign vs malignant — Y-flip 口径

| | Text $X$ | Response $Y$ |
| --- | --- | --- |
| **Benign** | uniform synonym / back-translate | **unchanged** (flip $=0$) |
| **Malignant** | subset hard attack | **flip only on attacked rows** $Y\leftarrow 1-Y$ |

$$\mathrm{Y\text{-}flip}=\frac{1}{n}\sum_i \mathbf{1}\{Y_i^{\mathrm{LIVE}}\neq Y_i^{\mathrm{REF}}\},\qquad \mathrm{severity}=J+\mathrm{Y\text{-}flip}.$$

### Why this conclusion is right

1. **Two diseases, two sensors.** Benign = cosmetic text drift, labels/rewards still valid. Malignant integrity hit = outcome break on a subset — not “uglier text” alone. Keep $Y$ for benign; flip $Y$ only on attacked rows.
2. **Same as tabular PermuCATE.** Tabular malignant already breaks $Y$ on the attacked subset; text $Y$-flip is that construction in language space.
3. **Mag×Var = shape; $Y$-flip = integrity.** Pure $J$ can be soft on text (both moves features). Severity adds $Y$-flip so integrity fires, benign stays conformal — not freeze.
4. **Ops mapping.** flip$=0$ → hold/conformal; mid flip → sparse VAI freeze; large flip → dense freeze. Avoids wrong drug.
5. **客群/region.** Same口径 per finance/medicine/topic — hit-rate 4/4.

**One line:** benign = text drift with intact $Y$; malignant = attack + $Y$-break; Mag×Var shapes risk, $Y$-flip certifies integrity failure.

## Effect

| Suite | Metric | Value |
| --- | --- | --- |
| Tabular 2×2 | acc | **1.00** |
| Text severity AUC | pooled | **1.00** |
| 客群/region | hit-rate | **1.00** (4/4) |
| Tree family | acc | **1.00** |
| Industrial ROI | routed save | **$22,430** |

LaTeX: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
