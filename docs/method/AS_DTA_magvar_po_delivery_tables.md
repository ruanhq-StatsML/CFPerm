# Mag×Var PermuCATE PO — delivery tables (integrated)

**Verdict:** effects match intent — ready to hand over.
**PO here** $=u=(\varphi-\hat\tau)^2$. Mag×Var + 客群/region on that $u$.

## PO

| Symbol | Definition |
| --- | --- |
| $W$ | REF=0 / LIVE=1 |
| $Y$ | response |
| $\varphi$ | DR pseudo-outcome (PermuCATE) |
| $\hat\tau$ | CATE RF |
| **$u$** | $(\varphi-\hat\tau)^2$ |
| $\mu,s$ | Mag / Var |
| severity | $J + Y$-flip |

## Effect headline

| Suite | Metric | Value |
| --- | --- | --- |
| Tabular 2×2 | acc | **1.00** |
| Token AUC | severity AUC | **1.00** |
| Decision tree | family acc | **1.00** |
| 客群/region | segment hit-rate | **1.00** (4/4) |
| Industrial ROI | cost saved (routed) | **$22,430** |

## 客群 / region Mag×Var

| Segment | Region | sev(ben) | sev(mal) | gap | OK? |
| --- | --- | --- | --- | --- | --- |
| general | hc3_all | 0.04 | 0.41 | 0.37 | yes |
| finance | hc3_finance | 0.09 | 0.45 | 0.36 | yes |
| medicine | hc3_medicine | 0.01 | 0.48 | 0.47 | yes |
| sci | news20_sci | -0.01 | 0.46 | 0.46 | yes |

## Business geometry

| Sensor | Action |
| --- | --- |
| μ↑, s stable | conformal |
| flat, no Y-break | hold |
| mid Y-flip | VAI freeze (segment-local) |
| large Y-flip | broader freeze |

LaTeX: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
