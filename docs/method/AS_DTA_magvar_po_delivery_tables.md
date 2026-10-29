# Mag×Var on PermuCATE PO — delivery tables

**PO here** $=u_i=(\varphi_i-\hat\tau_i)^2$ (PermuCATE RF). Not $1-\hat p$.
Mag×Var $(\mu,s)$ on that $u$ is the product triage proxy — cheap and business-aligned.

## 1. PO in one place

| Symbol | Definition |
| --- | --- |
| $W$ | REF=0 / LIVE=1 |
| $Y$ | response (reward or class label) |
| $X$ | covariates (text: TF-IDF→SVD) |
| $\hat\mu_w$ | RF $\mathbb{E}[Y\mid X,W=w]$ |
| $\hat\pi$ | RF / logistic $\mathbb{P}(W=1\mid X)$ |
| $\varphi$ | $(\hat\mu_1-\hat\mu_0)+\frac{(Y-\hat\mu_W)(W-\hat\pi)}{\hat\pi(1-\hat\pi)}$ |
| $\hat\tau$ | RF regress $\varphi\sim X$ |
| **$u_i$** | **unit PO** $(\varphi_i-\hat\tau_i)^2$ |
| $\mu,s$ | Mag$=\mathbb{E}[u]$; Var$=\mathrm{Std}(r_b)$ |
| $J$ | Mag×Var merge |

## 2. Business geometry

| Sensor | Typical shift | Cheap action |
| --- | --- | --- |
| $\mu\uparrow$, $s$ stable | Benign / label delay | conformal / wait |
| flat, no $Y$-break | Quiet | hold |
| mid $Y$-flip + Mag×Var | Malignant sparse | VAI freeze |
| large $Y$-flip | Malignant dense | broader freeze / retrain |

## 3. Effect — tabular 2×2

|  | pred benign | pred malignant |
| --- | --- | --- |
| true benign | 1 | 0 |
| true malignant | 0 | 1 |

Acc **1.00**. Malignant $J\!\approx\!114$ vs benign $J\!\approx\!-0.21$.

## 4. Effect — text token AUC (PermuCATE $u$)

| Board | AUC(severity) | AUC($J$) | AUC($Y$-flip) | AUC($\Delta\mu$) |
| --- | --- | --- | --- | --- |
| hc3_all | **1.00** | 0.61 | 1.00 | 0.69 |
| hc3_finance | **1.00** | 0.42 | 1.00 | 0.56 |
| news20_sci | **1.00** | 0.25 | 1.00 | 0.22 |
| **pooled** | **1.00** | 0.43 | 1.00 | 0.48 |

Severity $=J+Y$-flip. Mag alone is soft; Mag×Var + response integrity is the product score.

## 5. Effect — decision tree

| Board | family acc | branch acc |
| --- | --- | --- |
| hc3_all | 1.00 | 0.75 |
| hc3_finance | 1.00 | 0.75 |
| news20_sci | 1.00 | 0.75 |
| **mean** | **1.00** | 0.75 |

Family confusion (quick 4-way): quiet/benign safe block 6/6; sparse 3/3; dense 3/3.

## 6. Honesty

| Path | Unit risk | PermuCATE RF? |
| --- | --- | --- |
| Tabular Mag×Var | $(\varphi-\hat\tau)^2$ | **Yes** |
| Text token AUC | $(\varphi-\hat\tau)^2$ | **Yes** |
| Decision tree | $(\varphi-\hat\tau)^2$ | **Yes** |
| Token VAI | frozen $(\varphi-\hat\tau)^2$ | **Yes** |

LaTeX twin: `docs/method/AS_DTA_magvar_po_product_delivery.tex`
