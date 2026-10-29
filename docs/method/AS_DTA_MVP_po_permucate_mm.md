# Professor MVP — Mag×Var on PermuCATE RF + multimodal attribution

## Pitch (one screen)

> When mean PO-risk rises, **Magnitude × Variance** of unit PO across bootstrap
> sub-environments resolves **benign vs malignant shift** (2×2).  
> This causal MVP estimates unit risk with **your PermuCATE stack**:  
> **RF outcome models** $\hat\mu_0,\hat\mu_1$, RF propensity $\hat\pi$, RF $\hat\tau$ on DR $\varphi$,  
> $u_i=(\varphi_i-\hat\tau_i)^2$.  
> Multimodal attribution = group-permute / τ̂-VIMP mass on V/A/T — same “where is the PO mass” intuition.

## Estimator honesty

| path | PO unit | RF outcome? |
| --- | --- | --- |
| AIGC text Mag×Var boards (earlier) | proxy $1-\hat p(Y\mid x)$ | No (TF-IDF + logistic) |
| **This MVP** `as_dta_permucate_magvar.py` | $(\varphi-\hat\tau)^2$ | **Yes** — `RandomForestRegressor` μ₀/μ₁ + RF π + RF τ̂ |

Formulas align with `Python/src/permuCATE_vimp.py` (`po_risk` / `permuCATE_fit_nuisance`, `model_m='rf_regressor'`).

$$
\varphi=\hat\mu_1-\hat\mu_0+\frac{(Y-\hat\mu_W)(W-\hat\pi)}{\hat\pi(1-\hat\pi)},\quad
u=(\varphi-\hat\tau)^2,\quad
\mu=\mathbb{E}[u],\quad
s=\mathrm{Std}(r_b).
$$

## 2×2 (this run)

|  | pred benign | pred malignant |
| --- | --- | --- |
| true benign | 1 | 0 |
| true malignant | 0 | 1 |

acc = 1.0 · malignant $J$ strongly amplified by variance merge.

## Multimodal attribution — justify

1. Mag×Var answers **shift type** (uniform vs concentrated PO mass over sub-envs).  
2. Modality shares answer **channel** (which of V/A/T carries that mass):
   - **τ̂ VIMP mass** by feature group (stable primary)
   - **group-permute** Δ$R_{\mathrm{PO}}$ (destroy one modality among LIVE rows)

Same story as Mag×Var: risk concentrates on a structured subset — here a **modality block** instead of a row subset. Actuator: freeze the high-share modality; conformal the rest.

## Run

```bash
python3 scripts/run_as_dta_permucate_magvar.py
```

Artifacts: `results/as_dta_permucate_magvar/`
