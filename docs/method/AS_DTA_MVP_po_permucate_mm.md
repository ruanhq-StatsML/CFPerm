# PermuCATE RF Mag×Var + multimodal attribution MVP

## Estimator (important)

- **Yes — RandomForest outcome models** (μ̂₀, μ̂₁), RF propensity, RF τ̂ on DR φ.
- Unit PO-risk $u_i=(\varphi_i-\hat\tau_i)^2$ (same family as `permuCATE_vimp.po_risk`).
- Earlier AIGC text boards used supervised proxy $1-\hat p(Y\mid x)$ — different path.

## 2×2

|  | pred benign | pred malignant |
| --- | --- | --- |
| true benign | 1 | 0 |
| true malignant | 0 | 1 |

accuracy = 1.00

## Multimodal shares

| modality | group-permute | VIMP |
| --- | --- | --- |
| video | 0.000 | 0.492 |
| audio | 0.000 | 0.225 |
| text | 1.000 | 0.283 |

Mag×Var on PermuCATE RF unit risk (φ−τ̂)² resolves benign vs malignant shift (2×2); multimodal group-permute attributes the channel.

![](permucate_magvar_mm.png)

```bash
python3 scripts/run_as_dta_permucate_magvar.py
```


## Run (PermuCATE RF path)

```bash
python3 scripts/run_as_dta_permucate_magvar.py
```

Code: `Python/src/as_dta_permucate_magvar.py` (RF μ₀/μ₁ + π + τ̂; Mag×Var on (φ−τ̂)²; V/A/T group-permute).
