# Token Mag×Var AUC — PermuCATE unit PO (benign vs malignant)

Token-wise perturb→re-forward + Mag×Var AUC on PermuCATE unit PO (Y=response, W=REF/LIVE) — benign vs malignant (no multimodal).

## Estimator

- **PermuCATE RF** on text: $Y$=response, $W$=REF/LIVE, $X$=TF-IDF→SVD
- Unit $u_i=(\varphi_i-\hat\tau_i)^2$ — **not** $1-\hat p(Y\mid x)$
- Token VAI: drop → re-featurize → frozen $\hat\mu_w,\hat\pi,\hat\tau$

## Logic

1. Y = response label; W = REF(0)/LIVE(1); X = TF-IDF→SVD
2. PermuCATE RF: μ̂0, μ̂1, π̂ → φ → τ̂; unit u=(φ−τ̂)²
3. Magnitude μ, Variance s of bootstrap sub-env mean u
4. Token VAI: drop j → re-featurize → frozen models → |Δu|
5. Malignant oracle: subset attack + Y flip (outcome break)
6. Severity = J + Y_flip + rateΔ; AUC malignant vs benign

## Pooled AUC

- **severity**: 1.000
- **J**: 0.432
- **d_s**: 0.420
- **d_mu**: 0.485
- **label_flip_frac**: 1.000
- **token_top_share**: 0.583
- **J_plus_token**: 0.475

## Per board AUC(J)

| board | AUC(severity) | AUC(J) | AUC(Y-flip) | AUC(token) |
| --- | --- | --- | --- | --- |
| hc3_all | 1.000 | 0.611 | 1.000 | 0.361 |
| hc3_finance | 1.000 | 0.417 | 1.000 | 0.639 |
| news20_sci | 1.000 | 0.250 | 1.000 | 0.639 |

## Justify

Magnitude alone confuses uniform synonym with attack when both raise μ. Variance of PermuCATE unit risk + token-attr concentration rises under malignant subset triggers; AUC of J / d_s justifies separation.

![](token_magvar_auc.png)

```bash
python3 scripts/run_as_dta_token_magvar_auc.py --quick
python3 scripts/run_as_dta_token_magvar_auc.py
```
