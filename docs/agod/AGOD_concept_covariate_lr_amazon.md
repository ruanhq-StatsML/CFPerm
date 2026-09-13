# Concept↑ LR↑ / Covariate↑ LR↓ on Amazon Reviews

Accuracy-focused AGOD control on `jingxiang11111/amazon_reviews_for_rec`.

## Rule

Per modality \(m\), window \(t\) vs reference:

| Signal | Formula | LR effect |
|---|---|---|
| Covariate (\(P(X)\)) | \(\max(\mathrm{AUC}_m-0.5,0)\cdot(1+\mathrm{VIMP}_m)\) | **lower** LR |
| Concept (\(P(Y\|X)\)) | \(\mathrm{PO}_m\) | **raise** LR |

\[
\mathrm{score}_m=\lambda_c\cdot\mathrm{concept}_m-\lambda_v\cdot\mathrm{covariate}_m,\quad
\alpha=\mathrm{Softmax}(\mathrm{score}/\tau),\quad
\mathrm{LR}_m=\mathrm{lr}_0\big(\beta+(1-\beta)\alpha_m|M|\big)
\]

## Policies

- **B1**: equal LR
- **B2**: old AGOD — high total MSG → high LR
- **B3**: this rule — concept↑ / covariate↓

## Smoke result (3 local shards, 6 category windows)

| Policy | Mean Acc lift | Mean Acc post | Wins \(>0\) |
|---|---:|---:|---:|
| B1 | +0.058 | 0.539 | 3/6 |
| B2 | +0.091 | 0.556 | 4/6 |
| **B3** | **+0.072** | **0.549** | **4/6** |

- **B3 − B1** mean Acc lift = **+0.014**
- **B3 − B2** mean Acc lift = **−0.019**

On this stream, image carries most concept mass → B3 keeps image LR high (~1.4–1.7×) and damps text (~0.3–0.6×). That beats equal LR; total-shift→high-LR (B2) is still slightly stronger on Acc lift here.

## Run

```bash
python3 scripts/run_agod_amazon_concept_cov_lr.py
```

Artifacts: `results/agod_amazon/agod_amazon_concept_cov_lr.json`,
`AGOD_Amazon_ConceptCov_LR_Acc_Dashboard.png`.
