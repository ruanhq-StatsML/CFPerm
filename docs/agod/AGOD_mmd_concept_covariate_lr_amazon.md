# MMD-informed concept↑ / covariate↓ LR on Amazon

## Concrete takeaway

Yes — MMD logic can lift Acc, but **how** you assign concept vs covariate matters.

| Signal | Estimator | LR effect |
|---|---|---|
| Covariate `P(X)` | unbiased RBF **MMD²(X)** | ↓ LR |
| Concept `P(Y\|X)` | **PO** residual contrast (FSDS) | ↑ LR |

Pure MMD joint-excess as “concept” (B3) was too noisy on small windows.
**B5 = MMD covariate + PO concept** (FSDS coupling) is the Acc-winning rule.

## Why MMD helps the LR controller

1. **RF AUC×VIMP saturates.** After the first category jump, both modalities often look equally domain-separable → cov≈0.5/0.5, so “covariate↓LR” barely fires.
2. **MMD² stays graded.** Per-modality RKHS distance keeps moving across the stream (text cov often 0.4–0.8), so damping is informative.
3. **Concept still needs PO.** Joint MMD excess / residual on n≈60–80 is high-variance; PO is a stabler `P(Y|X)` proxy (same FSDS pairing: MMD for `P(X)`, PO for concept).
4. **Intensity gain alone (B4) overshoots.** Multiplying LR by `(1+κ·relu(con−cov))` on top of Softmax hurt Acc lift here.

## Smoke Acc (6 category windows)

| Policy | Rule | Mean Acc lift | Mean Acc post |
|---|---|---:|---:|
| B1 | equal LR | +0.058 | 0.539 |
| B2 | RF AUC/PO con−cov | +0.072 | 0.549 |
| B3 | pure MMD con−cov | +0.058 | 0.546 |
| B4 | MMD + intensity gain | +0.037 | 0.548 |
| **B5** | **MMD cov + PO concept** | **+0.087** | **0.559** |

- **B5 − B1** Acc lift = **+0.030**
- **B5 − B2** Acc lift = **+0.016**
- **B5 − B3** Acc lift = **+0.030**

Control law (B5):

```
covariate_m = MMD²(X_m^ref, X_m^cur)
concept_m   = PO_m
score_m     = λ_c · z(concept) − λ_v · z(covariate)
α           = Softmax(score / τ)   (EMA)
LR_m        = lr0 · (β + (1−β)·α_m·|M|)
```

```bash
python3 scripts/run_agod_amazon_mmd_lr.py
```
