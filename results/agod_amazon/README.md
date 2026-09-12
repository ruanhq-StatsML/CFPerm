# AGOD Amazon — modality MSG → per-modality LR

## Critical: not inference latency

Savings are **online adaptation FLOPs** (skip modality projection grads when
`α_m < θ`). Serving still runs text+image encoders.

## Logic upgrade vs fused-dim masking

| Draft sketch | This prototype |
|---|---|
| Drift on fused vector | MSG on modality blocks (`text_proj` / `image_proj`) |
| Mask fusion dimensions | Rescale **per-modality LR** + hard-gate grads |
| Assumed `jpg`+`txt` | Real schema: `item/user.json` + `patch.bin` + `label.json` |

## Control law

```
g_m = Normalize(AUC_m · VIMP_m + γ · PO_m)
α   = Softmax(g / τ)          # EMA-smoothed
LR_m = lr0 · (β + (1-β) · α_m · |M|)
if α_m < θ: zero ∂L/∂θ_m      # adaptation FLOPs ↓
```

## Smoke result (3 local shards, category stream)

| Policy | Rel. adapt FLOPs | Mean loss | Mean α_text | LR× text/image |
|---|---:|---:|---:|---:|
| B1 static | 1.00 | 0.692 | 0.50 | 1.00 / 1.00 |
| B2 AUC-only | 1.00 | 0.692 | 0.48 | 0.97 / 1.03 |
| B3 AGOD-LR | **0.50** | 0.692 | 0.30 | **0.66 / 1.34** |

Category shift from Tools & Home Improvement → Sports/Fashion/… is attributed
mostly to **image**; AGOD raises image LR and gates text updates.

## Accuracy rule: concept↑ LR / covariate↓ LR

```
covariate_m = max(AUC_m − 0.5, 0) · (1 + VIMP_m)   # P(X) shift → damp LR
concept_m   = PO_m                                   # P(Y|X) shift → raise LR
score_m     = λ_c · concept_m − λ_v · covariate_m
α           = Softmax(score / τ)   (EMA)
LR_m        = lr0 · (β + (1−β) · α_m · |M|)
```

Smoke (3 shards, 6 category windows), mean held-out Acc lift:

| Policy | Rule | Mean Acc lift | Mean Acc post |
|---|---|---:|---:|
| B1 | equal LR | +0.058 | 0.539 |
| B2 | high total shift → high LR | +0.091 | 0.556 |
| **B3** | **concept↑ / covariate↓** | **+0.072** | **0.549** |

B3 − B1 lift = **+0.014**; B3 − B2 lift = −0.019 on this smoke.

```bash
python3 scripts/run_agod_amazon_concept_cov_lr.py
```

## Run (base modality-LR prototype)

```bash
# put shards under data/amazon_reviews/shards/ (from HF), then:
python3 scripts/run_agod_amazon_modality_lr.py
```

HF dataset: [`jingxiang11111/amazon_reviews_for_rec`](https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec)
