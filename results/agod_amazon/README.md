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

## Run

```bash
# put shards under data/amazon_reviews/shards/ (from HF), then:
python3 scripts/run_agod_amazon_modality_lr.py
```

HF dataset: [`jingxiang11111/amazon_reviews_for_rec`](https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec)
