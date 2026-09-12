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

## System package (ML-infra)

Amazon is a **smoke client**. Controllers live in `agod/`:

| module | role |
|---|---|
| `agod.mmd` | bounded-cost MMD² sensor for \(P(X)\) |
| `agod.shift` | concept / covariate decomposition |
| `agod.lr_controller` | Softmax α → per-modality LR (actuator + EMA) |
| `agod.policies` | B1–B5 named policies; **default Acc = B5** |

Harness (IO/model only): `scripts/run_agod_amazon_mmd_lr.py`  
Infra write-up: [`docs/agod/AGOD_ml_infra_justification.md`](../../docs/agod/AGOD_ml_infra_justification.md)

```bash
PYTHONPATH=. python3 tests/test_agod_controller.py
PYTHONPATH=. python3 scripts/run_agod_amazon_mmd_lr.py
```

### Shipped Acc policy (B5)

```
covariate_m = MMD²(X_m^ref, X_m^cur)     # damp LR
concept_m   = PO_m                       # raise LR
score_m     = λ_c·z(concept) − λ_v·z(cov)
LR_m ∝ Softmax(score / τ)
```

| Policy | Mean Acc lift | vs B1 |
|---|---:|---:|
| B1 equal | +0.058 | — |
| B2 RF | +0.072 | +0.014 |
| B3 pure MMD | +0.058 | ~0 |
| B4 MMD+gain | +0.037 | −0.021 |
| **B5 MMD-cov+PO** | **+0.087** | **+0.030** |

## Legacy RF concept/cov script

```bash
python3 scripts/run_agod_amazon_concept_cov_lr.py
python3 scripts/run_agod_amazon_modality_lr.py
```

HF dataset: [`jingxiang11111/amazon_reviews_for_rec`](https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec)
