# FSDS × Interaction Value Online Adaptation

Localization names the sick modality; **clever covariates** treat at
inference without retrain:

```
score_adapted = score_raw / (1 + λ · W_shift)
```

W_shift mixes PO-risk, triangle inconsistency, and Interaction-Value drop.
`fsds_clever_cov` gates λ with typed hop `(ĉ, δ̂)` (locked TSS signs).

## Live numbers

| scenario | static | localize_only | clever_cov | **fsds_clever_cov** | lift vs localize |
| --- | ---: | ---: | ---: | ---: | ---: |
| healthy | 0.2597 | 0.2597 | 0.2567 | **0.2583** | 0.0014 |
| interaction_break | 0.3289 | 0.3289 | 0.3184 | **0.3095** | 0.0194 |
| alignment_attack | 0.3197 | 0.3197 | 0.2991 | **0.2928** | 0.0270 |
| concept_drift | 0.6768 | 0.6768 | 0.4968 | **0.4891** | 0.1877 |
| compound | 0.3418 | 0.3418 | 0.3091 | **0.3036** | 0.0383 |

## What each method does

- `static` — joint Ridge, no gate.
- `localize_only` — still computes π / Φ / triangle (diagnosis) but leaves preds alone.
- `clever_cov` — KS on PO-risk & alignment → inject W into scores.
- `fsds_clever_cov` — same W, λ from `(ĉ,δ̂,Δπ)` + triangle/Φ collapse.
- `oracle_penalty` — upper bound using true |error| as W.

## MSR-VTT hop probe

Source: `msrvtt` · n0=160 n1=160 · RF AUC=0.746

| π_text | π_video | π_audio | mean Φ_TI | triangle | mse_joint |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.107 | 0.634 | 0.259 | 0.0119 | 0.457 | 1.7146 |

## Run

```
python3 scripts/run_fsds_interaction_online_adapt.py
python3 scripts/run_fsds_interaction_online_adapt.py --quick
```

Code: `Python/src/fsds_interaction_online_adapt.py`

Product line: Pseudo-Outcome Risk / Interaction Shift as **clever covariates**
closes the multimodal loop — diagnose with FSDS, treat in milliseconds.
