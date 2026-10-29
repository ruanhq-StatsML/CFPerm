# Locked recipe: heatmap hop + modality-π + PO mask

## What to run

| setting | method | role |
| --- | --- | --- |
| **Amazon (text)** | **`hop_ridge`** | heatmap IW → closed-form Ridge — best online MSE (**1.378**) |
| **MSR-VTT (multi)** | **modality-π hop** | \(w_s=\exp(\gamma\sum_m\pi_m(\cos(\mu_s^m,\mu_{t-1}^m)-1))\) |
| **PO-risk VIMP** | **mask only** | zero / noise / missing on top-PO tokens or patches → ΔMSE |

Do **not** scale features by PO VIMP for training (amplify or downweight both lose to `hop_ridge`).

## Modality-π (contribution)

\(\pi_m\) from ``benchmark_feature_selection`` → ``modality_mass(vimp)``
(default = RF-Domain VIMP; swap the hook body to change selector).
Optional EWMA across hops (`ewma_pi≈0.3`). Causal reference is still \(\mu_{t-1}\).

Live MSR-VTT (video+audio → text): hop 0.0898 → hop+π_m **0.0896**; shares video≈0.71, audio≈0.26.

## Amazon mask attribution (top-16 PO)

| mask | ΔMSE |
| --- | --- |
| zero / missing | +0.126 |
| noise | +0.144 |

## Attribution → next-batch ensemble EWMA

`(ĉ, δ̂, Δπ)` → `attribution_ewma_rate` → λ_new for bank⊕Ridge (Amazon
`attr_adapter`) and for π persistence (MSR-VTT `ewma_pi="auto"`).

- covariate hop → small λ_new → keep past ensemble
- concept hop / π jump → large λ_new → refresh mix

See [`TUNING.md`](TUNING.md).

## Run

```bash
PYTHONPATH=Python/src python3 scripts/run_po_hop_attribution.py
PYTHONPATH=Python/src python3 scripts/run_attribution_adapter.py
```

API: `attribution_adapter.attribution_ewma_rate`, `modality_pi_shares`  
Tuning: [`TUNING.md`](TUNING.md)
