# Locked recipe: heatmap hop + modality-π + PO mask

## What to run

| setting | method | role |
| --- | --- | --- |
| **Amazon (text)** | **`hop_ridge`** | heatmap IW → closed-form Ridge — best online MSE (**1.378**) |
| **MSR-VTT (multi)** | **modality-π hop** | \(w_s=\exp(\gamma\sum_m\pi_m(\cos(\mu_s^m,\mu_{t-1}^m)-1))\) |
| **PO-risk VIMP** | **mask only** | zero / noise / missing on top-PO tokens or patches → ΔMSE |

Do **not** scale features by PO VIMP for training (amplify or downweight both lose to `hop_ridge`).

## Modality-π (contribution)

\(\pi_m\) = PO-risk VIMP mass on block \(m\), blended 50/50 with RF-Domain share on the same causal hop. Causal reference is still \(\mu_{t-1}\).

Live MSR-VTT (video+audio → text): hop 0.0898 → hop+π_m **0.0896**; shares video≈0.71, audio≈0.26.

## Amazon mask attribution (top-16 PO)

| mask | ΔMSE |
| --- | --- |
| zero / missing | +0.126 |
| noise | +0.144 |

## Run

```bash
PYTHONPATH=Python/src python3 scripts/run_po_hop_attribution.py
PYTHONPATH=Python/src python3 scripts/run_attribution_adapter.py
```

API: `attribution_adapter.modality_pi_shares` → `benchmark_feature_selection`  
Tuning: see [`TUNING.md`](TUNING.md)
