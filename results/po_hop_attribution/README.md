# PO-risk × heatmap hop-ridge

Iteration after the attribution-adapter board: keep **heatmap hop-ridge** as the
Amazon champion, then ask whether PO-risk feature selection and modality-specific
hop weights help further — and drill down with token/patch masks.

## Verdict (live)

| Amazon method (4 seeds, `n_per=240`, TF-IDF 128) | online MSE |
| --- | --- |
| **hop_ridge** (heatmap IW) | **1.378 (0.043)** |
| ridge_past | 1.402 (0.052) |
| hop_po_amp (amplify high-PO coords) | 1.448 (0.045) |
| hop_po_feat / hop_po_both (downweight high-PO) | 1.551 (0.065) |
| hop_po_mod (single text block ≈ hop_ridge) | 1.378 (0.043) |

**PO feature reweighting does not beat hop_ridge on Amazon.** PO-risk VIMP marks
*domain / hop* coordinates, not rating predictors; scaling features by that VIMP
(either amplify or TSS-style downweight) raises MSE. The useful PO use here is
**selection for mask attribution**, not a new training loss.

### Amazon token mask (top-16 PO coords, hop-ridge base)

| mask | ΔMSE vs unmasked |
| --- | --- |
| zero / missing (col-mean) | **+0.126 (0.028)** |
| noise | **+0.144 (0.026)** |

Masking top-PO tokens *hurts* online MSE → those coordinates carry predictive
mass under the hop-ridge fit. That is the drill-down attribution.

### MSR-VTT modality-π hop (cross-modal: video+audio → text)

Video-id one-hot is saturated on this 16-video head (acc≈1). Probe is online
Ridge reconstruction of the text block from video+audio.

| | online MSE |
| --- | --- |
| hop (uniform modality cosine) | 0.0898 |
| hop + π_m (PO⊕RF shares) | **0.0896** (lift +0.0002) |
| shares (last hop) | video ≈ 0.71, audio ≈ 0.26, text ≈ 0.03 |

Patch mask ΔMSE (top-32 PO coords in each *input* modality):

| modality | zero | noise |
| --- | --- | --- |
| video | +0.0009 | +0.0031 |
| audio | −0.0004 | +0.0017 |

Small effects on this head; directionally, video patches matter more under noise.

## Modality-specific hop weights (important)

Uniform heatmap hop:

\[
w_s = \exp\bigl(\gamma(\cos(\mu_s,\mu_{t-1})-1)\bigr)
\]

Modality-specific (this module):

\[
w_s = \exp\Bigl(\gamma \sum_m \pi_m \bigl(\cos(\mu_s^{(m)},\mu_{t-1}^{(m)})-1\bigr)\Bigr)
\]

with \(\pi_m\) = normalized PO-risk VIMP mass on block \(m\), blended 50/50 with
RF-Domain shares on the same causal hop. Causal reference remains \(\mu_{t-1}\).
On Amazon (one text block) this collapses to hop_ridge. On MSR-VTT it reweights
past batches by *which modality* the hop is using — the natural place for
heatmap + multimodal contribution.

## Run

```bash
PYTHONPATH=Python/src python3 scripts/run_po_hop_attribution.py
PYTHONPATH=Python/src python3 scripts/run_po_hop_attribution.py --quick
```

Code: `Python/src/po_hop_attribution.py`  
Artifacts: `po_hop_attribution.json`, `po_hop_attribution.png`, `PO_hop_attribution.tex`
