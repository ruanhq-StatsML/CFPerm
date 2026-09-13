# Attribution → next-batch training adapter

FSDS/TSS hop `(ĉ, δ̂)` becomes **training weights** for the next Amazon batch. Only score: online rating MSE.

## Live numbers (4 seeds, `n_per=240`, `X` is `(2160, 128)`)

| method | online MSE | vs bank |
| --- | --- | --- |
| **hop_ridge** | **1.378 (0.043)** | **−0.188** |
| ridge_past | 1.402 (0.052) | −0.164 |
| **dga_ridge** | **1.417 (0.054)** | **−0.150** |
| attr_adapter | 1.447 (0.058) | −0.119 |
| bank | 1.566 (0.063) | 0 |
| plateau | 1.590 (0.068) | +0.024 |
| river_pa | 2.462 (0.084) | +0.896 |

Hop-weighted Ridge is the win. DGA (gradient-alignment domain weights) beats
bank/plateau and sits next to `ridge_past`, but does **not** beat hop cosine IW
on this board. Typed bank⊕Ridge mix still beats bank/plateau. `river` PA is the
honest package baseline — it loses on this sparse TF-IDF board.

### DGA drop-in (Fan–Grangier–Ablin 2024, arXiv:2410.02498)

Density-ratio / embedding IS is the weak link DGA was built for. Mapping:
past categories → domains, \(D_{\mathrm{spe}}=B_{t-1}\), mirror descent + EMA →
Ridge sample weights (`dga_ridge`). Same adapter slot as hop IW; no new loss.
Amazon champion stays `hop_ridge`. Reach for DGA when heatmap/IS weights are
unstable (tiny specialized pocket / overfit domains). Note:
`docs/method/DGA_hop_note.tex`.

## What the adapter does

1. **Heatmap → sample weights** (Shimodaira IW surrogate):  
   `w_s = exp(γ (cos(μ_s, μ_{t−1}) − 1))` on past categories, then closed-form Ridge.
2. **DGA → sample weights** (Fan et al. 2024): gradient alignment + mirror descent + EMA over past batches as domains; `dga_ridge`.
3. **TSS signs → mixture** (attr_adapter): large `ĉ` / quiet `δ̂` → Wu bank up; large `δ̂` → hop Ridge up. EWMA stacks the two predictors' recent MSE.
4. **Rolling hop stats**: per hop log `(ĉ, δ̂, w_mean, predictor MSE)` — the online window statistics that drive the gate.
5. **Modality / VIMP (MSR-VTT path)**: same gate reads per-block heatmap hops and RF-Domain shares as *which-head budget*; Amazon is the scalar-text MSE board for the weight logic.

Follow-up locked recipe: Amazon stays **`hop_ridge`**. Multimodal uses
**modality-π hop**
(`w_s=exp(γ Σ_m π_m (cos(μ_s^m,μ_{t-1}^m)-1))`, API
`attribution_adapter.modality_hop_weights`). PO-risk VIMP is for **mask
attribution** only (zero/noise/missing), not feature reweighting — see
`results/po_hop_attribution/`.

## Run

```
python3 scripts/run_attribution_adapter.py
python3 scripts/run_attribution_adapter.py --quick
```

Code: `Python/src/attribution_adapter.py`  
Docs: `docs/method/Attribution_adapter_note.tex`, `docs/method/DGA_hop_note.tex`
