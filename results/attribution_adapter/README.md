# Attribution → next-batch training adapter

FSDS/TSS hop `(ĉ, δ̂)` becomes **training weights** for the next Amazon batch. Only score: online rating MSE.

## Live numbers (4 seeds, `n_per=240`, `X` is `(2160, 128)`)

| method | online MSE | vs bank |
| --- | --- | --- |
| **hop_ridge** | **1.378 (0.043)** | **−0.188** |
| ridge_past | 1.402 (0.052) | −0.164 |
| attr_adapter | 1.447 (0.058) | −0.119 |
| bank | 1.566 (0.063) | 0 |
| plateau | 1.590 (0.068) | +0.024 |
| river_pa | 2.462 (0.084) | +0.896 |

Hop-weighted Ridge is the win. Typed bank⊕Ridge mix still beats bank/plateau; pure bank is not enough once you reuse past batches with heatmap IW. `river` PA is the honest package baseline — it loses on this sparse TF-IDF board.

## What the adapter does

1. **Heatmap → sample weights** (Shimodaira IW surrogate):  
   `w_s = exp(γ (cos(μ_s, μ_{t−1}) − 1))` on past categories, then closed-form Ridge.
2. **TSS signs → mixture** (attr_adapter): large `ĉ` / quiet `δ̂` → Wu bank up; large `δ̂` → hop Ridge up. EWMA stacks the two predictors' recent MSE.
3. **Rolling hop stats**: per hop log `(ĉ, δ̂, w_mean, predictor MSE)` — the online window statistics that drive the gate.
4. **Modality / VIMP (MSR-VTT path)**: same gate reads per-block heatmap hops and RF-Domain shares as *which-head budget*; Amazon is the scalar-text MSE board for the weight logic.

Follow-up: PO-risk feature weights do **not** beat `hop_ridge` on Amazon; use PO for token/patch **mask attribution** and modality-π hop on MSR-VTT — see `results/po_hop_attribution/`.

## Run

```
python3 scripts/run_attribution_adapter.py
python3 scripts/run_attribution_adapter.py --quick
```

Code: `Python/src/attribution_adapter.py`  
Note: `docs/method/Attribution_adapter_note.tex`  
Figure: `attribution_adapter.png`
