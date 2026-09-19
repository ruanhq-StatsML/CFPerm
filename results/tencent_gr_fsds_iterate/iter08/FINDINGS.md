# iter08 — k-sweep + τ̂² row filter

**Framing (讲武德):** `W` = period. τ̂² row filter = keep period-shift mass (+ all positives) before selection — **not** sample weights for ATE.

## What we ran

1. **k ∈ {10, 12, 15, 18}** on `A_baseline_F`, `P_po_vimp_FSDS`, `Z_combined_PO_rare`
2. **`P_po_tau2_rows`** at k=15: τ̂²≥median rows (always keep converts) → PO-VIMP → FSDS

Helper: `filter_by_tau2_quantile` / `tau2_on_frame` in `po_risk_fsds.py` (stores `tau_model`).

## Mean W2 HGB (seeds 0/1/2)

| variant | k=10 | k=12 | k=15 | k=18 |
|---|---:|---:|---:|---:|
| A baseline | 0.559 | 0.648 | **0.716** | 0.714 |
| **P_po_vimp** | 0.472 | 0.553 | 0.722 | **0.724** |
| Z_PO_rare | 0.424 | 0.562 | **0.717** (σ=0.137) | 0.713 |
| P_po_τ̂²_rows | — | — | 0.714 | — |

## Takeaway

1. **k &lt; 15 hurts** — especially PO-guided paths (pool too aggressive / under-spec).  
2. **Default k=15** still the sweet spot for rare-π; **P_po** can take **k=18** for a tiny mean bump (0.724).  
3. **τ̂² row filter** ≈ plain `P_po_vimp` (no clear win on this grid; seed 0/1 identical). Keep as optional DS knob, not default.  
4. Seed σ still ≫ method gap — always report mean±std.  
5. Do not claim ATE.

Artifacts: `k_tau2_sweep.csv`, `k_tau2_summary.csv`.
