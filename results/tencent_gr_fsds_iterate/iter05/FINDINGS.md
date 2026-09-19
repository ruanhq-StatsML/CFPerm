# iter05 — PO-risk helps the data science loop

**Framing (讲武德):** `W` = period, not treatment. PO-risk = period-shift proxy  
`mean(τ̂²)`; fit **once**; feeds feature ranking for FSDS — **not** an ATE claim.

## What we shipped

Reusable helper: `scripts/tencent_gr/po_risk_fsds.py`

```text
fit_period_po / fit_po_on_windows  →  risk, τ̂², PO-VIMP
blend_cmean_po_scores              →  |δ| ⋈ VIMP ranking
```

Variants in the iterate harness:

| id | recipe |
|---|---|
| `P_po_vimp_FSDS` | PO-VIMP top-(k+3) → official FSDS |
| `Z_combined_PO` | cmean ⋈ PO-VIMP → π-stable F → FSDS |
| `Z_combined` | cmean → π → FSDS (no PO) |
| `A_baseline_F` | official FSDS only |

## Numbers (k=15, seeds 0/1/2)

| variant | mean W2 HGB | std | mean W2 AP |
|---|---:|---:|---:|
| A baseline FSDS | 0.716 | 0.154 | 0.0099 |
| Z combined | 0.720 | 0.157 | 0.0072 |
| **P PO-VIMP → FSDS** | **0.722** | 0.170 | 0.0070 |
| Z combined-PO | 0.711 | 0.171 | 0.0067 |

Per-seed W2 HGB:

| seed | A | P_po | Z | Z_PO |
|---:|---:|---:|---:|---:|
| 0 | 0.768 | 0.759 | **0.774** | 0.739 |
| 1 | 0.543 | 0.536 | 0.543 | 0.529 |
| 2 | 0.836 | **0.871** | 0.843 | 0.867 |

## PO feature story (seed 0)

Global PO-risk ≈ **1.15e-6**. Top VIMP: `ui_pop_mismatch` (**35.7%** share), then `u_n_exp`, covisit / item volume logs — shift-relevant graph feats a DS can read without claiming causality.

Artifacts: `po_feature_vimp_seed0.csv`, `po_risk_seed0.txt`, `sweep_seed_tight.csv`.

## Takeaway for every data science guy

1. Fit PO **once** on the localized support (period W).  
2. Use **PO-VIMP** (and/or τ̂²) as a **ranking prior** into FSDS.  
3. Keep official FSDS as the supervised gate (W1 fit; W2 holdout).  
4. Always average over seeds — PO does not remove rare-pos split noise.  
5. Do **not** say “PO proves treatment effects.”
