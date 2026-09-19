# PO-risk for every data science guy

讲武德 framing: **`W` = period (W1/W2), not treatment.**  
PO-risk = period-shift proxy `mean(τ̂²)`. Fit **once**. Feeds FSDS ranking — **not** an ATE claim.

## CLI (recommended)

```bash
PYTHONPATH=. python3 scripts/tencent_gr/po_help_fsds.py \
  --w1-grid results/tencent_gr_localize_fsds_time/W1_localized_grid.parquet \
  --w2-grid results/tencent_gr_localize_fsds_time/W2_localized_grid_sample.parquet \
  --select-k 15 --seed 0 \
  --out-dir results/tencent_gr_fsds_iterate/po_help_cli
```

## One-liner (Python)

```python
from po_risk_fsds import po_help_select

report = po_help_select(g_w1, g_w2, cols, k=15, seed=0)
# report["pool"]           → hand to official FSDS (Scaler→Var→SelectKBest)
# report["feature_table"]  → PO-VIMP ranks for notebooks
# report["blend_table"]    → |δ| ⋈ VIMP
# report["risk"]           → shift proxy (compare supports; not an effect size)
print(report["note"])
```

## Recipe (locked)

```
1. Fit period-PO once on localized support (W1-scaled X; W=period)
2. Rank by PO-VIMP and/or blend with cmean |δ|
3. Optional: rare-pos-capped π on the guided pool (lower σ)
4. Official FSDS on that pool (W1 fit only; W2 = temporal holdout)
5. Average over seeds — rare convert makes split noise ≫ method gap
```

Overnight evidence: `results/tencent_gr_fsds_iterate/ITERATION_LOG.md` (OVERNIGHT_SUMMARY).

## Helpers

| function | job |
|---|---|
| `fit_period_po` / `fit_po_on_windows` | risk, τ̂², VIMP |
| `blend_cmean_po_scores` | \|δ\| ⋈ VIMP |
| `bootstrap_pi_select` | stratified bootstrap π (rare-pos friendly) |
| `rare_pos_n_splits` | cap CV folds by `#positives` |
| `po_help_select` | DS entrypoint → `pool` + tables |

## Iterate variants

| id | recipe |
|---|---|
| `P_po_vimp_FSDS` | PO-VIMP top-(k+3) → FSDS |
| `Z_combined_PO` | cmean⋈PO → 5-fold π → FSDS |
| `P_po_boot_pi` | PO-help pool → bootstrap π → FSDS |
| `Z_combined_PO_rare` | cmean⋈PO → rare-pos-capped π → FSDS |

## Blend α

Default **α=0.5** (`|δ| ⋈ PO-VIMP`). On the current graph-feat width (d≈21, pool=k+3),
α∈{0.3,0.5,0.7} is a no-op — keep 0.5; do not over-tune. Prefer pool slack (k+3) over tight k.


## k and τ̂² rows

- Prefer **k=15** (or **18** for plain PO-VIMP). k∈{10,12} underperforms on this grid.
- Optional `filter_by_tau2_quantile` keeps high-τ̂² rows (+ all positives) before FSDS — optional knob; not better than plain PO-VIMP by default.


## Freezing a feature list

- Prefer **per-split** `P_po_vimp` + seed mean±std.
- If you need one frozen list: **leave-one-positive-out majority** beats seed-maj2 on σ (iter09). Do **not** ship seed-average VIMP top-k without FSDS reselect.

## Do / Don't

- Do use PO-VIMP as a **ranking prior** into FSDS  
- Do report seed mean±std under rare positives  
- Don't say “PO proves treatment effects”  
- Don't replace official FSDS with PO alone  
