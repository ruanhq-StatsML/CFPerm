# Hierarchical dissolve for multi-level drill — landability note

**讲武德:** `W`=period. Scores are shift proxies — **not** ATE.

## Naming (important)

| name | what it is |
|---|---|
| PyPI [`dissolve`](https://pypi.org/project/dissolve/) | API-deprecation migrator (`@replace_me`) — **wrong tool for attribution** |
| **This spike** | geopandas-style **rollup**: leaf score → `groupby(parent).agg` up `order → user → merchant` |

## What we tried

```text
W2 leaf score (shift_l2 or τ̂²)
    └─ dissolve → merchant ranking   vs  top-down MMD merchant
    └─ dissolve → user ranking       vs  MMD users (inside top merchants)
    └─ cascade: top-M dissolve merchants → dissolve users inside
```

Code:
- `scripts/tencent_gr/attribution_dissolve.py`
- `scripts/tencent_gr/run_drill_dissolve_smoke.py`
- Results: `results/tencent_gr_drill_dissolve/`

## Smoke numbers (localized grids)

**Merchant Jaccard (dissolve `shift_l2` vs MMD):**

| K | Jaccard |
|---:|---:|
| 5 | **0.67** |
| 10 | 0.54 |
| 20 | 0.60 |
| 50 | **0.79** |

**User Jaccard:** global dissolve ≈ **0**; cascade ≈ **0–0.03** vs MMD-in-merchants.  
Bottom-up ∑‖x−μ‖₂ ≠ entity-conditional MMD — different objective.

**τ̂² dissolve vs MMD merchants:** Jaccard ≈ 0–0.05 — poor match on this grid.

**FSDS-on-support AUC:** top-50 user supports have **0 W2 converts** → eval falls back; not informative under rare-pos.

## Landability verdict

| Question | Answer |
|---|---|
| Can we implement cheaply? | **Yes** — pandas only; plugs into existing W1/W2 grids |
| Replace three-step MMD drill? | **No** — user-level rankings diverge; keep MMD/cmean as the drill gate |
| Useful role? | **Yes as merchant prior / candidate gen** — shift_l2 dissolve tracks MMD tops well |
| Ship PyPI `dissolve`? | **No** |
| Claim ATE? | **No** |

### Recommended landing shape

```text
1. leaf shift_l2 on W2
2. dissolve → merchant TopM   (cheap prior; high overlap with MMD)
3. optional: cascade dissolve users inside TopM  (diagnostic only)
4. official path: cmean / MMD drill + stop table + FSDS / po_help_fsds
```

Do **not** replace the locked drill-stop protocol with dissolve alone.
