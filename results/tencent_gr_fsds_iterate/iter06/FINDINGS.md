# iter06 — rare-pos bootstrap π + PO-help for every DS guy

**Framing (讲武德):** `W` = period. PO-risk = `mean(τ̂²)` shift proxy.  
Fit once → ranking prior → official FSDS. **Not** an ATE.

Cookbook: [`docs/tencent_gr/PO_RISK_FOR_DS.md`](../../../docs/tencent_gr/PO_RISK_FOR_DS.md)

## What we shipped

| API | job |
|---|---|
| `po_help_select` | DS one-liner → `pool` + VIMP + |δ|⋈VIMP |
| `bootstrap_pi_select` | stratified bootstrap π (rare-pos friendly) |
| `rare_pos_n_splits` | cap CV folds by `#positives` |

New variants:

| id | recipe |
|---|---|
| `P_po_boot_pi` | PO-help pool → bootstrap π → FSDS |
| `Z_combined_PO_rare` | cmean⋈PO → rare-pos-capped π → FSDS |

Also: default `build_stable_pi_f` now uses `rare_pos_n_splits` (stops 5-fold warn when n_pos=3).

## Numbers (k=15, seeds 0/1/2)

| variant | mean W2 HGB | std | mean W2 AP |
|---|---:|---:|---:|
| **P PO-VIMP → FSDS** | **0.722** | 0.170 | 0.0070 |
| Z combined | 0.720 | 0.157 | 0.0072 |
| **Z combined-PO-rare** | **0.717** | **0.137** | 0.0054 |
| A baseline FSDS | 0.716 | 0.154 | 0.0099 |
| Z combined-PO (5-fold) | 0.711 | 0.171 | 0.0067 |
| P PO-boot-π | 0.701 | 0.156 | 0.0066 |

Per-seed W2 HGB:

| seed | A | P_po | Z | Z_PO | P_boot | Z_PO_rare |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.768 | 0.759 | **0.774** | 0.739 | 0.668 | 0.737 |
| 1 | 0.543 | 0.536 | 0.543 | 0.529 | 0.565 | **0.571** |
| 2 | 0.836 | **0.871** | 0.843 | 0.867 | 0.871 | 0.843 |

## Takeaway

1. **Default DS path:** `po_help_select` → hand `pool` to official FSDS (`P_po_vimp_FSDS` still best mean).  
2. **Variance story:** `Z_combined_PO_rare` cuts seed σ (0.137 vs ~0.15–0.17) and lifts the bad seed (0.571 vs ~0.54) — use when rare-pos split noise dominates.  
3. Bootstrap π alone (`P_po_boot_pi`) is **not** a free lunch (hurts seed 0); keep as ablation.  
4. Always report seed mean±std. Do **not** claim ATE.

Artifacts: `sweep_seed.csv`, `sweep_seed_summary.csv`, `variant_compare.csv`.
