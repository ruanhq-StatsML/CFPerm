# iter07 — α ablate + PO-VIMP×rare-π hybrid

**Framing (讲武德):** `W` = period. PO-risk = shift proxy. Not an ATE.

## Variants

| id | recipe |
|---|---|
| `P_po_vimp_FSDS` | PO-VIMP top-(k+3) → FSDS (no π) |
| `Z_combined_PO` / `_a03` / `_a07` | cmean⋈PO α∈{0.5,0.3,0.7} → rare-capped π → FSDS |
| `P_po_vimp_rare_pi` | PO-VIMP pool → rare-pos π → FSDS (no \|δ\|, no boot) |
| `Z_combined_PO_rare` | po_help (α=0.5) → rare π → FSDS |

## Numbers (k=15, seeds 0/1/2, pool=k+3)

| variant | mean W2 | std |
|---|---:|---:|
| **P_po_vimp_FSDS** | **0.722** | 0.170 |
| A baseline | 0.716 | 0.154 |
| Z_PO α=0.3/0.5/0.7 | **0.717** | **0.137** |
| P_po_vimp_rare_π | 0.717 | 0.137 |
| Z_combined_PO_rare | 0.717 | 0.137 |

α∈{0.3,0.5,0.7} are **bit-identical** on this grid (same selected sets / same W2).

## Tight-pool secondary (`alpha_pool_sweep.csv`)

| pool | mean W2 (Z_PO any α) |
|---|---:|
| k+3 (=18/21) | 0.717 |
| tight k (=15/21) | **0.556** (hurts) |

α still identical under tight pool — blend reorder does not change Top-15 enough on d=21.

## Takeaway

1. **Keep α=0.5** as the DS default; no need to tune α on this feature width.  
2. **Slack matters:** pool = k+3 (not k) so official FSDS can still move; tight k collapses W2.  
3. **P_po_vimp** (PO → FSDS, skip π) still best **mean**; rare-π family best **σ**.  
4. Hybrid `P_po_vimp_rare_pi` ≡ rare-π blended family here — π step dominates once the pool is near-full.  
5. Do not claim ATE.

Artifacts: `variant_compare.csv`, `alpha_pool_sweep.csv`, `sweep` via `iter07` / `_s1` / `_s2`.
