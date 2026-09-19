# iter09 — multi-seed / LOO-pos selection stability

**Framing (讲武德):** Freeze a feature set for DS handoff — still period-W PO, **not** ATE.

## Recipes

| tag | how the fixed set is built |
|---|---|
| `*__single` | per-seed selection (reference) |
| `*__maj2` | features in ≥2/3 seed selections |
| `*__intersect_fill` | seed intersection, fill to k by vote |
| `P_po_avgVIMP_topk` | seed-average PO-VIMP → top-k (no F reselect) |
| `P_po_avgVIMP_pool__maj2` | avg-VIMP pool → per-seed F → maj2 |
| `P_po__LOO_pos_maj` | leave-one-W1-positive-out ×4 → majority |

Eval: same fixed set, fit official FSDS on each seed’s W1-train → W2.

## Mean W2 (k=15, seeds 0/1/2)

| variant | mean W2 | σ |
|---|---:|---:|
| **P_po single** | **0.722** | 0.170 |
| Z_PO_rare single | 0.717 | 0.137 |
| A single | 0.716 | 0.154 |
| **P_po LOO-pos maj** | **0.708** | **0.114** |
| Z_PO_rare intersect/maj | ~0.68–0.69 | ~0.19 |
| P_po maj2 / intersect | ~0.67–0.68 | ~0.15 |
| seed-avg VIMP topk | 0.551 | 0.187 |
| A maj2 | 0.556 | 0.199 |

LOO raw intersection size: P_po seeds share **10/15**; Z_PO_rare **11/15**; A only **6/15**.

## Takeaway

1. **Don’t freeze seed-majority sets** as a free upgrade — maj2/intersect **lose mean** vs single-seed P_po.  
2. **LOO-pos majority** is the stability play: slightly lower mean (0.708) but **lowest σ (0.114)** — use when DS needs one frozen list.  
3. **Seed-avg PO-VIMP top-k alone** (skip FSDS reselect) **hurts** — always hand the pool to official FSDS.  
4. Production recipe stays: `po_help_select` / `P_po_vimp` per split, report seed mean±std; optional LOO freeze for docs.  
5. No ATE claim.

Artifacts: `stability_sweep.csv`, `stability_summary.csv`, `po_vimp_seed_avg.csv`, `fixed_*.csv`.
