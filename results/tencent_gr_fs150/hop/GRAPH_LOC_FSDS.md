# Recsys: graph localization → multi-layer FSDS

Tencent-GR. Graph first, then FSDS only drills what the graph localized.
PO-risk is **not** causal. RF-domain is a **portrait log**, not the score to optimize.

## Protocol

1. Left-window clk/cnv → user—item, user—merchant, session co-click, shop projection.
2. **Localize** (no Y): each graph family vs early/late clock W. Moved = RF-domain AUC ≥ 0.55.
3. **Drill** with Y=`y_post_clk_1d`:
   - L1 hops (funnel ∪ graph) LOGO
   - L2 families inside the localized graph hop
   - L3 LOCO columns inside that hop

n=12564 p=18. Clock = median `t_end`. W=1 later half.

## 1. Graph localization (no Y)

| family | n | RF-domain AUC | moved |
|---|---:|---:|---|
| user_connectivity | 2 | 0.605 | yes |
| item_connectivity | 2 | 0.523 | no |
| merchant_structure | 4 | 0.501 | no |

Localized families: `user_connectivity`.
Moved = this block's connectivity portrait differs early vs late. Not “this tower caused conversion”.

## 2. FSDS L1 — hops (LOGO)

RF-domain (all X) **0.781** · PO-risk **0.000178**

| hop | n | RF-mass | PO-mass | LOGO Δ | share |
|---|---:|---:|---:|---:|---:|
| graph_item | 2 | 0.004 | 0.044 | +0.00003 | 0.408 |
| funnel_order | 5 | 0.162 | 0.173 | +0.00002 | 0.353 |
| graph_merchant | 4 | 0.020 | 0.212 | +0.00002 | 0.240 |
| funnel_user | 5 | 0.662 | 0.521 | -0.00024 | 0.000 |
| graph_user | 2 | 0.152 | 0.051 | -0.00001 | 0.000 |

LOGO Δ>0: dropping the hop **lowers** R (tied to the early/late Y gap).
Δ≤0: this hop is not a concept-gap source on this clock. Read RF-mass for P(X).

Drill continues inside localized hops: `graph_user`.

On this clock that pairing is the point: **user_connectivity moved**, but L1 LOGO Δ on `graph_user` is ≤0.
Portrait movement ≠ Y-gap source. Funnel_user has the RF-mass (later people look different) and also LOGO Δ<0.
Small LOGO+ on `graph_item` / `funnel_order` / `graph_merchant` sits on PO-risk ~1e-4 — log, not a story.

## 3. FSDS L2 — families inside localized hops

Only one family inside the localized hop, so L2 LOGO is vacuous (nothing to drop).
L3 LOCO on those columns is the drill.

## 4. FSDS L3 — LOCO on localized hop columns

| feat | hop | LOCO ΔR |
|---|---|---:|
| `g_u_item_deg` | graph_user | +1.655635e-05 |
| `g_u_merch_deg` | graph_user | +1.209499e-05 |

LOCO ΔR>0: dropping the column lowers R. Still not a treatment effect.

## Read

- Localization answers **which graph block's portrait moved**.
- L1–L3 answer **among those blocks, what is tied to the Y-gap** (if anything).
- Funnel hops stay on the L1 board so graph is compared, not isolated.
- Do not turn LOGO share into a unique importance ranking.

`PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds_drill.py`
