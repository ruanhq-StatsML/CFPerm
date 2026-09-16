# Recsys: graph localization → multi-layer FSDS

Tencent-GR. Graph first (plus fake video/audio towers attached on shops),
then FSDS only drills what localization marked as moved.
PO-risk is **not** causal. RF-domain is a **portrait log**, not the score to optimize.
Video/audio are a messy DGP, not production embeddings.

## Protocol

1. Left-window clk/cnv → user—item, user—merchant, session co-click, shop projection.
2. Fake media DGP on the same graph:
   - `10000` video_ids ~ `randint(0, 1e10)`; `64`-d Student-t embeddings
   - `10000` audio_ids the same way; `64`-d Uniform[-1, 1]
   - Messy attach on each merchant: Zipf subset + 16 themes + t-bias (no W/Y coupling)
   - User mean-pool of merchants touched in the left window (visit-weighted)
3. **Localize** (no Y): each graph family ∪ `video_tower` ∪ `audio_tower` vs clock W.
   Moved = RF-domain AUC ≥ 0.55.
4. **Drill** with Y=`y_post_clk_1d`:
   - L1 hops (funnel ∪ graph ∪ towers) LOGO
   - L2 families inside the localized hops
   - L3 LOCO columns inside those hops (top-8 if a 64-d tower)

n=12564 p=146. Clock = median `t_end`. W=1 later half.
Shops attached=400 · users pooled=3808.
User-pool zero rate video early/late 0.035/0.047; mean L2-norm 2.79/3.07. Not a missingness artifact.

## 1. Localization (no Y)

| family | n | RF-domain AUC | moved |
|---|---:|---:|---|
| audio_tower | 64 | 0.853 | yes |
| video_tower | 64 | 0.829 | yes |
| user_connectivity | 2 | 0.605 | yes |
| item_connectivity | 2 | 0.523 | no |
| merchant_structure | 4 | 0.501 | no |

Localized families: `audio_tower, video_tower, user_connectivity`.
Moved = this block's portrait differs early vs late. Not “this tower caused conversion”.
Towers sit next to graph families on the same clock. Contrast is moved / not-moved.

## 2. FSDS L1 — hops (LOGO)

RF-domain (all X) **0.846** · PO-risk **0.000164**

| hop | n | RF-mass | PO-mass | LOGO Δ | share |
|---|---:|---:|---:|---:|---:|
| tower_audio | 64 | 0.308 | 0.318 | +0.00003 | 0.464 |
| graph_merchant | 4 | 0.001 | 0.078 | +0.00002 | 0.325 |
| tower_video | 64 | 0.308 | 0.310 | +0.00001 | 0.211 |
| funnel_order | 5 | 0.075 | 0.084 | -0.00000 | 0.000 |
| funnel_user | 5 | 0.252 | 0.201 | -0.00023 | 0.000 |
| graph_item | 2 | 0.000 | 0.002 | -0.00001 | 0.000 |
| graph_user | 2 | 0.056 | 0.008 | -0.00004 | 0.000 |

LOGO Δ>0: dropping the hop **lowers** R (tied to the early/late Y gap).
Δ≤0: this hop is not a concept-gap source on this clock. Read RF-mass for P(X).

Drill continues inside localized hops: `tower_audio, tower_video, graph_user`.
Moved families: `audio_tower, video_tower, user_connectivity`.
L1 on those hops: `tower_audio` LOGO Δ=+0.00003; `tower_video` LOGO Δ=+0.00001; `graph_user` LOGO Δ=-0.00004.
Portrait movement ≠ Y-gap source. Funnel hops stay on the L1 board so graph/towers are compared, not isolated.
PO-risk here is ~1e-4 scale — log, not a story.

## 3. FSDS L2 — families inside localized hops

RF-domain **0.842** · PO-risk **0.000407**

| family | n | RF-mass | PO-mass | LOGO Δ | share |
|---|---:|---:|---:|---:|---:|
| audio_tower | 64 | 0.486 | 0.493 | -0.00006 | 0.000 |
| user_connectivity | 2 | 0.090 | 0.072 | -0.00007 | 0.000 |
| video_tower | 64 | 0.424 | 0.435 | -0.00001 | 0.000 |

## 4. FSDS L3 — LOCO on localized hop columns

If a 64-d tower is in the drill, LOCO is top-8 by L1 RF-mass, not all 64 fits.

| feat | hop | LOCO ΔR |
|---|---|---:|
| `aud_30` | tower_audio | +2.941128e-05 |
| `vid_28` | tower_video | +6.122892e-06 |
| `vid_08` | tower_video | -3.957873e-07 |
| `g_u_item_deg` | graph_user | -4.480021e-06 |
| `g_u_merch_deg` | graph_user | -1.028376e-05 |
| `vid_00` | tower_video | -1.546148e-05 |
| `aud_55` | tower_audio | -1.664398e-05 |
| `aud_43` | tower_audio | -1.794291e-05 |

LOCO ΔR>0: dropping the column lowers R. Still not a treatment effect.

## Read

- Localization AUCs: audio_tower 0.853, video_tower 0.829, user_connectivity 0.605, item_connectivity 0.523, merchant_structure 0.501.
- Moved: `audio_tower, video_tower, user_connectivity`. Quiet: `item_connectivity, merchant_structure`.
- User-pooled 64-d shop vectors can move even when the **current order's** 4 merchant-graph scalars do not: the tower grain is who-went-where in the left window, not `g_m_*` on this row.
- DGP has no W/Y injection. A moved tower is mix/coverage geometry after aggregation, not “video caused conversion”.
- L1–L3 answer **among those blocks, what is tied to the Y-gap** (if anything). Tiny LOGO on ~1e-4 PO-risk stays a log.
- Do not turn LOGO share into a unique importance ranking.

Graph localization 的粒（实体 / 关系 / 半径 / 边类型 / 社区 / 图对图）见 `docs/reports/Graph_Localization_Logic.md`。localization 无 Y；FSDS 只接 moved。

`PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds_drill.py`
