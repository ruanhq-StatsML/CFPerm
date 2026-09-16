# Recsys: graph localization → two-step drill

These families are enough. No extra grains.
Use: **graph localization**, then a **two-step drill** only on the moved few.
PO-risk is **not** causal. RF-domain is a **portrait log**, not the score to optimize.
Video/audio are a messy DGP, not production embeddings.

## Protocol

1. Left-window clk/cnv → user—item, user—merchant, session co-click, shop projection.
2. Fake media DGP on the same graph (merchant attach → user mean-pool).
3. **Localize** (no Y) on five families:
   `user_connectivity` · `item_connectivity` · `merchant_structure` · `video_tower` · `audio_tower`.
   Moved = RF-domain AUC ≥ 0.55.
4. **Two-step drill** with Y=`y_post_clk_1d`, **only the moved families**:
   - Step 1: LOGO among those families
   - Step 2: LOCO columns (top-8 if a 64-d tower)

n=12564 p=136. Clock = median `t_end`. W=1 later half.
Shops attached=400 · users pooled=3808.
User-pool zero rate video early/late 0.035/0.047; mean L2-norm 2.79/3.07. Not a missingness artifact.

## 1. Graph localization (no Y)

| family | n | RF-domain AUC | moved |
|---|---:|---:|---|
| audio_tower | 64 | 0.853 | yes |
| video_tower | 64 | 0.829 | yes |
| user_connectivity | 2 | 0.605 | yes |
| item_connectivity | 2 | 0.523 | no |
| merchant_structure | 4 | 0.501 | no |

Localized families: `audio_tower, video_tower, user_connectivity`.
Moved = this block's portrait differs early vs late. Not “this tower caused conversion”.
Quiet families stop here. Drill does not open them.

## 2. Drill step 1 — LOGO on the moved few

RF-domain **0.842** · PO-risk **0.000407**

| family | n | RF-mass | PO-mass | LOGO Δ | share |
|---|---:|---:|---:|---:|---:|
| audio_tower | 64 | 0.486 | 0.493 | -0.00006 | 0.000 |
| user_connectivity | 2 | 0.090 | 0.072 | -0.00007 | 0.000 |
| video_tower | 64 | 0.424 | 0.435 | -0.00001 | 0.000 |

LOGO Δ>0: dropping the family **lowers** R (tied to the early/late Y gap).
Δ≤0: this family is not a concept-gap source on this clock. Read RF-mass for P(X).
Portrait movement ≠ Y-gap source. PO-risk ~1e-4 is a log, not a story.

## 3. Drill step 2 — LOCO on those columns

If a 64-d tower is in the moved set, LOCO is top-8 by step-1 RF-mass, not all 64 fits.

| feat | hop | LOCO ΔR |
|---|---|---:|
| `aud_33` | tower_audio | +2.409657e-06 |
| `vid_10` | tower_video | -2.886563e-05 |
| `aud_46` | tower_audio | -3.129216e-05 |
| `aud_38` | tower_audio | -4.347969e-05 |
| `g_u_item_deg` | graph_user | -5.623051e-05 |
| `g_u_merch_deg` | graph_user | -6.230930e-05 |
| `vid_24` | tower_video | -6.567858e-05 |
| `vid_16` | tower_video | -6.943418e-05 |

LOCO ΔR>0: dropping the column lowers R. Still not a treatment effect.

## Read

- Localization AUCs: audio_tower 0.853, video_tower 0.829, user_connectivity 0.605, item_connectivity 0.523, merchant_structure 0.501.
- Moved (drill these): `audio_tower, video_tower, user_connectivity`. Quiet (stop): `item_connectivity, merchant_structure`.
- Five families are enough. Do not add community / extra hops / funnel into this prototype.
- Use is two sentences: localize on the graph; two-step drill on the moved few.
- Do not turn LOGO share into a unique importance ranking.

`PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds_drill.py`
