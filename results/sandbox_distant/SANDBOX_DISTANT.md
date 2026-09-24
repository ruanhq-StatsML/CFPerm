# Sandbox distant bakeoff (away from AGOD methods)

Classical ML jobs — **no** excess / PO / soft-burn / claim router.

## A. Next-step stream forecast

- packs ok: 3
- best_rmse counts: `{'hgb': 1, 'naive_last': 1, 'ridge': 1}`
- mean HGB RMSE lift vs naive: **0.16083352137950715**

| dataset | best | naive RMSE | ridge RMSE | HGB RMSE | HGB lift |
|---|---|---:|---:|---:|---:|
| `metro_interstate` | `hgb` | 912.4 | 791.2 | 559.5 | 0.3868 |
| `beijing_pm25` | `naive_last` | 21.21 | 21.38 | 22.82 | -0.07566 |
| `waymo_proxy` | `ridge` | 0.2514 | 0.1251 | 0.2083 | 0.1713 |

## B. DiffusionDB prompt themes

- n_prompts: 4000
- best k: **10** (silhouette=0.011175084822925056)
- top terms (first 3 clusters):
  - C0: art, concept, concept art, digital, artstation, detailed
  - C1: film, movie, life, cinematic, retro, poster
  - C2: detailed, highly, highly detailed, lighting, cinematic, artstation

## C. Classical forecast → agent flywheel

**Headline:** 3 P0 flywheel opportunities; mean surprise_rate=0.003; 1 packs fell back to naive online

| id | pack | opportunity | priority |
|---|---|---|---|
| `FW-metro_interstate-hgb` | `metro_interstate` | specialize an HGB forecast agent | P0 |
| `FW-beijing_pm25-naive` | `beijing_pm25` | do NOT spend an HGB agent here | P1 |
| `FW-waymo_proxy-ridge` | `waymo_proxy` | cheap linear forecast agent | P0 |
| `FW-router` | `*` | pack→model router agent | P0 |

| pack | model | MAE | surprise_rate | retrains | fallback |
|---|---|---:|---:|---:|:---:|
| `metro_interstate` | `hgb` | 460.4 | 0 | 1 | Y |
| `beijing_pm25` | `naive_last` | 16.28 | 0.01 | 6 | N |
| `waymo_proxy` | `ridge` | 0.09833 | 0 | 6 | N |

## Effect reading (plain)

- Forecast: HGB lift ≫ 0 ⇒ learnable pack; ≈0 ⇒ random-walk (naive agent).
- Flywheel: opportunity is **pack→model routing + surprise/retrain/fallback**, not deeper nets.
- Themes: silhouette/terms only — context sticker for other agents.
- Not an AGOD efficiency / causal claim.

