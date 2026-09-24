# Sandbox distant bakeoff (away from AGOD methods)

Two unrelated classical ML jobs — **no** excess / PO / soft-burn / claim router.

## A. Next-step stream forecast

- packs ok: 3
- best_rmse counts: `{'hgb': 1, 'naive_last': 1, 'ridge': 1}`
- mean HGB RMSE lift vs naive: **0.18726140815649428**

| dataset | best | naive RMSE | ridge RMSE | HGB RMSE | HGB lift |
|---|---|---:|---:|---:|---:|
| `metro_interstate` | `hgb` | 854.6 | 717.4 | 486.9 | 0.4303 |
| `beijing_pm25` | `naive_last` | 24.06 | 25.2 | 24.68 | -0.02583 |
| `waymo_proxy` | `ridge` | 0.2642 | 0.1263 | 0.2226 | 0.1573 |

## B. DiffusionDB prompt themes

- n_prompts: 4000
- best k: **10** (silhouette=0.011175084822925056)
- top terms (first 3 clusters):
  - C0: art, concept, concept art, digital, artstation, detailed
  - C1: film, movie, life, cinematic, retro, poster
  - C2: detailed, highly, highly detailed, lighting, cinematic, artstation

## Effect reading (plain)

- Forecast: if HGB lift ≫ 0 on a pack, lag+X features beat last-value; if ≈0, series is near random-walk.
- Themes: silhouette picks a usable k; terms are descriptive clusters only.
- Explicitly **not** an AGOD efficiency / causal / transfer claim.

