# DiffusionDB window-scheme sweep

Window scheme/size is part of the estimand: equal_count vs equal_time vs fixed width change n_by_T, calendar span, ΔȲ and selected tokens. Concatenating cfg/step/sampler into X surfaces hyperparam drift in VIMP.

| config | ΔȲ | HGB AUC | hp_vimp | n_by_T | span_h | top_blend |
|---|---:|---:|---:|---|---|---|
| `equal_count_k2` | 0.0331 | 0.602 | 0.000 | `{0: 1000, 1: 1000}` | `{0: 142.38333333333333, 1: 184.0}` | the, art, rutkowski, greg rutkowski, greg |
| `equal_count_k4` | 0.0681 | 0.616 | 0.000 | `{0: 500, 1: 500, 2: 500, 3: 500}` | `{0: 75.1, 1: 67.03333333333333, 2: 83.9, 3: 99.9}` | focus, artstation, the, sharp focus, render |
| `equal_time_k2` | 0.0389 | 0.624 | 0.000 | `{0: 1126, 1: 874}` | `{0: 163.11666666666667, 1: 163.21666666666667}` | greg, the, focus, greg rutkowski, digital art |
| `equal_time_k4` | 0.0613 | 0.605 | 0.000 | `{0: 527, 1: 599, 2: 490, 3: 384}` | `{0: 81.36666666666666, 1: 81.45, 2: 81.53333333333333, 3: 81.58333333333333}` | focus, artstation, painting of, the, render |
| `width_48h` | 0.0744 | 0.572 | 0.000 | `{0: 302, 1: 333, 2: 381, 3: 300, 4: 273, 5: 248, 6: 163}` | `{0: 47.85, 1: 47.63333333333333, 2: 47.9, 3: 47.81666666666667, 4: 46.483333333333334, 5: 47.6, 6: 38.483333333333334}` | 4k, render, portrait, artstation, high |
| `equal_count_k2_hp` | 0.0331 | 0.602 | 0.103 | `{0: 1000, 1: 1000}` | `{0: 142.38333333333333, 1: 184.0}` | hp_cfg, hp_step, greg, the, art |

## Read
- `equal_count_*`: balanced n, unequal calendar width
- `equal_time_*`: equal calendar span, unequal n
- `width_*`: fixed hours; early/late = first/last occupied bin
- `*_hp`: TF-IDF ⊕ cfg/step/sampler — check `hp_vimp_share` and top_cov

