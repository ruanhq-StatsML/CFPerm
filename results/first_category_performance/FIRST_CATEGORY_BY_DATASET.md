# First-category performance by dataset (multiseed)

Seeds [2026, 2027, 2028], n_first=8 per novel category (earliest stream indices), hold out 35% of that window.

| dataset | uniform | weighted | both_fused | ΔW | ΔF | P(F>U) |
| --- | --- | --- | --- | --- | --- | --- |
| har | 0.556±0.091 | 0.741±0.105 | 0.741±0.105 | +0.185 | +0.185 | 67% |
| mice | 0.424±0.023 | 0.269±0.055 | 0.269±0.055 | -0.155 | -0.155 | 0% |
| cnae9 | 0.356±0.075 | 0.356±0.140 | 0.356±0.140 | -0.000 | -0.000 | 33% |
| jvowel | 0.528±0.142 | 0.444±0.104 | 0.444±0.104 | -0.083 | -0.083 | 0% |
| wine | 0.667±0.272 | 0.889±0.157 | 0.889±0.157 | +0.222 | +0.222 | 33% |
| cancer | 1.000±0.000 | 0.778±0.314 | 0.778±0.314 | -0.222 | -0.222 | 0% |
| news20b | 0.056±0.079 | 0.056±0.079 | 0.056±0.079 | +0.000 | +0.000 | 33% |
| covtype | 0.111±0.157 | 0.056±0.079 | 0.056±0.079 | -0.056 | -0.056 | 0% |

mean Δ weighted **-0.014**, mean Δ fused **-0.014**

Core API: `first_indices_per_category` → `record_first_category_performance`.
