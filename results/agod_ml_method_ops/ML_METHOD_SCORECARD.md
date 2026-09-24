# ML method ops scorecard (20-min RSI)

**Headline:** excess above null under current evidence — evidential only, not a driver ID

Composes with `transfer_null` excess / null — **no core rewrite**.

## Partial excess (confounder partialling)

| case | excess_raw | excess_partial | Δexcess | confounder R² |
|---|---:|---:|---:|---:|
| volume-driven score | 0.345 | -0.00161 | 0.346 | 0.995 |
| real signal ⊥ volume | 0.214 | 0.24 | -0.0261 | 0.225 |

## Excess learning curve (real-signal case)

- reading: excess softens / noisy with n — check rare labels / confounders
- excess@20% → full: 0.243 → 0.2
- gain: -0.0436

## Page–Hinkley on excess stream

- alarm: **True** (index=6)
- reading: PH statistic crossed threshold at pair index 6 (skill-drop *candidate*; confirm before acting)

## Effect (what improved)

1. False skill from volume is **quantified** (Δexcess), not hand-waved.
2. Sample hunger of excess estimate is visible (curve).
3. Skill-drop is detectable online (PH) — dual to data-drift gates.

