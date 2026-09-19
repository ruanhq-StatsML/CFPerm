# FSDS iteration `iter04` — k-sweep × seed-sweep × MI

Ablation on combined recipe vs baseline FSDS / MI-combined.  
`k ∈ {8,10,12,15,18}`, `seed ∈ {0,1,2}` (45 runs). W2 never selects.

## Headline findings

1. **Seed dominates k.** At default k=15, W2 HGB swings **0.54 → 0.84** across seeds for both A and Z (GroupShuffleSplit + n_pos≈3). Any single-seed crowning is noise.
2. **Mean over seeds:** Z (cmean+π→FSDS) still edges A at k=15 (**0.720 vs 0.716**); gap tiny vs seed σ≈0.16.
3. **k sweet band ≈ 15–18** for mean W2; k≤10 is weak/unstable.
4. **MI combined hurts** at k=15 (mean W2 **0.592**, Δ−0.12 vs A). Only competitive at k=18. Stick to **F** in SelectKBest for this rare-pos pack (matches `auto_feats` F default).
5. Soft-corr / J*-only remain out of the combined chain (iter02–03).

## Mean W2 HGB by k (3 seeds)

| k | A baseline | Z combined | Z-MI |
|---:|---:|---:|---:|
| 8 | 0.574±0.25 | **0.616±0.29** | 0.602 |
| 10 | **0.559** | 0.557 | 0.544 |
| 12 | **0.648** | 0.637 | 0.593 |
| 15 | 0.716±0.15 | **0.720±0.16** | 0.592 |
| 18 | 0.714 | 0.714 | **0.720** |

## k=15 seed detail

| seed | A W2 | Z W2 | MI W2 |
|---:|---:|---:|---:|
| 0 | 0.768 | **0.774** | 0.540 |
| 1 | 0.543 | 0.543 | 0.416 |
| 2 | 0.836 | **0.843** | 0.819 |

## Method takeaway (讲武德)

Keep **Z = cmean guidance → π-stable F → official FSDS**.  
Next iters should **report mean±std over seeds** (and/or fix split with stratified rare-pos handling), not chase one seed. Prefer **k≈15**, F not MI.

Artifacts: `sweep_clean.csv`, `sweep_by_k.csv`.
