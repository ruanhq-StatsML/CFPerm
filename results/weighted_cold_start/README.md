# Weighted cold-start under category shift (v2)

Cold-start v3 = detect/attribute → representation separation on S → dynamic head expand (buffer≥N* or footprint stable) → optional TTA on S. Weights still fire at expand; sep is the instant-learn step.

## Protocol (four axes)

- novel_acc (sparse discovery)
- known retention / forgetting on pinned REF
- attribution: feat_hit + vimp_top_mass
- cold-start efficiency n* (budget → target novel_acc)
- boundary: incomparable rate vs label overlap
- repr-sep: gap_lift on S before head expand

### Modes

| mode | meaning |
| --- | --- |
| weighted | auto fused/sample + private-live + known replay |
| membership | subspace membership only (diagnostic) |
| uniform | flat partial_fit |
| repr_sep | attribute → margin sep on S → expand head (+ optional TTA) |

## CIFAR-100 C-GCD proxy

| mode | known | novel | forgetting | feat_hit | vimp_top_mass |
| --- | --- | --- | --- | --- | --- |
| weighted | 0.649 | 0.393 | 0.332 | 0.125 | 0.422 |
| membership | 0.673 | 0.196 | 0.279 | 0.125 | 0.422 |
| uniform | 0.610 | 0.196 | 0.337 | — | — |

Δ novel (weighted − uniform): **0.196**

Δ novel (membership − uniform): **0.000**

Δ forgetting (uniform − weighted; &gt;0 ⇒ weighted forgets less): **0.006**

## Cold-start efficiency (n* novel samples → novel_acc ≥ 0.60)

| mode | n* | curve (budget:acc) |
| --- | --- | --- |
| weighted | 16 | 4:0.25, 8:0.33, 12:0.50, 16:0.83, 24:0.75 |
| membership | 16 | 4:0.25, 8:0.17, 12:0.50, 16:0.92, 24:0.67 |
| uniform | 16 | 4:0.25, 8:0.33, 12:0.42, 16:0.67, 24:0.92 |
| repr_sep | 16 | 4:0.33, 8:0.50, 12:0.58, 16:0.92, 24:1.00 |

efficiency winner: **membership**

## Office-Home PDA proxy

| overlap | shared acc | incomparable rate | adapt rate |
| --- | --- | --- | --- |
| 0.1 | 0.767 | 1.000 | 0.000 |
| 0.4 | 0.883 | 0.000 | 1.000 |
| 0.8 | 0.808 | 0.000 | 1.000 |

Pipeline v3: detect/attribute → repr-sep on S → expand head (weights) → optional TTA on S.
See `docs/method/Weighted_cold_start_logic.md` and `scripts/run_repr_sep_cold_start.py`.

```bash
python3 scripts/run_weighted_cold_start.py --quick
python3 scripts/run_repr_sep_cold_start.py --quick
```
