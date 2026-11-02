# Weighted cold-start under category shift (v2)

Directed adapt = gate + boundary + upweight localized/private samples. C-GCD reports novel/forget/feat_hit + n* efficiency; PDA reports incomparable rate rising as overlap falls.

## Protocol (four axes)

- novel_acc (sparse discovery)
- known retention / forgetting on pinned REF
- attribution: feat_hit + vimp_top_mass
- cold-start efficiency n* (budget → target novel_acc)
- boundary: incomparable rate vs label overlap

### Modes

| mode | meaning |
| --- | --- |
| weighted | membership ∪ private-live labels |
| membership | subspace membership only |
| uniform | flat partial_fit |

## CIFAR-100 C-GCD proxy

| mode | known | novel | forgetting | feat_hit | vimp_top_mass |
| --- | --- | --- | --- | --- | --- |
| weighted | 0.638 | 0.268 | 0.367 | 0.125 | 0.422 |
| membership | 0.619 | 0.143 | 0.322 | 0.125 | 0.422 |
| uniform | 0.610 | 0.196 | 0.337 | — | — |

Δ novel (weighted − uniform): **0.071**

Δ novel (membership − uniform): **-0.054**

Δ forgetting (uniform − weighted; &gt;0 ⇒ weighted forgets less): **-0.029**

## Cold-start efficiency (n* novel samples → novel_acc ≥ 0.60)

| mode | n* | curve (budget:acc) |
| --- | --- | --- |
| weighted | 16 | 4:0.25, 8:0.33, 12:0.50, 16:0.92, 24:1.00 |
| membership | 16 | 4:0.25, 8:0.25, 12:0.25, 16:0.83, 24:0.67 |
| uniform | 16 | 4:0.25, 8:0.33, 12:0.42, 16:0.67, 24:0.92 |

efficiency winner: **weighted**

## Office-Home PDA proxy

| overlap | shared acc | incomparable rate | adapt rate |
| --- | --- | --- | --- |
| 0.1 | 0.767 | 1.000 | 0.000 |
| 0.4 | 0.911 | 0.000 | 1.000 |
| 0.8 | 0.774 | 0.000 | 1.000 |

Directed adapt = detect gate + Jaccard boundary + `w=1+λ·membership` + `partial_fit`.

```bash
python3 scripts/run_weighted_cold_start.py --quick
```
