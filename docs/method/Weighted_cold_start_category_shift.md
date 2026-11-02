# Weighted cold-start under category shift (v2)

Cold-start = gate + boundary + auto(fused|sample) weights + known replay. Membership diagnoses; fused sample_w×subset_w adapts when concentrated.

## Protocol (four axes)

- novel_acc (sparse discovery)
- known retention / forgetting on pinned REF
- attribution: feat_hit + vimp_top_mass
- cold-start efficiency n* (budget → target novel_acc)
- boundary: incomparable rate vs label overlap

### Modes

| mode | meaning |
| --- | --- |
| weighted | auto fused/sample + private-live + known replay |
| membership | subspace membership only (diagnostic) |
| uniform | flat partial_fit |

## CIFAR-100 C-GCD proxy

| mode | known | novel | forgetting | feat_hit | vimp_top_mass |
| --- | --- | --- | --- | --- | --- |
| weighted | 0.613 | 0.393 | 0.328 | 0.125 | 0.422 |
| membership | 0.616 | 0.071 | 0.297 | 0.125 | 0.422 |
| uniform | 0.610 | 0.196 | 0.337 | — | — |

Δ novel (weighted − uniform): **0.196**

Δ novel (membership − uniform): **-0.125**

Δ forgetting (uniform − weighted; &gt;0 ⇒ weighted forgets less): **0.010**

## Cold-start efficiency (n* novel samples → novel_acc ≥ 0.60)

| mode | n* | curve (budget:acc) |
| --- | --- | --- |
| weighted | 16 | 4:0.25, 8:0.58, 12:0.50, 16:1.00, 24:0.83 |
| membership | 24 | 4:0.25, 8:0.33, 12:0.50, 16:0.58, 24:0.75 |
| uniform | 16 | 4:0.25, 8:0.33, 12:0.42, 16:0.67, 24:0.92 |

efficiency winner: **weighted**

## Office-Home PDA proxy

| overlap | shared acc | incomparable rate | adapt rate |
| --- | --- | --- | --- |
| 0.1 | 0.767 | 1.000 | 0.000 |
| 0.4 | 0.911 | 0.000 | 1.000 |
| 0.8 | 0.731 | 0.000 | 1.000 |

Directed adapt = detect gate + Jaccard boundary + `w=1+λ·membership` + `partial_fit`.

```bash
python3 scripts/run_weighted_cold_start.py --quick
```

## Multi-seed (2026,2027,2028)

- mean Δnovel (weighted−uniform): **0.131**
- mean Δforget (uniform−weighted): **0.015**
- efficiency winners: weighted, None, membership
