# Representation-first cold-start

repr_sep = attribute → margin separation on S → expand head when buffer≥N* / footprint stable. Instant learn = sep; learn = expand.

| mode | known | novel | forgetting | feat_hit |
| --- | --- | --- | --- | --- |
| repr_sep | 0.579 | 0.518 | 0.385 | 0.125 |
| weighted | 0.649 | 0.393 | 0.332 | 0.125 |
| uniform | 0.610 | 0.196 | 0.337 | — |

Δ novel (repr − uniform): **0.321**
Δ novel (repr − weighted): **0.125**

## repr_sep per-stage

| stage | novel | known | phase | gap_lift | expanded |
| --- | --- | --- | --- | --- | --- |
| 0 | nan | 0.488 | expand | 0.309 | True |
| 1 | 0.750 | 0.613 | expand | 0.346 | True |
| 2 | 0.286 | 0.636 | expand | 0.318 | True |

```bash
python3 scripts/run_repr_sep_cold_start.py --quick
```
