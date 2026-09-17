# Group attribution — permute T, not X

Serving gates ask whether P(Y|X) hopped in time. This is the other object:
which X_j carry a **group** difference. T is the group. CFPerm permutes T.
The 13 audit x_* are attribution coordinates for the group contrast, not for a time hop.

```bash
PYTHONPATH=. python3 scripts/prototype_group_attribution.py
```

Decision (same as `R/CFPerm_vimp.R`): fit group x X importance, permute T, B times,
feature p-value = fraction of nulls ≥ observed, across-feature threshold from the
0.95 quantile of those upper tails, reject iff ≥ 1 features clear it.

| Stream | T | groups | reject | hits | top |
|---|---|---:|---|---|---|
| HH helpful vs harmless | T=0 helpful queue, T=1 harmless queue | 2 | yes | x_n_toks, x_n_chars | x_n_toks, x_n_chars, x_hedge |
| BeaverTails vs ToxicChat | T=0 BeaverTails is_safe, T=1 ToxicChat human | 2 | no | — | x_i_count, x_newlines, x_refuse |
| sparse vs dense usable | T=0 sparse-only usable, T=1 dense-only usable | 2 | no | — | x_rrf_top1_mass, x_rank_corr, x_fuse_uniq |
| multi-step hop 0 / 1 / 2+ | T=0 first hop, T=1 second, T=2 later | 3 | yes | x_n_toks, x_n_chars | x_n_chars, x_n_toks, x_bang |

Post-hoc localization: pull subset indices, look at the **mean** (already computed),
then pairwise **MMD** and **PO-risk**. All three are already in the code.

## Conditional mean (already computed — look at this)

| Stream | subset | n | mean Y | mean top x |
|---|---|---:|---:|---:|
| HH helpful vs harmless | T0 | 1200 | 0.491 | 1.67 |
| HH helpful vs harmless | T1 | 1200 | 0.491 | 1.42 |
| HH helpful vs harmless | Q0 | 623 | 0.520 | 0.339 |
| HH helpful vs harmless | Q1 | 577 | 0.390 | 0.816 |
| HH helpful vs harmless | Q2 | 600 | 0.498 | 1.49 |
| HH helpful vs harmless | Q3 | 600 | 0.550 | 3.54 |
| BeaverTails vs ToxicChat | T0 | 1200 | 0.427 | 0.015 |
| BeaverTails vs ToxicChat | T1 | 1200 | 0.780 | 0.11 |
| BeaverTails vs ToxicChat | Q0 | 1853 | 0.584 | 0 |
| BeaverTails vs ToxicChat | Q1 | 547 | 0.667 | 0.273 |
| sparse vs dense usable | T0 | 1200 | 0.671 | 0.109 |
| sparse vs dense usable | T1 | 1200 | 0.358 | 0.109 |
| sparse vs dense usable | Q0 | 826 | 0.391 | 0.105 |
| sparse vs dense usable | Q1 | 564 | 0.566 | 0.106 |
| sparse vs dense usable | Q2 | 980 | 0.578 | 0.107 |
| sparse vs dense usable | Q3 | 30 | 0.900 | 0.347 |
| multi-step hop 0 / 1 / 2+ | T0 | 380 | 0.395 | 0.21 |
| multi-step hop 0 / 1 / 2+ | T1 | 379 | 0.417 | 0.232 |
| multi-step hop 0 / 1 / 2+ | T2 | 441 | 0.370 | 0.196 |
| multi-step hop 0 / 1 / 2+ | Q0 | 301 | 0.000 | 0.0418 |
| multi-step hop 0 / 1 / 2+ | Q1 | 303 | 0.000 | 0.105 |
| multi-step hop 0 / 1 / 2+ | Q2 | 297 | 0.582 | 0.204 |
| multi-step hop 0 / 1 / 2+ | Q3 | 299 | 0.997 | 0.5 |

Read the mean first. Pairwise MMD and PO-risk are the significance next to it.

## Pairwise subset MMD / PO-risk

| Stream | pair | n | MMD | MMD p | PO-risk | PO p | mean Y |
|---|---|---|---:|---:|---:|---:|---|
| HH helpful vs harmless | T0 vs T1 | 1200/1200 | 0.00492 | 0.0385 | 0.013 | 0.0385 | 0.491/0.491 |
| BeaverTails vs ToxicChat | T0 vs T1 | 1200/1200 | 0.279 | 0.0385 | 0.00546 | 0.0385 | 0.427/0.780 |
| sparse vs dense usable | T0 vs T1 | 1200/1200 | 0 | 1 | 0.00852 | 0.0385 | 0.671/0.358 |
| multi-step hop 0 / 1 / 2+ | T0 vs T1 | 380/379 | 0.0563 | 0.0385 | 2.57e-09 | 1 | 0.395/0.417 |
| multi-step hop 0 / 1 / 2+ | T0 vs T2 | 380/441 | 0.271 | 0.0385 | 0.000169 | 1 | 0.395/0.370 |
| multi-step hop 0 / 1 / 2+ | T1 vs T2 | 379/441 | 0.121 | 0.0385 | 0.000169 | 1 | 0.417/0.370 |

## Synthetic check

| DGP | reject | top |
|---|---|---|
| planted CATE on x1 in group 2 | yes | x1, x0, x3 |
| null: Y depends on X, not T | no | x1, x3, x5 |

Planted: three groups, only group 2 depends on `x1`. Null: Y depends on X, not on T.

Not this object: refusal rate, HH chosen, episode success, single-channel Recall.

