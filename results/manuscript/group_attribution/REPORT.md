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

Post-hoc localization: pull subset indices (T groups, quartiles of the top X_j),
then pairwise **MMD** and **PO-risk**. Conditional means stay in the JSON; they are not the test.

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

