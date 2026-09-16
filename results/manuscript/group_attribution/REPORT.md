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

## Synthetic check

| DGP | reject | top |
|---|---|---|
| planted CATE on x1 in group 2 | yes | x1, x0, x3 |
| null: Y depends on X, not T | no | x1, x3, x5 |

Planted: three groups, only group 2 depends on `x1`. Null: Y depends on X, not on T.

Not this object: refusal rate, HH chosen, episode success, single-channel Recall.

