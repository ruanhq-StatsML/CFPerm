# Drill dissolve smoke — landability

**讲武德:** `W`=period. Leaf scores / dissolve sums are shift proxies — **not** ATE.

## Naming

PyPI [`dissolve`](https://pypi.org/project/dissolve/) = API-deprecation migrator. **Not used.**
Here **dissolve** = geopandas-style rollup: leaf score → `groupby(parent).sum/mean`.

## Setup

- leaves (W2 rows): 11963
- PO-risk: 1.14753e-06
- Top merchants=30, Top users=50

## Overlap vs top-down MMD (Jaccard@K)

### Merchant (dissolve shift_l2 vs MMD)

```
 k  jaccard  n_dissolve  n_ref
 5 0.666667           5      5
10 0.538462          10     10
20 0.600000          20     20
50 0.785714          50     50
```

### Merchant (dissolve τ̂² vs MMD)

```
 k  jaccard  n_dissolve  n_ref
 5 0.000000           5      5
10 0.052632          10     10
20 0.025641          20     20
50 0.052632          50     50
```

### User (global dissolve shift_l2 vs MMD-in-merchants)

```
 k  jaccard  n_dissolve  n_ref
 5      0.0           5      5
10      0.0          10     10
20      0.0          20     20
50      0.0          50     50
```

### User (cascade dissolve: top-M merchants → users vs MMD)

```
 k  jaccard  n_dissolve  n_ref
 5 0.000000           5      5
10 0.000000          10     10
20 0.025641          20     20
50 0.020408          50     50
```

### User (dissolve τ̂² vs MMD)

```
 k  jaccard  n_dissolve  n_ref
 5      0.0           5      5
10      0.0          10     10
20      0.0          20     20
50      0.0          50     50
```

## FSDS on dissolved user support

```
                   tag                             note  W2_hgb_auc    W2_ap  n_users  n_te   ok                                                 top5
     dissolve_l2_users eval_fallback_full_W2 n_users=50    0.601299 0.001639       50 11963 True e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec
dissolve_cascade_users eval_fallback_full_W2 n_users=50    0.601299 0.001639       50 11963 True e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec
   dissolve_tau2_users eval_fallback_full_W2 n_users=50    0.601299 0.001639       50 11963 True e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec
             mmd_users eval_fallback_full_W2 n_users=50    0.601299 0.001639       50 11963 True e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec
     random_head_users eval_fallback_full_W2 n_users=50    0.601299 0.001639       50 11963 True e_n_exp,u_n_events,u_n_exp,u_n_uniq_items,u_span_sec
```

## Landability verdict

1. **Implementable:** pure pandas groupby — no new heavy deps; plugs into existing grids.
2. **Merchant-level dissolve (shift_l2) tracks top-down MMD** reasonably (see Jaccard@K).
3. **Global user dissolve ≠ cascaded MMD users** (Jaccard≈0); use **cascade dissolve** (merchants→users) to match drill shape.
4. **τ̂² leaf dissolve** poorly matches MMD tops here — keep as optional PO side-channel, not the drill gate.
5. **Useful as:** cheap bottom-up prior / candidate gen before official cmean·MMD·FSDS; not a replacement.
6. **Do not** ship PyPI `dissolve` (API migrator); do not claim ATE.

Artifacts under `/workspace/results/tencent_gr_drill_dissolve`.
