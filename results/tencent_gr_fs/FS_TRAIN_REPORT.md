# TencentGR auto-1000 → feature-select → train

- users: **6000**
- generated dims: **1000** (target 1000)
- used after leak-drop: **1000**
- selected: **128** via `f`
- label: `future_cnv` (pos rate=0.157)

## Metrics

```json
{
  "hgb": {
    "auc": 0.8214498360104281,
    "ap": 0.48053643185183337,
    "acc": 0.8533333333333334,
    "sec": 0.24283838272094727,
    "n_selected": 128
  },
  "logreg": {
    "auc": 0.81399209486166,
    "ap": 0.48865699537733503,
    "acc": 0.8573333333333333,
    "sec": 2.5623631477355957,
    "n_selected": 128
  }
}
```

## Top-20 selected

```json
[
  {
    "name": "life_log1p_n_cnv",
    "score": 1090.419620914414
  },
  {
    "name": "log1p_pay_cnt",
    "score": 1090.419620914414
  },
  {
    "name": "log1p_abs_pay_cnt",
    "score": 1090.419620914414
  },
  {
    "name": "trans_cnv_to_exp",
    "score": 1086.6931851891777
  },
  {
    "name": "trans_exp_to_cnv",
    "score": 1074.140827566234
  },
  {
    "name": "30d_log1p_n_cnv",
    "score": 1054.2762717589294
  },
  {
    "name": "item_entropy_cnv",
    "score": 1052.248090203517
  },
  {
    "name": "sess_cnv_sess_rate",
    "score": 1001.0841463781658
  },
  {
    "name": "dec_hl7d_dec_cnv",
    "score": 1001.0097728556357
  },
  {
    "name": "sess_depth_cnv_mean",
    "score": 973.4289935569008
  },
  {
    "name": "life_cnv_share",
    "score": 955.1276348195488
  },
  {
    "name": "30d_n_cnv",
    "score": 942.3901172986189
  },
  {
    "name": "n_uniq_cnv_item",
    "score": 903.9045001562869
  },
  {
    "name": "life_n_cnv",
    "score": 902.891119820815
  },
  {
    "name": "pay_cnt",
    "score": 902.891119820815
  },
  {
    "name": "attr_cnv_wo_prior_clk_cnt",
    "score": 902.4805085175046
  },
  {
    "name": "ui_only_exp_share",
    "score": 894.205593113959
  },
  {
    "name": "ui_exp_per_item_mean",
    "score": 891.2403609264442
  },
  {
    "name": "log1p_abs_life_ctcvr",
    "score": 873.3539457009244
  },
  {
    "name": "x_pay_cnt__sess_bounce_rate",
    "score": 864.1658129790109
  }
]
```

## Recipe

1. **Generate**: combinatorial windows × funnel × decay × session × attribution × ARPU/deal × Markov × TOD × crosses ≈ 1000.
2. **Select**: VarianceThreshold → SelectKBest(F / MI) → top-k.
3. **Train**: HistGradientBoosting + LogisticRegression holdout AUC/AP.

Leakage note: for `pay_user`, raw `pay_cnt` / `life_n_cnv` / `arpu_sum` are dropped before selection.
