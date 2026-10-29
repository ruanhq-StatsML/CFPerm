# After subset localization — strategy adjustment

Localization alone is not enough: map first-stage F → **tiers**, then
change how you aggregate and where you act.

| tier | rule | statistical action |
| --- | --- | --- |
| `NO_ACT` | F < 10 | local-retire IV; hold; do not quote IV ATE for decisions |
| `ESTIMATE_ONLY` | 10 ≤ F < 20 | publish local IV with WATCH; do not scale |
| `ACT` | F ≥ 20 | use local IV ATE; eligible for F- or precision-weights |

## Population shares

| tier | share |
| --- | --- |
| NO_ACT | 31.2% |
| ESTIMATE_ONLY | 15.6% |
| ACT | 53.1% |

| coverage | mean |
| --- | --- |
| identification (F≥10) | 68.8% |
| action (F≥20) | 53.1% |
| policy conflict (local NO_ACT vs global IV quote) | 31.2% |
| F heterogeneity ratio | 2169.63× |

## Aggregator |ATÊ − τ| (lower better)

| strategy | abs_err |
| --- | --- |
| global_iv | 0.020 (0.008) |
| naive_ols | 0.221 (0.050) |
| equal_pool_all | 0.181 (0.120) |
| equal_pool_id | 0.044 (0.055) |
| f_weighted | 0.012 (0.007) |
| precision_wt | 0.012 (0.007) |
| adaptive_act | 0.011 (0.007) |

## Insights

- best_aggregator=adaptive_act (|err|=0.011) vs global_iv=0.020 — weight by local F/precision, don't flat-pool NO_ACT
- tier_mix: NO_ACT=31% ESTIMATE_ONLY=16% ACT=53%
- coverage mean: identify=69% act=53%; policy_conflict=31% of mass would get a global IV quote while local F<10
- equal_pool_all |err|=0.181 vs f_weighted=0.012 — including weak segments pollutes the pooled ATE
- strength_gap: Q1_low F=3.4 (NO_ACT) vs Q4_high F=1246.8 (ACT)
- coverage: identify=75%  act=50%  hold_out=25%
- global_iv competitive here; still report tiers so NO_ACT mass is not oversold
- policy_conflict: 25% of population is NO_ACT locally but would get a global IV quote

ROI mapping (outside pure stats): ACT mass → scalable bets;
ESTIMATE_ONLY → learn more / small probes; NO_ACT → do not spend
on an IV-identified uplift claim.

## Integration with uplift monitoring

Cross-check meta-learner LOCO--AUUC / LIVE quintiles against these IV tiers
before REALLOCATE that claims identified uplift. Full actionable-items
iteration (tables + business gates B1--B5):

- `docs/method/Monitor_Uplift_FSDS_actionables.tex` (paste after Overall Recommendation)
- `docs/method/Monitor_Uplift_FSDS.tex` (full MVP writeup with section appended)
- `docs/method/Monitor_Uplift_FSDS_actionables.md`

## Run

```bash
python3 scripts/run_iv_selection_synth.py
```

Code: `Python/src/iv_selection_synth.py` (`subset_adaptive_strategy`, `strategy_suite`)
