# FSDS-empowered RL closed loop

Three modules driven by FSDS typed drift `(c, δ)` — **locked signs**:
covariate large → shrink η; concept large → raise η / refresh reward heads.

## Modules

| id | role |
| --- | --- |
| meta_adapt | few-shot policy inner-loop conditioned on drift descriptor |
| reward_refresh | mix weights over reward components from typed drift |
| region_lr | TSS η per feature group |
| meta_iv | PO-risk uplift monitor + Important-Variable selection gates δ / η |

## Boards (proxies — full sims not installed here)

| board | stands in for | drift |
| --- | --- | --- |
| bandit_cov | Meta-World task params | covariate |
| bandit_concept | reward / map change | concept |
| bandit_mixed | both | both |
| maze_proxy | D4RL support mismatch | covariate |
| causal_proxy | CausalWorld mass/friction + objective flip | both |
| procgen_proxy | Procgen layout/style OOD + gem reward flip | both |
| meta_vision | Meta-Dataset domain hop + class map rotation | both |
| recsys_proxy | RecSys popularity cov + CTR concept | both |
| trading_proxy | finance vol/liquidity + risk↔return flip | both |
| dialogue_proxy | dialogue domain vocab + safety rubric | both |
| mm_hop_proxy | multimodal V/A/T hop (MSR-VTT-like) | both |

## Metrics (PO-risk empowered)

| metric | meaning |
| --- | --- |
| value | offline self-norm IW policy value ↑ |
| w_err / mix_mse | reward-mix fidelity vs true component weights ↓ |
| n_iv | # Important Variables from PO-risk VIMP |
| R_PO | meta-learner uplift monitor $E[\hat\tau^2]$ |
| ATE | $E[\mathrm{po}]$ batch-effect proxy |
| ESS | Kish ESS fraction of IW weights ↑ |
| regret | batch-mean reward − IW value ↓ |
| IV Jaccard | hop-to-hop stability of selected IVs ↑ |
| IV entropy | entropy of IV mass across feature groups |

## Results (offline IW value ↑ better)

### `bandit_cov`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 2.2124 (0.3031) | 0.2160 | 0.4073 | nan | nan | nan | 0.825 | -0.4051 | nan |
| fsds_full | 2.4500 (0.2677) | 0.0158 | 0.0089 | nan | nan | nan | 0.652 | -0.6427 | nan |
| fsds_meta_iv | 2.4395 (0.2537) | 0.0158 | 0.0089 | 7.0 | 0.0250 | 0.0064 | 0.655 | -0.6323 | 0.193 |
| naive_boost | 2.4760 (0.2739) | 0.2160 | 0.4073 | nan | nan | nan | 0.598 | -0.6687 | nan |

### `bandit_concept`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 2.6309 (0.1964) | 0.5753 | 3.1388 | nan | nan | nan | 0.837 | -0.3461 | nan |
| fsds_full | 3.0117 (0.2061) | 0.0410 | 0.0295 | nan | nan | nan | 0.620 | -0.7269 | nan |
| fsds_meta_iv | 2.9881 (0.2067) | 0.0410 | 0.0295 | 7.0 | 0.0393 | 0.0364 | 0.654 | -0.7033 | 0.266 |
| naive_boost | 3.0299 (0.2088) | 0.5753 | 3.1388 | nan | nan | nan | 0.620 | -0.7451 | nan |

### `bandit_mixed`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 2.7941 (0.4365) | 0.5753 | 3.4288 | nan | nan | nan | 0.830 | -0.3866 | nan |
| fsds_full | 3.1448 (0.4156) | 0.0421 | 0.0285 | nan | nan | nan | 0.590 | -0.7374 | nan |
| fsds_meta_iv | 3.1140 (0.4072) | 0.0421 | 0.0285 | 7.0 | 0.0396 | 0.0393 | 0.620 | -0.7065 | 0.228 |
| naive_boost | 3.1759 (0.4274) | 0.5753 | 3.4288 | nan | nan | nan | 0.556 | -0.7685 | nan |

### `maze_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 0.7857 (0.0517) | 0.4546 | 0.3027 | nan | nan | nan | 0.802 | -0.2630 | nan |
| fsds_full | 0.9558 (0.0528) | 0.0926 | 0.0043 | nan | nan | nan | 0.605 | -0.4332 | nan |
| fsds_meta_iv | 0.9347 (0.0381) | 0.0926 | 0.0043 | 6.0 | 0.0049 | 0.0065 | 0.653 | -0.4121 | 0.291 |
| naive_boost | 0.9676 (0.0524) | 0.4546 | 0.3027 | nan | nan | nan | 0.588 | -0.4450 | nan |

### `causal_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | -0.4332 (0.0233) | 0.3342 | 0.1177 | nan | nan | nan | 0.919 | -0.0259 | nan |
| fsds_full | -0.3602 (0.0263) | 0.3620 | 0.3085 | nan | nan | nan | 0.650 | -0.0988 | nan |
| fsds_meta_iv | -0.3395 (0.0258) | 0.3620 | 0.3085 | 6.0 | 0.0036 | 0.0126 | 0.565 | -0.1196 | 0.354 |
| naive_boost | -0.3595 (0.0252) | 0.3342 | 0.1177 | nan | nan | nan | 0.672 | -0.0995 | nan |

### `procgen_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | -0.0189 (0.0145) | 0.3894 | 0.0321 | nan | nan | nan | 0.941 | -0.0127 | nan |
| fsds_full | 0.0444 (0.0236) | 0.1605 | 0.0130 | nan | nan | nan | 0.796 | -0.0760 | nan |
| fsds_meta_iv | 0.0456 (0.0238) | 0.1605 | 0.0130 | 7.0 | 0.0011 | 0.0003 | 0.801 | -0.0772 | 0.279 |
| naive_boost | 0.0541 (0.0193) | 0.3894 | 0.0321 | nan | nan | nan | 0.776 | -0.0856 | nan |

### `meta_vision`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 3.9315 (0.7155) | 0.5435 | 8.8115 | nan | nan | nan | 0.718 | -0.4545 | nan |
| fsds_full | 4.0920 (0.8386) | 0.0663 | 0.1346 | nan | nan | nan | 0.474 | -0.6151 | nan |
| fsds_meta_iv | 4.1030 (0.8371) | 0.0663 | 0.1346 | 7.0 | 0.0488 | 0.0559 | 0.502 | -0.6261 | 0.233 |
| naive_boost | 4.1402 (0.7746) | 0.5435 | 8.8115 | nan | nan | nan | 0.502 | -0.6632 | nan |

### `recsys_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 2.3592 (0.6105) | 0.3450 | 0.8719 | nan | nan | nan | 0.686 | -0.5962 | nan |
| fsds_full | 2.5362 (0.6279) | 0.0551 | 0.0393 | nan | nan | nan | 0.490 | -0.7731 | nan |
| fsds_meta_iv | 2.5256 (0.6197) | 0.0551 | 0.0393 | 7.0 | 0.0246 | 0.0181 | 0.521 | -0.7625 | 0.211 |
| naive_boost | 2.5290 (0.6468) | 0.3450 | 0.8719 | nan | nan | nan | 0.483 | -0.7659 | nan |

### `trading_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | -0.2684 (0.0232) | 0.3342 | 0.0870 | nan | nan | nan | 0.942 | -0.0044 | nan |
| fsds_full | -0.2174 (0.0229) | 0.0905 | 0.0115 | nan | nan | nan | 0.803 | -0.0553 | nan |
| fsds_meta_iv | -0.2176 (0.0212) | 0.0905 | 0.0115 | 6.0 | 0.0022 | -0.0153 | 0.811 | -0.0551 | 0.240 |
| naive_boost | -0.2173 (0.0223) | 0.3342 | 0.0870 | nan | nan | nan | 0.805 | -0.0554 | nan |

### `dialogue_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 2.9624 (0.1044) | 0.2931 | 1.6713 | nan | nan | nan | 0.759 | -0.4912 | nan |
| fsds_full | 3.1149 (0.1183) | 0.0802 | 0.2124 | nan | nan | nan | 0.540 | -0.6437 | nan |
| fsds_meta_iv | 3.1096 (0.1255) | 0.0802 | 0.2124 | 7.0 | 0.0279 | 0.0274 | 0.558 | -0.6384 | 0.274 |
| naive_boost | 3.1662 (0.0867) | 0.2931 | 1.6713 | nan | nan | nan | 0.553 | -0.6950 | nan |

### `mm_hop_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO | ATE | ESS ↑ | regret ↓ | IV-J ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| static | 2.8679 (0.3493) | 0.2677 | 1.1505 | nan | nan | nan | 0.836 | -0.3341 | nan |
| fsds_full | 3.1800 (0.4503) | 0.0539 | 0.0766 | nan | nan | nan | 0.494 | -0.6462 | nan |
| fsds_meta_iv | 3.1809 (0.4475) | 0.0539 | 0.0766 | 8.0 | 0.0258 | 0.0446 | 0.534 | -0.6471 | 0.176 |
| naive_boost | 3.0934 (0.3295) | 0.2677 | 1.1505 | nan | nan | nan | 0.501 | -0.5595 | nan |

## Run

```bash
python3 scripts/run_fsds_rl_closed_loop.py
python3 scripts/run_fsds_rl_closed_loop.py --quick
python3 scripts/run_fsds_rl_closed_loop.py --compare-meta-iv
python3 scripts/run_fsds_rl_closed_loop.py --diverse --compare-meta-iv
python3 scripts/run_fsds_rl_closed_loop.py --boards all --compare-meta-iv
```

Code: `Python/src/fsds_rl_closed_loop.py`

## Opportunity notes

- Install `gymnasium` + Meta-World / Procgen / D4RL / CausalWorld to replace proxies.
- Diverse boards (RecSys / trading / dialogue / multimodal) stress PO-risk IV
  under preference flips, regime risk, safety rubrics, and V/A/T hops.
- `fsds_meta_iv` reuses the repo meta-learner (PO-risk uplift) + IV selection
  to focus TSS / descriptor on shift-driving features.
- Do **not** use `naive_boost` in production — it raises η on covariate shift.
