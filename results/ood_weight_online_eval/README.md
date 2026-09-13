# OOD-risk as weights — online quantification

Compare **OOD / proximity / gradient-alignment weights** as sample weights
for closed-form Ridge on streaming boards.

Metrics: **online MSE**, **cum regret vs uniform**, **BWT** (batch-0 MSE),
**oracle gap** (vs best method on that seed).

> Waymo / Interstate / raw stock feeds are **not** in this repo.
> `interstate` and `stock` are **regime-shift surrogates**;
> `diabetes` is the real source→target table from `datasets/datasets.zip`.

## Methods

| id | weight rule |
| --- | --- |
| uniform | equal weights on all past (`ridge_past`) |
| hop | heatmap cosine IW |
| dga | Fan–Grangier–Ablin gradient alignment + EMA |
| mahal | Mahalanobis proximity to last-batch mean |
| attr | TSS-gated bank⊕hop adapter |

## Results

### `amazon`

| method | online MSE | regret vs uniform | BWT | oracle gap |
| --- | --- | --- | --- | --- |
| uniform | 1.4023 (0.0517) | +0.0000 (0.0000) | 0.6959 (0.0600) | 0.1976 (0.1535) |
| hop | 1.3784 (0.0429) | -0.1918 (0.1567) | 0.7711 (0.0605) | 0.0058 (0.0071) |
| dga | 1.4167 (0.0540) | +0.1147 (0.1482) | 1.0269 (0.0628) | 0.3123 (0.0967) |
| mahal | 1.3994 (0.0473) | -0.0231 (0.2753) | 0.8800 (0.0846) | 0.1746 (0.1282) |
| attr | 1.4469 (0.0576) | +0.3562 (0.1573) | 0.6959 (0.0600) | 0.5538 (0.1376) |

### `synthetic`

| method | online MSE | regret vs uniform | BWT | oracle gap |
| --- | --- | --- | --- | --- |
| uniform | 3.3666 (0.5138) | +0.0000 (0.0000) | 1.9154 (0.4453) | 9.8352 (2.5197) |
| hop | 3.0668 (0.5064) | -2.3981 (0.6418) | 2.0333 (0.4184) | 7.4371 (1.9224) |
| dga | 3.0365 (0.4742) | -2.6411 (0.5757) | 3.8973 (0.7638) | 7.1941 (2.1107) |
| mahal | 2.1372 (0.4849) | -9.8352 (2.5197) | 5.7020 (0.6418) | 0.0000 (0.0000) |
| attr | 2.7171 (0.4127) | -5.1958 (0.9483) | 1.9154 (0.4453) | 4.6394 (2.6873) |

### `diabetes`

| method | online MSE | regret vs uniform | BWT | oracle gap |
| --- | --- | --- | --- | --- |
| uniform | 1.9618 (0.0074) | +0.0000 (0.0000) | 1.9842 (0.0371) | 0.0266 (0.0187) |
| hop | 1.9596 (0.0098) | -0.0179 (0.0333) | 1.9876 (0.0394) | 0.0087 (0.0150) |
| dga | 1.9614 (0.0073) | -0.0033 (0.0091) | 1.9862 (0.0371) | 0.0233 (0.0129) |
| mahal | 1.9603 (0.0112) | -0.0118 (0.0446) | 1.9872 (0.0388) | 0.0148 (0.0264) |
| attr | 1.9795 (0.0041) | +0.1414 (0.0510) | 1.9842 (0.0371) | 0.1680 (0.0516) |

### `interstate`

| method | online MSE | regret vs uniform | BWT | oracle gap |
| --- | --- | --- | --- | --- |
| uniform | 2.2370 (0.0944) | +0.0000 (0.0000) | 2.0884 (0.0795) | 0.0000 (0.0000) |
| hop | 2.2867 (0.0963) | +0.3977 (0.0339) | 2.1189 (0.0814) | 0.3977 (0.0339) |
| dga | 2.2375 (0.0945) | +0.0040 (0.0020) | 2.0874 (0.0796) | 0.0040 (0.0020) |
| mahal | 2.2486 (0.0946) | +0.0927 (0.0151) | 2.0919 (0.0819) | 0.0927 (0.0151) |
| attr | 2.2989 (0.1249) | +0.4951 (0.2501) | 2.0884 (0.0795) | 0.4951 (0.2501) |

### `stock`

| method | online MSE | regret vs uniform | BWT | oracle gap |
| --- | --- | --- | --- | --- |
| uniform | 2.7016 (0.1417) | +0.0000 (0.0000) | 2.0580 (1.7298) | 0.2539 (0.1277) |
| hop | 2.8228 (0.1728) | +0.9690 (0.5566) | 2.1185 (1.6989) | 1.2230 (0.4884) |
| dga | 2.7123 (0.1479) | +0.0855 (0.1464) | 2.0580 (1.7298) | 0.3394 (0.1184) |
| mahal | 2.7633 (0.1935) | +0.4937 (0.4420) | 2.0581 (1.7297) | 0.7477 (0.4097) |
| attr | 2.6699 (0.1466) | -0.2539 (0.1277) | 2.0580 (1.7298) | 0.0000 (0.0000) |

## Run

```bash
python3 scripts/run_ood_weight_online_eval.py
python3 scripts/run_ood_weight_online_eval.py --quick
python3 scripts/run_ood_weight_online_eval.py --boards amazon,diabetes,interstate
```

Code: `Python/src/ood_weight_online_eval.py`
Skill: `.cursor/skills/ood-risk-weights-online/SKILL.md`
