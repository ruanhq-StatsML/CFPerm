# smart_bandit_router_v1 — when is it good enough to go online?

Router actions: `Keep` / `Route_to_Human` / `Rollback`.
Learning unit: **Bayesian Normal-mean posterior inside each granularity leaf**
(Thompson sample per subgroup). Global LinTS is not the go-live unit.

## Evaluator verdict ladder

| verdict | meaning | traffic |
| --- | --- | --- |
| **READY** | All hard + soft gates pass, incl. Bayesian subgroup mass ≥55% confident | open online |
| **CONDITIONAL** | Hard gates pass; Bayesian leaf confidence / stability incomplete | shadow or limited canary only |
| **HOLD** | Any hard gate fails (lift / fairness / granularity / human share) | do not go online |

### Hard gates (must pass)

- `lift_vs_keep` — bandit cumulative **business** reward ≥ Keep
- `fairness_floor` — protected Keep ≥ 80% of global Keep
- `granularity_support` — leaves under min_support + MI
- `human_share_sane` — Keep ≥ 20%, human ≤ 70%

### Soft gates (needed for READY)

- `lift_nonnegative` — strict lift ≥ 0 (2% exploration noise → CONDITIONAL)
- `learning_stability` — late batches not collapsing vs early
- `oracle_gap_bounded` — regret vs oracle bounded, or clear lift
- `bayesian_subgroup_confidence` — ≥2 mature leaves and ≥55% traffic mass
  where the leaf posterior confidently prefers Keep (P≥0.55 stable) or a
  challenger that beats Keep (P≥0.60)

## Suite board

| dataset | verdict | bandit | Keep | lift | bayes conf mass | fairness |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| synthetic_credit | **READY** | -96242 | -1136050 | +1039808 | 99% | Y |
| breast_cancer | **CONDITIONAL** | 13670 | 13880 | -210 | 100% | Y |
| credit_g | **READY** | 2998 | -2310 | +5308 | 100% | Y |
| adult | **READY** | 122555 | 86055 | +36500 | 100% | Y |
| phoneme | **READY** | 27795 | 6365 | +21430 | 96% | Y |

## Run

```
python3 scripts/run_smart_bandit_router.py
python3 scripts/run_smart_bandit_router.py --quick
```

Skill: `skills/smart_bandit_router_v1.yaml`
Code: `Python/src/smart_bandit_router.py`

## synthetic_credit — evaluator note

- verdict: **READY** — All hard + soft gates pass, including Bayesian subgroup coverage. Safe to open online traffic with the existing fairness floor.

- `lift_vs_keep` [hard]: **PASS** — Bandit business reward ≥ Keep, or within 2% when Keep share ≥75%
- `lift_nonnegative` [soft]: **PASS** — Strict lift ≥ 0 required for READY (else CONDITIONAL if hard gates pass)
- `fairness_floor` [hard]: **PASS** — Protected Keep rate ≥ 80% of global Keep
- `granularity_support` [hard]: **PASS** — Context leaves grown under min_support + MI
- `human_share_sane` [hard]: **PASS** — Keep ≥ 20%, human ≤ 70%; human≥1% unless nearly-full automation
- `learning_stability` [soft]: **PASS** — Late-batch reward not collapsing vs early (sign-safe)
- `oracle_gap_bounded` [soft]: **PASS** — Regret vs oracle < 70% of |oracle|, or clear lift vs Keep
- `bayesian_subgroup_confidence` [soft]: **PASS** — ≥2 mature leaves and ≥55% traffic mass in leaves where the Normal-mean posterior confidently prefers Keep or a challenger (P≥0.55 Keep-stable / P≥0.60 beats Keep)
