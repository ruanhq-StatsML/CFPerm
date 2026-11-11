# Smart Bandit Router × FSDS/CFPerm — when is software good enough to go online?

**Audience:** evaluators / reviewers deciding whether the routing policy can take live traffic.  
**Software board (open-source + synthetic):** READY on `synthetic_credit`, `credit_g`, `adult`, `phoneme`; CONDITIONAL on `breast_cancer` (tiny exploration tax vs Keep).  
**Code:** `Python/src/smart_bandit_router.py` · **Skill:** `skills/smart_bandit_router_v1.yaml` · **Run:** `python3 scripts/run_smart_bandit_router.py`

---

## 一句话（中文）

子群颗粒度直接复用 **FSDS / CFPerm** 的子集定位（min_support + typed / CFPerm 坐标）；子群内部用 **Bayesian** 后验做 Thompson sampling；上线标准是 **READY / CONDITIONAL / HOLD**，不是 MSE。

---

## Architecture (three layers)

```
FSDS / CFPerm subset localization          ← granularity (who is in which bucket)
        ↓  leaf / named subset id
Bayesian Normal-mean Thompson Sampling     ← learning (per subset, not global LinTS)
        ↓  Keep | Route_to_Human | Rollback
Go-live gates → READY / CONDITIONAL / HOLD ← evaluator decision
```

We do **not** invent a separate bandit partitioner. Subset granularity is the same object FSDS already produces for monitoring and uplift localization.

---

## Layer 1 — Subset granularity = FSDS / CFPerm

| FSDS / CFPerm product | Role in the router |
| --- | --- |
| Model score / domain score | where the deployed scorer puts the unit |
| Domain / OOD intensity | covariate-loud axis (ĉ-style) |
| CFPerm residual / shift score | CFPerm attribution of shift along coords |
| Typed `(ĉ, δ̂)` proxies | covariate vs concept loudness (TSS-locked signs) |
| LOCO / LIVE quintiles / IV tiers *(production plug-in)* | **named** subsets from the existing FSDS monitor |
| `min_support` + MI gain | same support discipline as post-hoc localization / IV-subset tiers |

Software board helper (same axes, tree leaves under min_support):

```text
partition_fsds_cfperm_subsets(score, ood, cfperm_shift, y,
                              typed_c=…, typed_delta=…, min_support=…)
→ leaf_ids   # Bandit subgroup = FSDS/CFPerm subset
```

Production path: skip the tree and pass **named** FSDS subsets (LOCO top group, LIVE AUUC quintile, typed bin, CFPerm residual band, IV `ACT`/`NO_ACT` tier) as `group_ids`. Bayesian learning and go-live gates stay identical.

**Fairness floor** is defined on minority-mass FSDS subsets (protected Keep ≥ 80% of global Keep), not on “high-OOD ⇒ always Human,” so fairness does not fight correct Human routing on shift-loud mass.

---

## Layer 2 — Bayesian learning *inside* each subset

Each FSDS/CFPerm leaf keeps an independent `NormalMeanPosterior` per action (`Keep`, `Route_to_Human`, `Rollback`), warm-started from a hierarchical global prior.

- **Choose:** Thompson sample μ per action inside the leaf (cold leaf → global prior).
- **Update:** business PnL on the chosen action (evaluator currency).
- **Report:** per-leaf `P(best beats Keep)` and traffic-mass coverage — this is the soft gate `bayesian_subgroup_confidence`.

Global LinTS is intentionally **not** the go-live unit. The evaluator asks: *do mature FSDS subsets have confident posteriors?*

---

## Layer 3 — Evaluator verdict ladder

| Verdict | Meaning | Traffic |
| --- | --- | --- |
| **READY** | All hard + soft gates pass, incl. Bayesian subgroup mass ≥55% confident | open online |
| **CONDITIONAL** | Hard gates pass; soft gates still open | shadow / limited canary only |
| **HOLD** | Any hard gate fails | do not go online |

### Hard gates (must pass)

| Gate | Rule |
| --- | --- |
| `lift_vs_keep` | Bandit business reward ≥ Keep, or within 2% when Keep share ≥75% |
| `fairness_floor` | Protected Keep ≥ 80% of global Keep (traffic-weighted) |
| `granularity_support` | ≥2 FSDS/CFPerm leaves, min_support ≥ 20 |
| `human_share_sane` | Keep ≥ 20%, human ≤ 70% |

### Soft gates (needed for READY)

| Gate | Rule |
| --- | --- |
| `lift_nonnegative` | Strict lift ≥ 0 (2% exploration noise → CONDITIONAL, not HOLD) |
| `learning_stability` | Late-batch reward not collapsing vs early |
| `oracle_gap_bounded` | Regret vs oracle bounded, or clear lift vs Keep |
| `bayesian_subgroup_confidence` | ≥2 mature leaves and ≥55% mass where posterior confidently prefers Keep (P≥0.55) or a challenger that beats Keep (P≥0.60) |

---

## Software board (full suite)

| dataset | verdict | bandit | Keep | lift | bayes conf mass | fairness |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| synthetic_credit | **READY** | -19322 | -1136050 | +1116728 | 100% | Y |
| breast_cancer | **CONDITIONAL** | 13670 | 13880 | -210 | 100% | Y |
| credit_g | **READY** | 2998 | -2310 | +5308 | 100% | Y |
| adult | **READY** | 122555 | 86055 | +36500 | 100% | Y |
| phoneme | **READY** | 27795 | 6365 | +21430 | 96% | Y |

**Read for evaluators**

- Four of five boards are **READY** — safe to open online traffic under the fairness floor.
- `breast_cancer` is **CONDITIONAL** only because strict lift is slightly negative (~1.5% exploration tax on an already strong Keep baseline). Hard gates pass → shadow/canary, not HOLD.
- Bayesian subgroup confidence is high across the board (96–100% mass) — FSDS/CFPerm leaves are mature enough to trust local posteriors.

Reward board / learning path: `results/smart_bandit_router/smart_bandit_router.png`.

---

## How this plugs into the existing FSDS stack

| Existing FSDS artifact | Router use |
| --- | --- |
| `typed_shift_stepsize.py` `(ĉ, δ̂)` | typed axes for subset features / TSS signs |
| CFPerm / LOCO / MMD localization | named subset ids as `group_ids` |
| Monitor uplift LIVE quintiles / IV tiers | production granularity; `ACT` vs `NO_ACT` can gate Rollback / Human |
| `fsds_rl_closed_loop` reward refresh | optional shaped-reward mix (proxy / label-acq / survival) |
| min_support / ESS / SRM gates | same fail-closed support discipline before acting on a leaf |

Cross-links: `docs/method/Monitor_Uplift_FSDS_actionables.md`, `docs/method/IV_subset_strategy.md`, `docs/method/FSDS_RL_closed_loop.md`.

---

## Run / test

```bash
python3 -m pytest tests/test_smart_bandit_router.py -q
python3 scripts/run_smart_bandit_router.py          # full board
python3 scripts/run_smart_bandit_router.py --quick  # smoke
```

Auto suite table (regenerated): `results/smart_bandit_router/README.md`.  
This file is the **evaluator-facing** method note; it is not overwritten by the runner.
