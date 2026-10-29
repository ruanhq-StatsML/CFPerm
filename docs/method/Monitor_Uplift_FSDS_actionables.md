# Iterated actionable items: meta-learner localization × IV subset tiers

Append section for the MVP writeup *Minimal Viable Prototype: FSDS for monitoring the uplift model* (after “Overall Recommendation”). LaTeX: `Monitor_Uplift_FSDS_actionables.tex`.

## Why this iteration

1. **IV subset localization** — first-stage strength is not portable across strata; local \(F\) → action tiers before any spend claim.
2. **Meta-learner localization** (LOCO–AUUC / LIVE quintiles) — names where ranking broke; those names must be **cross-validated** against IV identification before a business move that assumes an identified incremental effect.

## IV subset tiers (synthetic extreme F-het, 8×n=2000)

| Tier | Local \(F\) | Share | Statistical action |
| --- | --- | --- | --- |
| `NO_ACT` | \(F<10\) | 31.2% | local-retire IV; hold; do not quote IV ATE |
| `ESTIMATE_ONLY` | \(10\le F<20\) | 15.6% | publish + WATCH; do not scale |
| `ACT` | \(F\ge 20\) | 53.1% | use local IV; F-/precision-weight eligible |

Coverage: identify 68.8% · act 53.1% · **policy conflict 31.2%** (local `NO_ACT` would still get a global IV quote).

| Aggregator | \|ATÊ − τ\| |
| --- | --- |
| naive_ols | 0.221 |
| equal_pool_all | 0.181 |
| global_iv | 0.020 |
| f_weighted / precision_wt | **0.012** |
| adaptive_act | **0.011** |

## Cross-validation matrix (meta × IV)

| Meta / AUUC on \(S\) | IV tier(\(S\)) | Joint read | Default action |
| --- | --- | --- | --- |
| AUUC healthy; \(\bar{\hat\tau}\approx\hat\tau^{\mathrm{obs}}\) | `ACT` | Rank OK + ID OK | Continue / scale; F-weight |
| AUUC healthy | `NO_ACT` | Sorts; not IV-identified | Monitor ranking; no IV uplift claim; RCT / margin only |
| AUUC red (cap) | `ACT` | Rank broken; ID OK | Reallocate on ranking; keep IV as relearn target |
| AUUC red | `NO_ACT` | Double red | Hold spend; RefreshRef / better \(Z\) or RCT |
| LOCO flags \(g\); mixed AUUC | `ESTIMATE_ONLY` | Soft ID | Watch / probe; owner sign-off |
| High \(\|\Delta\hat\tau^{\mathrm{obs}}\|\); PO accept | any | Descriptive ATE | Monitor only — never REALLOCATE on Layer-2 ATE alone |
| \(\mathrm{AUUC}_{ovlp}\gg\mathrm{AUUC}_{live}\) | ACT on overlap, NO_ACT on tail | Mix shift | NarrowSupport; do not global-IV the tail |

## Business gates (owners must enter)

- **B1** AUUC-green × IV-`NO_ACT` — scale vs RCT vs hold; do not brand \(\hat\tau^{\mathrm{obs}}\) as incremental profit.
- **B2** AUUC-red × IV-`ACT` — ranking cap vs exiting identified mass (P&L sign-off); prefer component refit before killing ACT mass.
- **B3** Instrument legitimacy — KEEP/WATCH/RETIRE is statistical; usable \(Z\) is business/legal; subset F-het ≥3× escalates.
- **B4** `ESTIMATE_ONLY` + P1 REALLOCATE — dual sign-off or downgrade to probe.
- **B5** Layer-2 ATE pair vs IV tier disagree — forbid Layer-2 shift if higher-ATE slice is `NO_ACT`.

## Actionable items list (v2)

| Pri | Code | Item |
| --- | --- | --- |
| P0 | GateSRM | Fail closed if SRM on \(T\) breaks |
| P0 | GateESS | Low ESS / domain AUC≈1 → RefreshRef only; defer LOCO cuts & IV weights |
| P1 | ReallocAUUC | Cap worst / continue best AUUC slice; attach IV tier on ticket |
| P1 | NarrowSupport | Deploy REF \(\hat\tau\) only on overlap when ovlp ≫ live |
| P1 | IV-NoActHold | `NO_ACT`: no IV ATE quote; no identified-uplift scale |
| P2 | HoldFeature | Keep LOCO top group \(g\) in prod; down-rank extreme-\(g\) segments |
| P2 | IV-ActWeight | On `ACT`, F-/precision-weight (not equal-pool) |
| P2 | IV-WatchProbe | `ESTIMATE_ONLY`: probe only; dual sign-off before P1 |
| P2 | CrossCheck | Every P1 REALLOCATE shows matrix cell + B-gate if conflict |
| P3 | Relearn | After trial failure under PO / ovlp SLA; prioritize IV-`ACT` mass for eval |
| P3 | RetireZ | Global/local RETIRE screens; WATCH if subset F-ratio ≥3× |
| P4 | MonitorATE | Layer-2 ATE pairs diagnostic unless PO reject **and** IV `ACT` |

**Justification.** Meta LOCO–AUUC localizes ranking failure; IV tiers localize identification failure; business gates resolve objective and instrument legitimacy. Either alone over-claims.

```bash
python3 scripts/run_iv_selection_synth.py
python3 scripts/run_fsds_rl_closed_loop.py --compare-meta-iv
```
