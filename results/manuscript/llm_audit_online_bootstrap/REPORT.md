# LLM audit — Palm–Nagler online bootstrap prototype

Two gates on the same `(Y, X, batch)` tables. HH chosen is not Y.

## Alignment

| | Frozen-ref online AR-bootstrap | Last-two `hop_fires` |
|---|---|---|
| Probe | Fit once on `D_ref` (batches `< n_ref`) | Refit on `B_{t-1}` every hop |
| Score `s_t` | MSE of frozen lstsq on `B_t` (concept-board probe) | 0-1 OOS error `e_now` |
| Null | \(\mu_{\mathrm{ref}}\) = mean lstsq MSE on `D_ref` mini-batches | previous hop `e_prev` |
| Statistic | \(\Delta_t = s_t - \mu_{\mathrm{ref}}\) | `e_now / e_prev` |
| UQ | Palm & Nagler AR-bootstrap CI | none (hard gate) |
| Fire | \(\mathrm{CI}_{\mathrm{lo}} > 0\) | ratio \(\ge \gamma\) and `e_prev ≥ e_floor` |
| Needs | ≥2 trail updates before a CI exists | first hop always quiet (`e_prev is None`) |

Defaults: `n_ref_batches=4`, labeled cut at batch `4`, `γ=1.5`, `n_boot=200`, `β=√2−1`, hop overlay `p=0.92`.

Frozen-ref answers: is current auditor error *significantly above* the reference window?
Last-two answers: did `P(Y|X)` hop between adjacent windows?

## Results

| Stream | Regime | pass | μ_ref | mean Δ | boot fires | first CI_lo>0 | delay | hop@cut | n hop fires |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| HH helpful | consistent | 0.49 | 0.171 | 0.037 | 1 | 8 | 4.000 | no | 0 |
| HH helpful | hop | 0.52 | 0.171 | 0.279 | 10 | 5 | 1.000 | yes | 2 |
| HH harmless | consistent | 0.49 | 0.170 | 0.045 | 9 | 5 | 1.000 | no | 1 |
| HH harmless | hop | 0.49 | 0.170 | 0.300 | 9 | 6 | 2.000 | yes | 2 |
| BeaverTails | native | 0.43 | 0.228 | 0.023 | 7 | 8 | — | — | 0 |
| BeaverTails | hop | 0.54 | 0.228 | 0.078 | 10 | 5 | 1.000 | no | 0 |
| WildGuard | native | 0.92 | 0.084 | -0.027 | 0 | — | — | — | 0 |
| WildGuard | hop | 0.23 | 0.084 | 0.747 | 9 | 5 | 1.000 | yes | 4 |
| ToxicChat | native | 0.78 | 0.208 | 0.012 | 4 | 9 | — | — | 4 |
| ToxicChat | hop | 0.24 | 0.208 | 0.119 | 10 | 5 | 1.000 | yes | 3 |

Hop regime: expect a *level* shift in Δ (CI sits above 0) and last-two fire at the cut. Consistent / native: Δ near 0. On 11 trail batches a 0.04 in-sample gap can still make `n_fires` look large — read mean Δ and the plot, not the raw fire count.

## What the labels did

- **HH hop vs consistent.** Frozen Δ jumps from ~0.04 to ~0.28–0.30. Last-two fires at the cut only on hop.
- **WildGuard** is the clean real-label pair: native Δ≈0 (0 bootstrap fires); Y-flip hop Δ≈0.75 and last-two fires at the cut.
- **BeaverTails** style-X is a weak map for `is_safe` (μ_ref≈0.23, near Bernoulli variance). Hop overlay moves Δ from 0.02 to 0.08 so the bootstrap can call excess, but last-two ratio at cut is 1.35 < γ=1.5.
- **ToxicChat** native is already non-quiet on last-two (human toxicity on these X is not a stationary map). Hop overlay still lifts frozen Δ (0.01 → 0.12).
- Bootstrap delay is at least 1 batch because a CI needs two trail updates. Last-two can fire on the first post-cut window.

## Cut-window last-two log

### HH helpful — `consistent`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.175 | 0.200 | 1.143 |
| 3 | no | 0.200 | 0.150 | 0.750 |
| 4 | no | 0.150 | 0.200 | 1.333 |
| 5 | no | 0.200 | 0.213 | 1.063 |
| 6 | no | 0.213 | 0.312 | 1.471 |

### HH helpful — `hop`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.175 | 0.200 | 1.143 |
| 3 | no | 0.200 | 0.150 | 0.750 |
| 4 | yes | 0.150 | 0.775 | 5.167 |
| 5 | no | 0.775 | 0.225 | 0.290 |
| 6 | yes | 0.225 | 0.425 | 1.889 |

### HH harmless — `consistent`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.200 | 0.150 | 0.750 |
| 3 | no | 0.150 | 0.225 | 1.500 |
| 4 | no | 0.225 | 0.238 | 1.056 |
| 5 | no | 0.238 | 0.200 | 0.842 |
| 6 | no | 0.200 | 0.250 | 1.250 |

### HH harmless — `hop`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.200 | 0.150 | 0.750 |
| 3 | no | 0.150 | 0.225 | 1.500 |
| 4 | yes | 0.225 | 0.738 | 3.278 |
| 5 | no | 0.738 | 0.250 | 0.339 |
| 6 | no | 0.250 | 0.325 | 1.300 |

### BeaverTails — `hop`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.512 | 0.412 | 0.805 |
| 3 | no | 0.412 | 0.463 | 1.121 |
| 4 | no | 0.463 | 0.625 | 1.351 |
| 5 | no | 0.625 | 0.525 | 0.840 |
| 6 | no | 0.525 | 0.625 | 1.190 |

### WildGuard — `hop`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.750 | 0.512 | 0.683 |
| 3 | no | 0.512 | 0.037 | 0.073 |
| 4 | yes | 0.037 | 0.912 | 24.333 |
| 5 | no | 0.912 | 0.113 | 0.123 |
| 6 | no | 0.113 | 0.062 | 0.556 |

### ToxicChat — `hop`

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.375 | 0.325 | 0.867 |
| 3 | no | 0.325 | 0.325 | 1.000 |
| 4 | yes | 0.325 | 0.488 | 1.500 |
| 5 | no | 0.488 | 0.137 | 0.282 |
| 6 | no | 0.137 | 0.050 | 0.364 |

## Serving

| Gate | Quiet | Fire |
|---|---|---|
| AR-bootstrap `CI_lo>0` | current error is compatible with `D_ref` | frozen auditor is significantly worse; do not treat current judge as gold |
| last-two hop | adjacent windows still the same map | policy-pack / judge swap between the last two batches; Top-k re-review |

```bash
PYTHONPATH=. python3 scripts/llm_audit_online_bootstrap_prototype.py
```

