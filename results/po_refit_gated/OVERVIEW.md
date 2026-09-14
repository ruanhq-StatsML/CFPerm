# Gated PO-risk reweighting — overview

**The method is online RFPerm + PO-risk.** The probe is the shallow
online RF on 上一批 (\(T=0\)); the score is instance `po_risk0`;
the weight is \(w=\sqrt{\texttt{po\_risk0}}\) on \(T=1\). Always-on
√PO IPTW and always-on DRE lose to last-two uniform on quiet and
real consecutive streams. Matching that fact: **periodic annealing**
— stay uniform, heat only when the consecutive OOS probe jumps,
then cool back.

This note is the overview. Numbered boards live in
`AGOD_po_refit_gated.md` (synth) and `AGOD_po_refit_real.md`
(batch-200 clocks). Method LaTeX (PO-risk reweighting, not AGOD):
`docs/po_risk/PO_risk_reweight_method.tex`. Tables:
`AGOD_performance_tables.tex`. Code:
`agod.online_rfperm.run_online_rfperm`.

## What we ship

Probe is the same shallow RF as the IPTW stream (`n_estimators=20`,
`max_depth=4`) on 上一批 as \(T=0\). Instance

\[
\texttt{po\_risk0}_i \;=\; |Y_i - \mu_0(X_i)|
\]

mixed with the batch gap (`instance_po_risk`, mix \(0.5\)). Consecutive
OOS, no K-fold (trees overfit an in-sample \(e_1>e_0\) and would fire
every hop):

\[
e_{\mathrm{now}}  = \mathrm{err}(\mu_0 \text{ fit } B_{t-1} \to B_t),\quad
e_{\mathrm{prev}} = \mathrm{err}(\mu_0 \text{ fit } B_{t-2} \to B_{t-1}).
\]

Fire iff \(e_{\mathrm{now}} / e_{\mathrm{prev}} \ge \gamma\) (default
\(1.5\)). Skip the first hop: there is no previous OOS.

- **Quiet** — last two batches, \(w=1\). This is the default.
- **Fire** — same rows. \(T=0\) stays \(1\); \(T=1\) gets
  \(w=\sqrt{\texttt{po\_risk0}}\) (mean \(1\)). Reweight, do not
  subset, do not drop the old batch. **Every observation** on the
  last two batches has its own `po_risk0_i` and `w_i`
  (`quantify_last_two`; hop history stores p10–p90).

Two levers, not one:

| lever | when the gate fires | what it does |
|---|---|---|
| **reweight** (`rfperm`) | last two, √PO on \(T=1\) | soft: upweight the new rows |
| **drop-old** (`resid`) | train on the new batch only | hard: throw away the stale map |

PO-tail subset localization (`q=0.3`) is off the board. At batch
size 200 that is ~60 rows and it loses to both levers.

The next-batch model is RF / XGBoost (not Ridge). Score is always
the **next** batch.

## Why PO-risk is significant when OOD is large

Squared loss splits OOS risk into noise plus the regression gap
\(\Delta_t^2=\mathbb{E}_t[(\mu_t^\star-\mu_0)^2]\). The gate ratio is
\(\approx 1+\Delta_t^2/\sigma^2\), so it fires only when the
\(P(Y\mid X)\) hop is large relative to noise. Instance
\(\sqrt{\texttt{po\_risk0}}\approx\sigma\) when \(\Delta_t\) is small
(mean-1 weights collapse to uniform) and
\(\approx|\mu_t^\star-\mu_0|(X)\) when \(\Delta_t\) is large (weights
track where the map flipped). The batch gap \(\delta\) tilts \(T=1\)
versus \(T=0\); \(r_i\) ranks within \(T=1\). That is the similar /
covariate wash versus the concept-cut gain. DRE cannot see \(\Delta_t\):
on concept it matches uniform (1.584 \(\approx\) 1.565) while
\(\sqrt{\mathrm{PO}}\) pays (1.466); on covariate it hurts
(0.728 \(\to\) 0.874). A global map flip makes \(\sqrt{\mathrm{PO}}\)
almost constant on \(T=1\); drop-old is then the harder lever.

## Annealing

The ratio is a one-shot heat. After a fire, \(e_{\mathrm{prev}}\)
on the next hop *is* that large OOS error, so

\[
\frac{e_{\mathrm{now}}}{e_{\mathrm{prev}}} \ll \gamma
\]

and the stream cools back to last-two uniform. That is the periodic
annealing schedule:

1. Quiet clock → \(w=1\).
2. Distribution hop → heat (√PO on the new batch).
3. Next hop → anneal (uniform again), unless another hop arrives.

A **one-shot** concept cut (synth \(B_4\)) heats once. A **periodic**
clock (occupancy day/night, PM2.5 episodes) reheats on a schedule.
Reweighting is built for (2) then (3). It is not built to ride a
cycle.

In-sample batch PO (\(e_1>e_0\) on the same trees) is not a gate:
it never anneals.

## Synth (expected)

RF, 4 seeds, 8 batches × 100. Quiet methods match to numerical
identity when fire\(=0\). DRE is always-on logistic \(p(x)\) on the
\emph{same} last-two rows.

| scene | uniform | DRE | reweight | drop-old | \(\Delta\) vs DRE |
|---|---:|---:|---:|---:|---:|
| similar | 0.540 | 0.541 | 0.540 | 0.540 | 0.001 |
| covariate | 0.728 | 0.874 | **0.720** | 0.852 | **0.154** |
| concept | 1.565 | 1.584 | 1.466 | **1.332** | **0.118** |
| mixed | 1.889 | 2.020 | 1.818 | **1.657** | **0.202** |

\(\Delta=\) DRE \(-\) reweight. Positive = \(\sqrt{\mathrm{PO}}\) wins.
DRE cannot see a \(P(Y\mid X)\) hop with stable \(P(X)\): on concept
it sits on uniform (1.584 \(\approx\) 1.565) while reweight pays.
On covariate it reweights the wrong axis and *hurts* (0.728 \(\to\) 0.874).

Concept path, cut at \(B_4\), RF, next-batch MSE:

| train \(t\) | test | uniform | DRE | reweight | drop-old |
|---|---|---:|---:|---:|---:|
| 3 | \(B_4\) ← cut | 5.054 | 5.051 | 5.054 | 5.054 |
| 4 | \(B_5\) | 2.149 | 2.254 | **1.556** | **0.755** |
| 5 | \(B_6\) | 0.503 | 0.506 | 0.503 | 0.503 |

The cut hop itself is unforecastable. The hop *after* the cut is
where annealing pays: reweight is the soft lever (2.15 → 1.56),
drop-old is the hard one (→ 0.76). DRE is *worse* than uniform on
that hop (2.25). Then fire goes dark and PO/uniform sit together.
XGB is the same shape.

Covariate hops should stay quiet. Drop-old overreacts (0.73 → 0.85).
DRE overreacts more (→ 0.87). Reweight barely moves (0.728 → 0.720).

## Real clocks, batch = 200

No shuffle. Row order is the clock. 24 consecutive batches. Local
Affec / Tencent / COCO packs were not on this machine.

RF RMSE (↓) / Acc (↑):

| dataset | clock | uniform | reweight | drop-old | fire |
|---|---|---:|---:|---:|---:|
| interstate | time | 714.8 | 714.8 | 698.1 | 0.00 |
| nyc_taxi | time | 2.586 | 2.585 | 2.584 | 0.05 |
| bike_hour | time | 57.06 | 56.93 | 57.83 | 0.05 |
| beijing_pm25 | time | **78.5** | 82.6 | 81.2 | 0.14 |
| electricity | time | 0.804 | 0.804 | 0.801 | 0.23 |
| airlines | time | 0.673 | 0.673 | 0.673 | 0.00 |
| occupancy | time | **0.859** | 0.831 | 0.859 | 0.27 |
| diabetes_readmit | shift | 0.638 | 0.638 | 0.638 | 0.00 |
| california | spatial | 0.693 | 0.698 | 0.697 | 0.05 |

Quiet clocks (airlines, interstate, readmit, and most taxi/bike hops)
match uniform. That is the intended anneal-to-default. Small
reweight moves on taxi / bike / electricity are noise-scale. XGB
says the same.

Drop-old helps interstate a little (715 → 698) and nowhere else
systematically.

## Periodic reheating (the remaining case)

Occupancy RFPerm fires at hops \(t=5,8,12,17,19,20\) (rate 0.27).
After each fire the next ratio collapses (e.g. \(t=5\) fire →
\(t=6\) ratio \(0.20\)) — that is annealing. Then the day/night
occupancy map flips again and the probe reheats. Several of those
ratios are \(10^7\)-scale because \(e_{\mathrm{prev}}\approx 0\)
on a constant-label stretch (empty rooms); \(\gamma=1.5\) is then
vacuous.

Beijing PM2.5 fires at \(t=2,12,22\) (rate 0.14): episode-scale,
not every batch, still enough to make √PO hurt the next hour-block
(78.5 → 82.6 RMSE).

Electricity is in between (rate 0.23): reweight is a wash (Acc
0.804 → 0.804 RF, 0.810 → 0.813 XGB).

So: **annealing works**. The leftover failure mode is a clock whose
\(P(Y\mid X)\) is itself periodic at the batch scale. Raising
\(\gamma\), requiring two consecutive jumps, or flooring
\(e_{\mathrm{prev}}\) would be the only follow-up with new
information. It is not a reason to go back to always-on IPTW.

## What we tried and dropped

- Always-on √PO / DRE — overreacts; uniform wins on quiet and on
  these real clocks.
- In-sample PO-ratio / in-sample \(e_1>e_0\) — fires every hop on
  trees; never anneals.
- K-fold residual on 上一批 — trees overfit train residuals; replaced
  by consecutive OOS.
- PO-tail subset (`q=0.3`) — too thin at batch 200; cut-hop MSE
  worse than full-batch reweight.
- Fire + drop the old batch **and** √PO (`resid_po`) — wash vs
  drop-old alone. The gain on a global map flip is the drop, not
  the weights.
- Ridge as the board predictor — swapped for RF / XGB.

## Default

Last-two uniform. Reweight only on a consecutive OOS jump, then
anneal. Drop-old is a separate, harder switch for a one-shot
concept cut. Do not subset. Do not IPTW every batch.
