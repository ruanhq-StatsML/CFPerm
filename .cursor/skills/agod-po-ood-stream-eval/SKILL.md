---
name: agod-po-ood-stream-eval
description: Evaluate AGOD PO-risk as streaming OOD score vs density-ratio (DRE) on real stream packs. Use when choosing datasets, sizing streams (n/d/batch), defining MSE/regret metrics, characterizing gradual-shift properties, or running/visualizing metro/stocks/beijing/waymo/affec PO-IPTW experiments.
---

# AGOD PO-OOD Stream Evaluation

Use this skill to **pick / size / characterize** streaming datasets and to run the green **√PO vs DRE** protocol.

## Core claim

On gradual-shift streams, treat PO-risk as a **label-aware OOD score**:

```python
w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))
model.fit(X, y, sample_weight=w / w.mean())
```

Compare against logistic **DRE** (`w ∝ p(cur|x)/p(ref|x)` on X only). Default soft mode is `sqrt` (not hard `prop`).

## When to run

- User mentions Waymo / Metro Interstate / stocks / Beijing PM2.5 / Affec / stream packs
- Asks for MSE, regret, OOD score, IPTW weights, or PO vs density-ratio
- Asks what dataset **size / properties** matter for the claim

## Dataset properties to record (required card)

For every pack, compute and log:

| property | definition | why it matters |
|---|---|---|
| `n` | usable ordered rows | need ≥ `batch_size * n_batches` |
| `d` | feature dim (after optional PCA) | RF/DRE stability; PCA→32 is fine for smoke |
| `batch_size` | default **100** | matches streaming AGOD protocol |
| `n_batches` | default **30–40** | enough for cum-regret curves |
| `freq` | hourly / daily / 0.1s / block | sets natural drift timescale |
| `target` | continuous y name | MSE/MAE regime (not Acc) |
| `proxy` | true if synthetic stand-in | disclose Waymo kinematics proxy |
| `y_std` | std of target | scale for relative metrics |
| `x_mean_drift_avg` | mean ‖μ_t − μ_{t−1}‖ across batches | covariate drift magnitude |
| `y_mean_drift_avg` | mean |ȳ_t − ȳ_{t−1}| | label/concept drift proxy |
| `gradual_shift_score` | `x_mean_drift_avg / std(X)` | primary X-axis for “when √PO wins” plots |

**Sizing rule of thumb**

- Minimum smoke: `n ≥ 2000`, `batch_size=100`, `n_batches=20`
- Preferred: `n ≥ 4000`, `n_batches=40`
- If class-sorted or non-temporal order → **fix order** (PC1 / timestamp) before streaming
- Prefer **real** packs; proxy only when dump missing (label `proxy: true`)

## Packed loaders in this repo

| key | path | notes |
|---|---|---|
| `metro_interstate` | `data/stream_packs/metro_interstate/` | traffic volume, hourly |
| `beijing_pm25` | `data/stream_packs/beijing_pm25/` | PM2.5, hourly |
| `stocks_SPY` / `QQQ` / `AAPL` | `data/stream_packs/stocks/` | next-day return |
| `waymo_proxy` | `data/stream_packs/waymo_proxy/` | kinematics proxy if real Waymo absent |
| `affec` | `results/affec_fsds/*xyw*cache.npz` | multimodal affect stream |

Load via `agod.stream_packs.LOADERS[name](root, max_n=...)`.

## Metrics protocol (do not invent new weight families)

Modes: `uniform`, `prop`, `sqrt`, `inv`, `dre`.

For t = 1…T−1:

1. Probe from batch t−1 predicts batch t → instance PO (residual, optional batch-PO mix)
2. Fit batch t with weights → evaluate on batch **t+1**
3. Log next-batch **MSE**, **MAE**
4. Regret:
   - `cum_regret_vs_uniform = Σ (MSE_mode − MSE_uniform)`
   - `cum_regret_vs_dre = Σ (MSE_mode − MSE_dre)` (negative ⇒ mode beats DRE)
   - `win_rate_vs_dre = mean(MSE_mode < MSE_dre)`

Always report **√PO vs DRE head-to-head win count**, plus **relative MSE** `MSE/uniform` so metro-scale does not dominate plots.

## Commands

```bash
PYTHONPATH=. python3 scripts/run_agod_po_ood_stream_metrics.py \
  --batch-size 100 --n-batches 40 \
  --datasets metro_interstate beijing_pm25 stocks_SPY stocks_QQQ waymo_proxy affec \
  --out results/agod_po_ood_metrics
```

Outputs under `results/agod_po_ood_metrics/`:

- `summary.json` — trajectories + dataset property cards
- `PO_OOD_METRICS_REPORT.md`
- plots: `mse_rel_uniform_bars.png`, `cum_regret_sqrt_vs_dre(_norm).png`, `shift_vs_po_advantage.png`, per-dataset traj

## Visualization checklist

1. **Relative MSE bars** (`MSE/uniform`) — cross-dataset comparable
2. **Cum regret √PO vs DRE** — preferably also a **normalized** version per dataset
3. **Shift score vs √PO advantage** — property plot for the paper/story
4. Per-dataset MSE trajectories

## Interpretation guide

- √PO should beat DRE on gradual / label-risk drift (Affec, traffic, PM2.5, Waymo proxy)
- Uniform can still win absolute MSE; the claim is **vs DRE**, not “beats everything”
- Hard `prop` often overfits the current batch — keep `sqrt` as default
- `inv` can win on noisy returns (stocks) — report honestly
- DRE ignores Y → can explode when concept drifts with mild covariate shift

## Anti-patterns

- Do not add density-ratio / DGA variants as the main method
- Do not shuffle away temporal order
- Do not treat high-cardinality IDs as regression targets without coarsening
- Do not claim real Waymo if using `waymo_proxy`
- Do not expand the weight taxonomy beyond uniform/prop/sqrt/inv/dre without an explicit ask
