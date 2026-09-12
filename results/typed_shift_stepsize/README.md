# Typed shift stepsize (TSS)

Statistical method scaffold, not a product board. Plug in another labelled stream via `stream_from_bundle` / `TypedStream(X, y, batch)`.

## Estimators

- Covariate intensity \(c_m\): \(|d|\) of the modality coordinate-mean under consecutive batches \(W\). Absolute, not a simplex share.
- Concept intensity \(\delta_m\): two-fold excess 0-1 risk of a unimodal ridge head after mean-aligning \(X\). A location shift in \(P(X)\) is removed; a change in \(P(Y\mid X)\) is not.
- Signed map: \(s=\sqrt{n_{\mathrm{iter}}}\), \(\eta^\star=\eta_0(1+\beta\hat\delta_m/s)/(1+\lambda\hat c_m s)\). Video **缓升** toward the target; drops can be faster. Both quiet \(\Rightarrow\) hold last \(\eta_m\).

## What to compare

Scalar clocks (constant, cosine, plateau) copy one η onto every head. Per-head clocks: `plateau_m`, `restart_m` (SGDR reset on δ_m), `polyak_m` (loss-proportional). Typed maps (`tss`, `inv_c`, `fsds_pi`) set η for the **next** epoch.

Pre-specified metrics: **BWT** under covariate-only; **post-change accuracy** under concept-only / both. Look at \(\bar\eta_v,\bar\eta_a,\bar\eta_t\) to see whether a method is actually per-head.

On the oracle stream (6 seeds): global cosine/plateau keep \(\eta_v=\eta_a=\eta_t\). `plateau_m` shrinks the video head under covariate-only \((0.047,0.070,0.078)\) and stays competitive after concept drift. `polyak_m` raises video LR when unimodal CE rises, which is the wrong sign for pure covariate shift. TSS tracks `inv_c` here because the ridge \(\hat\delta\) jump is small (\(0.010\to 0.034\)); a stronger concept estimator is the next knob, not a different loop.

## Run

```
python3 scripts/run_typed_shift_stepsize.py
python3 scripts/run_typed_shift_stepsize.py --quick
```

MSR-VTT windows are an illustration of \(\hat c\) and the quiet text freeze (`tss_msrvtt_illustration.json`). A second labelled extract is required before claiming a downstream-task magnitude. Grafted video-id rotation (`tss_msrvtt_grafted_concept.json`) is a controlled concept check on real \(X\), not a naturalistic task.

Method writeup: `docs/method/MSRVTT_typed_shift_stepsize.tex`. Table: `TSS_vs_schedulers.tex`.
