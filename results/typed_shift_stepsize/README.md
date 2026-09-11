# Typed shift stepsize (TSS)

Statistical method scaffold, not a product board. Plug in another labelled stream via `stream_from_bundle` / `TypedStream(X, y, batch)`.

## Estimators

- Covariate intensity \(c_m\): \(|d|\) of the modality coordinate-mean under consecutive batches \(W\). Absolute, not a simplex share.
- Concept intensity \(\delta_m\): two-fold excess 0-1 risk of a unimodal ridge head after mean-aligning \(X\). A location shift in \(P(X)\) is removed; a change in \(P(Y\mid X)\) is not.
- Signed map: \(\eta_m=\eta_0(1+\beta\hat\delta_m)/(1+\lambda\hat c_m)\), quiet freeze if both channels are off. Covariate **lowers** the step; concept **raises** it. Do not invert \(\pi_m=c_m/\sum c\).

## What to compare

Scalar LR schedulers (constant, step, cosine, inv-time, plateau) cannot type the shift. \(\eta\propto\pi\) uses the same \(c_m\) with the wrong sign. Oracle TSS substitutes the known DGP channels.

Pre-specified metrics: **BWT** (accuracy on batch 0) under covariate-only; **post-change accuracy** under concept-only / both.

## Run

```
python3 scripts/run_typed_shift_stepsize.py
python3 scripts/run_typed_shift_stepsize.py --quick
```

MSR-VTT windows are an illustration of \(\hat c\) and the quiet text freeze (`tss_msrvtt_illustration.json`). A second labelled extract is required before claiming a downstream-task magnitude. Grafted video-id rotation (`tss_msrvtt_grafted_concept.json`) is a controlled concept check on real \(X\), not a naturalistic task.

Method writeup: `docs/method/MSRVTT_typed_shift_stepsize.tex`. Table: `TSS_vs_schedulers.tex`.
