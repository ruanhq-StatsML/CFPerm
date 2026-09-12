# Typed auxiliary losses: balance, correlation, effective rank

Not a scheduler board. Three objects next to TSS, on the same \((\hat c,\hat\delta)\) pair.

## How this is evaluated (two layers)

**Identification** (does the regularizer type the hop?). Oracle covariate-only vs concept-only. Inverse-CE must not match \(\hat c\) on a covariate hop. Typed \(L_{\mathrm{corr}}\) must fire on covariate and stay off on concept. \(\Delta\mathrm{erank}(\mathrm{Cov})\) must fail to rank a mean hop. Amazon leave-one-out \(\Delta\mathrm{erank}(R)\) must rank Gift / SubBox with the \(1-\cos\) hops. Table: `Typed_aux_losses.tex`. Figure: `typed_aux_losses.png`.

**Risk** (does that typing move the metric the hop is supposed to move?). Covariate-only: BWT / Brier on batch 0. Concept-only: post-change accuracy / CE. Amazon: online rating MSE and BWT MSE, plateau as the clock that already won that board. Do not pool into one multimodal score. Table: `Typed_aux_risk.tex`. Figure: `typed_aux_risk.png`.

Balance is allowed to hurt BWT. Typed \(L_{\mathrm{corr}}\) is an EWC-style leak penalty, not a joint-CE minimizer. Hop-Gram erank extra-shrink is stronger TSS, not a new intensity.

## Run

```
python3 -m pytest tests/test_typed_aux_losses.py
python3 scripts/run_typed_aux_losses.py
python3 scripts/run_typed_aux_losses.py --skip-risk
```

Method: `docs/method/Typed_aux_losses_note.tex`.
