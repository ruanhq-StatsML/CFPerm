# Typed auxiliary losses: balance, correlation, effective rank

Not a scheduler board. Three objects next to TSS, on the same \((\hat c,\hat\delta)\) pair.

## What they are

- **Balance / inverse-CE / OGM-GE** is the usual multimodal regularizer: up-weight the lagging unimodal head. Symmetric. On a covariate hop that is the Polyak sign and puts clip-level text first.
- **Corr\((z^{(m)},W)\)** is Cohen's \(d\) in correlation units. It duplicates \(\hat c_m\). Do not sell it as a new intensity.
- **Typed \(L_{\mathrm{corr}}\)** is the redesign: \(\hat c_m(1-\hat c_{m'}/\sum_k\hat c_k)\,C_{mm'}^2\). Penalize a moving head that drags a quiet head. Large on the covariate hop; near zero on concept (because \(\hat c\) is quiet).
- **Effective rank** is a ranking, not a third shift type. \(\Delta\mathrm{erank}(\mathrm{Cov}\,X^{(m)})\) is location-invariant and does not rank the synthetic mean hop. Leave-one-out \(\Delta\mathrm{erank}(R)\) of the Amazon cosine Gram ranks Gift Cards and Subscription Boxes. A \(2\times 2\) hop Gram is a monotone rewrite of \(1-\cos\).

## Run

```
python3 tests/test_typed_aux_losses.py
python3 scripts/run_typed_aux_losses.py
```

Lead figure: `typed_aux_losses.png`. Table: `Typed_aux_losses.tex`. Method: `docs/method/Typed_aux_losses_note.tex`.
