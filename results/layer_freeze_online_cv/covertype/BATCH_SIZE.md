# MA of PO-risk — covertype (geographic order, class 2 vs rest)

Raw PO-risk jitters. Causal moving average is the stability readout.
No online-bootstrap (repeated MLP inference cannot be afforded).
MA below 2× ref-split baseline → all-layer backprop.
n_ref=10000.

| n_new | batches | MA window | raw frac large | MA frac large | MA max / baseline | all-layer backprop |
|---:|---:|---:|---:|---:|---:|---|
| 500 | 60 | 5 | 0.00 | 0.00 | 1.53 | yes |
| 1000 | 30 | 5 | 0.20 | 0.07 | 2.26 | no |
| 2000 | 15 | 5 | 0.27 | 0.13 | 2.68 | no |
| 5000 | 6 | 5 | 0.17 | 0.00 | 1.42 | yes |
