# MA of PO-risk — electricity (NSW, ordered in time)

Raw PO-risk jitters. Causal moving average is the stability readout.
No online-bootstrap (repeated MLP inference cannot be afforded).
MA below 2× ref-split baseline → all-layer backprop.
n_ref=10000.

| n_new | batches | MA window | raw frac large | MA frac large | MA max / baseline | all-layer backprop |
|---:|---:|---:|---:|---:|---:|---|
| 500 | 60 | 5 | 0.00 | 0.00 | 0.39 | yes |
| 1000 | 30 | 5 | 0.00 | 0.00 | 0.90 | yes |
| 2000 | 15 | 5 | 0.13 | 0.00 | 1.68 | yes |
| 5000 | 6 | 5 | 0.33 | 0.50 | 4.37 | no |
