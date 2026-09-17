# MA of PO-risk — airlines (flight delay, time-ordered)

Raw PO-risk jitters. Causal moving average is the stability readout.
No online-bootstrap (repeated MLP inference cannot be afforded).
MA below 2× ref-split baseline → all-layer backprop.
n_ref=10000.

| n_new | batches | MA window | raw frac large | MA frac large | MA max / baseline | all-layer backprop |
|---:|---:|---:|---:|---:|---:|---|
| 20 | 200 | 50 | 0.00 | 0.00 | 0.42 | yes |
| 50 | 80 | 20 | 0.00 | 0.00 | 0.72 | yes |
| 100 | 40 | 10 | 0.00 | 0.00 | 0.95 | yes |
| 500 | 8 | 5 | 0.00 | 0.00 | 1.34 | yes |
