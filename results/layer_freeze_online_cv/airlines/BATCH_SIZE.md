# Batch size is not when to update — airlines (flight delay, time-ordered)

n_ref=10000, no online-bootstrap, no repeated MLP inference.
n_new=20: PO-risk is a seismograph (cv≈1.2), never a stable hop. frac large = 0.

| n_new | batches | frac large | mean PO-risk | std |
|---:|---:|---:|---:|---:|
| 20 | 200 | 0.00 | 2.35e-07 | 2.84e-07 |
| 50 | 80 | 0.00 | 2.67e-07 | 3.46e-07 |
| 100 | 40 | 0.00 | 2.50e-07 | 3.46e-07 |
| 500 | 8 | 0.00 | 3.73e-07 | 6.20e-07 |
