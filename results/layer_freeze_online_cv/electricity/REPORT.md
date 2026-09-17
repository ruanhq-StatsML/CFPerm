# Layer-freeze online CV — electricity (NSW, ordered in time)

k layers → k+1 models (`model_0` … `model_k`). `model_i` trains the top i layers.
Incoming batch is T=1. Reference is T=0, n_ref=10000.
CV statistic = streaming PO-risk of that model's μ. i* = argmin_i PO-risk.

- hidden_dims = `[64, 32]`
- k = 3 → models 0…3
- n_ref = 10000, n_batches = 12
- **recommend train top 2 layer(s)** (median i*). Freeze below that.

| t | stream hop | model_0 | model_1 | model_2 | model_3 | i* |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.0012 | 0.000191 | 0.000192 | 0.000182 | 0.000178 | 3 |
| 1 | 0.000877 | 0.000318 | 0.000273 | 0.000223 | 0.000208 | 3 |
| 2 | 0.000886 | 0.000944 | 0.00106 | 0.000967 | 0.000872 | 3 |
| 3 | 0.00149 | 0.00127 | 0.0011 | 0.00088 | 0.00084 | 3 |
| 4 | 0.00102 | 0.000877 | 0.000842 | 0.00116 | 0.00117 | 1 |
| 5 | 0.00364 | 0.00443 | 0.00487 | 0.00723 | 0.00746 | 0 |
| 6 | 0.0127 | 0.00689 | 0.00588 | 0.00484 | 0.00504 | 2 |
| 7 | 0.00416 | 0.00456 | 0.00446 | 0.00533 | 0.00563 | 1 |
| 8 | 0.000445 | 0.000724 | 0.000359 | 0.000339 | 0.000417 | 2 |
| 9 | 0.00214 | 0.00215 | 0.00243 | 0.00308 | 0.00306 | 0 |
| 10 | 0.00296 | 0.000228 | 0.000204 | 0.000568 | 0.000607 | 1 |
| 11 | 0.00124 | 1.34e-05 | 0.000129 | 3.16e-05 | 2.03e-05 | 0 |

Read: if `model_0` (frozen) PO-risk stays high while a shallow `model_i` drops, unfreeze that far.
If a deep i spikes, you over-updated and washed the reference P(Y|X) — freeze those bottom layers.

