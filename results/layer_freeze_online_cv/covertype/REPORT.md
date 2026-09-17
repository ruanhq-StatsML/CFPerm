# Layer-freeze online CV — covertype (geographic order, class 2 vs rest)

k layers → k+1 models (`model_0` … `model_k`). `model_i` trains the top i layers.
Incoming batch is T=1. Reference is T=0, n_ref=10000.
CV statistic = streaming PO-risk of that model's μ. i* = argmin_i PO-risk.

- hidden_dims = `[64, 32]`
- k = 3 → models 0…3
- n_ref = 10000, n_batches = 10
- **recommend train top 3 layer(s)** (median i*). Freeze below that.

| t | stream hop | model_0 | model_1 | model_2 | model_3 | i* |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.000381 | 0.000355 | 0.000356 | 0.000361 | 0.00034 | 3 |
| 1 | 0.0004 | 0.000698 | 0.000719 | 0.000646 | 0.000647 | 2 |
| 2 | 0.000297 | 0.000447 | 0.000483 | 0.000509 | 0.00044 | 3 |
| 3 | 0.000531 | 0.000569 | 0.000706 | 0.000851 | 0.000613 | 0 |
| 4 | 0.000414 | 0.000583 | 0.000788 | 0.0007 | 0.000382 | 3 |
| 5 | 0.000654 | 0.000967 | 0.00124 | 0.000762 | 0.000421 | 3 |
| 6 | 0.00082 | 0.00144 | 0.00168 | 0.000996 | 0.000546 | 3 |
| 7 | 0.00145 | 0.00226 | 0.00239 | 0.00162 | 0.000771 | 3 |
| 8 | 0.00127 | 0.00192 | 0.00214 | 0.0013 | 0.000917 | 3 |
| 9 | 0.00101 | 0.00135 | 0.00144 | 0.000887 | 0.000629 | 3 |

Read: if `model_0` (frozen) PO-risk stays high while a shallow `model_i` drops, unfreeze that far.
If a deep i spikes, you over-updated and washed the reference P(Y|X) — freeze those bottom layers.

