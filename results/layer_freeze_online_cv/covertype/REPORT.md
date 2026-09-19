# PO × MSE board — covertype (geographic order, class 2 vs rest)

RF PO-risk should not collapse first; serving MSE is the series that breaks.
MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.
PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.
T=1 on the incoming batch. n_ref=10000, n_new=10000.
po_base=0.0003324. mse_base=0.08222. mmd_base=0.047598851181268165.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 8, n_watch = 0, n_x_shift = 1, n_tricky = 1, n_freeze = 0
- **X shift (MMD)**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---|---|
| 0 | 0.000199 | 0.000199 | 0.182 | 0.182 | 0.0571 | 0.0571 | tricky | none — keep training every layer |
| 1 | 0.000251 | 0.000225 | 0.157 | 0.169 | 0.214 | 0.136 | X shift (MMD) | none — keep training every layer |
| 2 | 0.000319 | 0.000256 | 0.134 | 0.157 | 0.188 | 0.153 | keep training | none — keep training every layer |
| 3 | 0.000444 | 0.000303 | 0.077 | 0.137 | 0.183 | 0.161 | keep training | none — keep training every layer |
| 4 | 0.000389 | 0.00032 | 0.127 | 0.135 | 0.207 | 0.17 | keep training | none — keep training every layer |
| 5 | 0.000381 | 0.000357 | 0.219 | 0.143 | 0.213 | 0.201 | keep training | none — keep training every layer |
| 6 | 0.000499 | 0.000406 | 0.198 | 0.151 | 0.195 | 0.197 | keep training | none — keep training every layer |
| 7 | 0.000685 | 0.00048 | 0.14 | 0.152 | 0.177 | 0.195 | keep training | none — keep training every layer |

Read the two trends. Freeze only when both break.

