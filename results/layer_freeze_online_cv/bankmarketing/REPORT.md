# PO × MSE board — bank-marketing (campaign order)

RF PO-risk should not collapse first; serving MSE is the series that breaks.
MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.
PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.
OnlineRFPerm marks the shift-onset (frozen RF, T=MSE−E_ref, last-two hop).
T=1 on the incoming batch. n_ref=10000, n_new=5000.
po_base=5.304e-07. mse_base=0.02286. mmd_base=0.08254973501278862.
onset_hat=4, onset_rank=None, labeled onset=None.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 6, n_watch = 0, n_x_shift = 2, n_tricky = 0, n_freeze = 0, n_rfperm_hop = 1
- **X shift (MMD)**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | RFPerm T | hop | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| 0 | 1.68e-07 | 1.68e-07 | 0.0325 | 0.0325 | 0.288 | 0.288 | 0.0117 |  | keep training | none — keep training every layer |
| 1 | 7.5e-09 | 8.8e-08 | 0.0383 | 0.0354 | 0.477 | 0.383 | 0.0167 |  | keep training | none — keep training every layer |
| 2 | 1.17e-06 | 4.5e-07 | 0.0429 | 0.0379 | 0.647 | 0.471 | 0.024 |  | keep training | none — keep training every layer |
| 3 | 2.13e-08 | 3.43e-07 | 0.0431 | 0.0392 | 0.47 | 0.471 | 0.0182 |  | keep training | none — keep training every layer |
| 4 | 4.35e-08 | 2.83e-07 | 0.146 | 0.0606 | 0.426 | 0.462 | 0.123 | yes | X shift (MMD) | none — keep training every layer |
| 5 | 1.63e-07 | 2.82e-07 | 0.0775 | 0.0696 | 0.0853 | 0.421 | 0.0598 |  | X shift (MMD) | none — keep training every layer |

Read the two trends. Freeze only when both break.

