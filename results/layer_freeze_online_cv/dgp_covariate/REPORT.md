# PO × MSE board — DGP gradual covariate (bump; μ walks after batch 2)

RF PO-risk should not collapse first; serving MSE is the series that breaks.
MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.
PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.
OnlineRFPerm marks the shift-onset (frozen RF, T=MSE−E_ref, last-two hop).
T=1 on the incoming batch. n_ref=10000, n_new=2000.
po_base=4.025e-05. mse_base=0.1935. mmd_base=-0.0002545509573698146.
onset_hat=None, onset_rank=None, labeled onset=2.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 10, n_watch = 0, n_x_shift = 0, n_tricky = 0, n_freeze = 0, n_rfperm_hop = 0
- **keep training**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | RFPerm T | hop | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| 0 | 8.64e-06 | 8.64e-06 | 0.199 | 0.199 | -0.000318 | -0.000318 | 0.00924 |  | keep training | none — keep training every layer |
| 1 | 1.29e-05 | 1.08e-05 | 0.194 | 0.197 | 0.000667 | 0.000175 | 0.00349 |  | keep training | none — keep training every layer |
| 2 | 2.02e-05 | 1.39e-05 | 0.19 | 0.194 | 0.0054 | 0.00192 | 0.000426 |  | keep training | none — keep training every layer |
| 3 | 8.39e-06 | 1.25e-05 | 0.178 | 0.19 | 0.0221 | 0.00695 | -0.0129 |  | keep training | none — keep training every layer |
| 4 | 1.02e-05 | 1.21e-05 | 0.156 | 0.184 | 0.0677 | 0.0191 | -0.0288 |  | keep training | none — keep training every layer |
| 5 | 4.55e-06 | 1.12e-05 | 0.133 | 0.17 | 0.125 | 0.0442 | -0.0531 |  | keep training | none — keep training every layer |
| 6 | 9.44e-06 | 1.06e-05 | 0.0926 | 0.15 | 0.198 | 0.0837 | -0.0895 |  | keep training | none — keep training every layer |
| 7 | 4.03e-06 | 7.33e-06 | 0.0746 | 0.127 | 0.241 | 0.131 | -0.108 |  | keep training | none — keep training every layer |
| 8 | 1.76e-06 | 6e-06 | 0.0514 | 0.102 | 0.316 | 0.189 | -0.134 |  | keep training | none — keep training every layer |
| 9 | 1.51e-06 | 4.26e-06 | 0.0283 | 0.0759 | 0.397 | 0.255 | -0.156 |  | keep training | none — keep training every layer |

Read the two trends. Freeze only when both break.

