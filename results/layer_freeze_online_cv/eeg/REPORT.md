# PO × MSE board — eeg-eye-state (time-ordered EEG)

RF PO-risk should not collapse first; serving MSE is the series that breaks.
MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.
PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.
OnlineRFPerm marks the shift-onset (frozen RF, T=MSE−E_ref, last-two hop).
T=1 on the incoming batch. n_ref=6000, n_new=1500.
po_base=0.0008644. mse_base=0.1643. mmd_base=0.03844862250978898.
onset_hat=2, onset_rank=None, labeled onset=None.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 5, n_watch = 0, n_x_shift = 3, n_tricky = 0, n_freeze = 0, n_rfperm_hop = 1
- **X shift (MMD)**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | RFPerm T | hop | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| 0 | 0.000774 | 0.000774 | 0.275 | 0.275 | 0.0943 | 0.0943 | 0.206 |  | keep training | none — keep training every layer |
| 1 | 0.00045 | 0.000612 | 0.198 | 0.236 | 0.116 | 0.105 | 0.0731 |  | keep training | none — keep training every layer |
| 2 | 0.000513 | 0.000579 | 0.644 | 0.372 | 0.114 | 0.108 | 0.35 | yes | X shift (MMD) | none — keep training every layer |
| 3 | 0.000301 | 0.00051 | 0.375 | 0.373 | 0.177 | 0.125 | 0.179 |  | X shift (MMD) | none — keep training every layer |
| 4 | 0.000639 | 0.000535 | 0.321 | 0.362 | 0.0818 | 0.117 | 0.179 |  | X shift (MMD) | none — keep training every layer |

Read the two trends. Freeze only when both break.

