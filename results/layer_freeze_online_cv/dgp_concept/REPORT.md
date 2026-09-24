# PO × MSE board — DGP gradual concept (β flips after batch 2)

RF PO-risk should not collapse first; serving MSE is the series that breaks.
MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.
PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.
OnlineRFPerm marks the shift-onset (frozen RF, T=MSE−E_ref, last-two hop).
T=1 on the incoming batch. n_ref=10000, n_new=2000.
po_base=4.29e-05. mse_base=0.1321. mmd_base=-0.0002545509573698146.
onset_hat=None, onset_rank=None, labeled onset=2.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 10, n_watch = 5, n_x_shift = 0, n_tricky = 0, n_freeze = 0, n_rfperm_hop = 0
- **watch (concept?)**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | RFPerm T | hop | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| 0 | 1.92e-05 | 1.92e-05 | 0.135 | 0.135 | -0.000318 | -0.000318 | 0.00381 |  | keep training | none — keep training every layer |
| 1 | 1.11e-05 | 1.51e-05 | 0.136 | 0.136 | 0.000667 | 0.000175 | 0.00985 |  | keep training | none — keep training every layer |
| 2 | 1.65e-05 | 1.56e-05 | 0.14 | 0.137 | 0.000328 | 0.000226 | 0.0122 |  | keep training | none — keep training every layer |
| 3 | 4.04e-05 | 2.18e-05 | 0.169 | 0.145 | -0.000625 | 1.29e-05 | 0.0325 |  | keep training | none — keep training every layer |
| 4 | 0.000155 | 4.85e-05 | 0.198 | 0.156 | 0.000576 | 0.000125 | 0.0554 |  | keep training | none — keep training every layer |
| 5 | 0.000551 | 0.000155 | 0.236 | 0.176 | 1.52e-05 | 0.000192 | 0.0966 |  | watch (concept?) | none — keep training every layer |
| 6 | 0.00114 | 0.00038 | 0.256 | 0.2 | 0.000453 | 0.000149 | 0.136 |  | watch (concept?) | none — keep training every layer |
| 7 | 0.00265 | 0.000905 | 0.277 | 0.227 | 0.00195 | 0.000473 | 0.198 |  | watch (concept?) | none — keep training every layer |
| 8 | 0.00385 | 0.00167 | 0.263 | 0.246 | -0.000778 | 0.000443 | 0.232 |  | watch (concept?) | none — keep training every layer |
| 9 | 0.0049 | 0.00262 | 0.248 | 0.256 | -0.000794 | 0.000169 | 0.262 |  | watch (concept?) | none — keep training every layer |

Read the two trends. Freeze only when both break.

