# PO × MSE board — electricity (NSW, ordered in time)

RF PO-risk should not collapse first; serving MSE is the series that breaks.
MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.
PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.
T=1 on the incoming batch. n_ref=10000, n_new=5000.
po_base=0.000159. mse_base=0.121. mmd_base=0.026323457816864404.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 7, n_watch = 0, n_x_shift = 0, n_tricky = 0, n_freeze = 0
- **keep training**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---|---|
| 0 | 3.06e-05 | 3.06e-05 | 0.119 | 0.119 | 0.0166 | 0.0166 | keep training | none — keep training every layer |
| 1 | 0.000205 | 0.000118 | 0.205 | 0.162 | 0.0226 | 0.0196 | keep training | none — keep training every layer |
| 2 | 5.48e-09 | 7.87e-05 | 0.232 | 0.185 | 0.121 | 0.0533 | keep training | none — keep training every layer |
| 3 | 3.54e-09 | 5.9e-05 | 0.159 | 0.179 | 0.124 | 0.0711 | keep training | none — keep training every layer |
| 4 | 1.86e-08 | 4.72e-05 | 0.188 | 0.181 | 0.155 | 0.0878 | keep training | none — keep training every layer |
| 5 | 1.98e-08 | 4.11e-05 | 0.27 | 0.211 | 0.0476 | 0.094 | keep training | none — keep training every layer |
| 6 | 3.33e-09 | 1.01e-08 | 0.122 | 0.194 | 0.00996 | 0.0914 | keep training | none — keep training every layer |

Read the two trends. Freeze only when both break.

