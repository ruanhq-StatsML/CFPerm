# PO × MSE board — electricity (NSW, ordered in time)

PO+MSE both break → freeze that layer's training. PO only → watch.
MSE broken, PO quiet → not concept drift; read MMD of X (covariate shift).
Both quiet → keep training. No online-bootstrap.
T=1 on the incoming batch. n_ref=10000, n_new=5000.
po_base=0.0005145. mse_base=0.1217. mmd_base=0.026323457816864404.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 7, n_watch = 4, n_x_shift = 0, n_tricky = 0, n_freeze = 0
- **watch (concept?)**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---|---|
| 0 | 9.59e-05 | 9.59e-05 | 0.119 | 0.119 | 0.0191 | 0.0191 | keep training | none — keep training every layer |
| 1 | 0.000964 | 0.00053 | 0.193 | 0.156 | 0.0275 | 0.0233 | keep training | none — keep training every layer |
| 2 | 0.00166 | 0.000906 | 0.22 | 0.177 | 0.121 | 0.0559 | keep training | none — keep training every layer |
| 3 | 0.00597 | 0.00217 | 0.16 | 0.173 | 0.128 | 0.0739 | watch (concept?) | none — keep training every layer |
| 4 | 0.000905 | 0.00192 | 0.193 | 0.177 | 0.157 | 0.0906 | watch (concept?) | none — keep training every layer |
| 5 | 9.41e-06 | 0.0019 | 0.268 | 0.207 | 0.0575 | 0.0983 | watch (concept?) | none — keep training every layer |
| 6 | 6.08e-05 | 0.00172 | 0.123 | 0.193 | 0.00569 | 0.094 | watch (concept?) | none — keep training every layer |

Read the two trends. Freeze only when both break.

