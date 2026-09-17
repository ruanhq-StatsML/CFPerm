# PO × MSE board — covertype (geographic order, class 2 vs rest)

PO+MSE both break → freeze that layer's training. PO only → watch.
MSE broken, PO quiet → not concept drift; read MMD of X (covariate shift).
Both quiet → keep training. No online-bootstrap.
T=1 on the incoming batch. n_ref=10000, n_new=10000.
po_base=0.0003472. mse_base=0.08163. mmd_base=0.047598851181268165.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 8, n_watch = 0, n_x_shift = 1, n_tricky = 1, n_freeze = 0
- **X shift (MMD)**

| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | action | freeze training |
|---:|---:|---:|---:|---:|---:|---:|---|---|
| 0 | 0.000255 | 0.000255 | 0.178 | 0.178 | 0.0659 | 0.0659 | tricky | none — keep training every layer |
| 1 | 0.000189 | 0.000222 | 0.158 | 0.168 | 0.2 | 0.133 | X shift (MMD) | none — keep training every layer |
| 2 | 0.000523 | 0.000322 | 0.124 | 0.153 | 0.193 | 0.153 | keep training | none — keep training every layer |
| 3 | 0.00102 | 0.000497 | 0.0861 | 0.136 | 0.184 | 0.161 | keep training | none — keep training every layer |
| 4 | 0.000267 | 0.000451 | 0.135 | 0.136 | 0.214 | 0.171 | keep training | none — keep training every layer |
| 5 | 0.000222 | 0.000445 | 0.21 | 0.142 | 0.208 | 0.2 | keep training | none — keep training every layer |
| 6 | 0.000344 | 0.000476 | 0.204 | 0.152 | 0.206 | 0.201 | keep training | none — keep training every layer |
| 7 | 0.000159 | 0.000403 | 0.146 | 0.156 | 0.176 | 0.198 | keep training | none — keep training every layer |

Read the two trends. Freeze only when both break.

