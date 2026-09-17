# PO × MSE board — covertype (geographic order, class 2 vs rest)

PO broken and MSE broken → freeze. PO broken, MSE holds → watch. Both quiet → keep training.
No online-bootstrap. Freeze-depth PO_Dict / MSE_Dict only on freeze hops.
T=1 on the incoming batch. n_ref=10000, n_new=10000.
po_base=0.0003472. mse_base=0.0819.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 8, n_watch = 0, n_freeze = 0
- **keep training**

| t | PO-risk | PO MA | MSE | MSE MA | action |
|---:|---:|---:|---:|---:|---|
| 0 | 0.000255 | 0.000255 | 0.181 | 0.181 | keep training |
| 1 | 0.000189 | 0.000222 | 0.158 | 0.169 | keep training |
| 2 | 0.000523 | 0.000322 | 0.124 | 0.154 | keep training |
| 3 | 0.00102 | 0.000497 | 0.0819 | 0.136 | keep training |
| 4 | 0.000267 | 0.000451 | 0.128 | 0.134 | keep training |
| 5 | 0.000222 | 0.000445 | 0.226 | 0.144 | keep training |
| 6 | 0.000344 | 0.000476 | 0.197 | 0.151 | keep training |
| 7 | 0.000159 | 0.000403 | 0.141 | 0.155 | keep training |

Read the two trends. Freeze only when both break.

