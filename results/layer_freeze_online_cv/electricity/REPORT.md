# PO × MSE board — electricity (NSW, ordered in time)

PO broken and MSE broken → freeze. PO broken, MSE holds → watch. Both quiet → keep training.
No online-bootstrap. Freeze-depth PO_Dict / MSE_Dict only on freeze hops.
T=1 on the incoming batch. n_ref=10000, n_new=5000.
po_base=0.0005145. mse_base=0.121.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 7, n_watch = 4, n_freeze = 0
- **watch**

| t | PO-risk | PO MA | MSE | MSE MA | action |
|---:|---:|---:|---:|---:|---|
| 0 | 9.59e-05 | 9.59e-05 | 0.118 | 0.118 | keep training |
| 1 | 0.000964 | 0.00053 | 0.204 | 0.161 | keep training |
| 2 | 0.00166 | 0.000906 | 0.229 | 0.184 | keep training |
| 3 | 0.00597 | 0.00217 | 0.158 | 0.177 | watch |
| 4 | 0.000905 | 0.00192 | 0.196 | 0.181 | watch |
| 5 | 9.41e-06 | 0.0019 | 0.262 | 0.21 | watch |
| 6 | 6.08e-05 | 0.00172 | 0.124 | 0.194 | watch |

Read the two trends. Freeze only when both break.

