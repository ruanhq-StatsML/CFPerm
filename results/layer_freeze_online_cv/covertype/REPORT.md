# 从哪一层开始 update — covertype (geographic order, class 2 vs rest)

看板直接读 PO-risk。i* = argmin_i PO-risk：从这一层开始 update。
T=1 on the incoming batch. n_ref=10000, n_new=10000.
n_new is large so we do not online-bootstrap.

- hidden_dims = `[64, 32]`
- k = 3 → models 0…3
- n_batches = 8
- **start updating from model_3** (median i*)

| t | model_0 | model_1 | model_2 | model_3 | start from |
|---:|---:|---:|---:|---:|---|
| 0 | 0.000432 | 0.000445 | 0.000315 | 0.000291 | model_3 |
| 1 | 0.000559 | 0.000764 | 0.000483 | 0.000317 | model_3 |
| 2 | 0.00132 | 0.00148 | 0.000688 | 0.000613 | model_3 |
| 3 | 0.00195 | 0.00208 | 0.00144 | 0.00118 | model_3 |
| 4 | 0.00154 | 0.0017 | 0.00156 | 0.000945 | model_3 |
| 5 | 0.000944 | 0.00141 | 0.00237 | 0.00184 | model_0 |
| 6 | 0.00201 | 0.00249 | 0.00199 | 0.0012 | model_3 |
| 7 | 0.00106 | 0.00259 | 0.00258 | 0.00162 | model_0 |

Read the PO-risk. Nothing else.

