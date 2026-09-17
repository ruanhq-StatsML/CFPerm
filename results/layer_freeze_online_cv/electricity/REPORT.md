# 从哪一层开始 update — electricity (NSW, ordered in time)

看板直接读 PO-risk。i* = argmin_i PO-risk：从这一层开始 update。
T=1 on the incoming batch. n_ref=10000, n_new=5000.
n_new is large so we do not online-bootstrap.

- hidden_dims = `[64, 32]`
- k = 3 → models 0…3
- n_batches = 7
- **start updating from model_1** (median i*)

| t | model_0 | model_1 | model_2 | model_3 | start from |
|---:|---:|---:|---:|---:|---|
| 0 | 0.000118 | 0.000112 | 9.15e-05 | 8.83e-05 | model_3 |
| 1 | 0.00107 | 0.00102 | 0.000763 | 0.000753 | model_3 |
| 2 | 0.00203 | 0.00225 | 0.00374 | 0.00364 | model_0 |
| 3 | 0.00598 | 0.00484 | 0.00503 | 0.0052 | model_1 |
| 4 | 0.00185 | 0.00126 | 0.00213 | 0.00234 | model_1 |
| 5 | 0.000132 | 0.000101 | 0.000202 | 0.000174 | model_1 |
| 6 | 0.000144 | 0.000195 | 0.000289 | 0.000243 | model_0 |

Read the PO-risk. Nothing else.

