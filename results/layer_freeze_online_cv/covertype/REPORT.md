# PO-risk board — covertype (geographic order, class 2 vs rest)

No large deviation → all trainable. Large deviation → from which layer to freeze.
Tabular PO-risk keeps a separate outcome model and a separate propensity model.
T=1 on the incoming batch. n_ref=10000, n_new=10000. baseline=0.0003472.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 8, n_large = 1
- **large deviation present, freeze prototype still says all trainable**

| t | PO-risk | baseline | large | action |
|---:|---:|---:|---|---|
| 0 | 0.000255 | 0.000347 |  | all trainable |
| 1 | 0.000189 | 0.000347 |  | all trainable |
| 2 | 0.000523 | 0.000347 |  | all trainable |
| 3 | 0.00102 | 0.000347 | yes | all trainable |
| 4 | 0.000267 | 0.000347 |  | all trainable |
| 5 | 0.000222 | 0.000347 |  | all trainable |
| 6 | 0.000344 | 0.000347 |  | all trainable |
| 7 | 0.000159 | 0.000347 |  | all trainable |

Large-deviation batches, PO-risk conditional on freeze-depth:

| t | model_0 | model_1 | model_2 | model_3 | freeze from |
|---:|---:|---:|---:|---:|---|
| 3 | 0.000917 | 0.000911 | 0.000724 | 0.000652 | all trainable |

Read the PO-risk. Nothing else.

