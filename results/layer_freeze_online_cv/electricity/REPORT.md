# PO-risk board — electricity (NSW, ordered in time)

No large deviation → all trainable. Large deviation → from which layer to freeze.
Tabular PO-risk keeps a separate outcome model and a separate propensity model.
T=1 on the incoming batch. n_ref=10000, n_new=5000. baseline=0.0005076.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 7, n_large = 2
- **on large-deviation batches, freeze from model_0**

| t | PO-risk | baseline | large | action |
|---:|---:|---:|---|---|
| 0 | 0.000116 | 0.000508 |  | all trainable |
| 1 | 0.000936 | 0.000508 |  | all trainable |
| 2 | 0.00183 | 0.000508 | yes | freeze from model_0 |
| 3 | 0.00598 | 0.000508 | yes | freeze from model_0 |
| 4 | 0.000905 | 0.000508 |  | all trainable |
| 5 | 1.07e-05 | 0.000508 |  | all trainable |
| 6 | 9.31e-05 | 0.000508 |  | all trainable |

Large-deviation batches, PO-risk conditional on freeze-depth:

| t | model_0 | model_1 | model_2 | model_3 | freeze from |
|---:|---:|---:|---:|---:|---|
| 2 | 0.00184 | 0.00192 | 0.00201 | 0.00204 | freeze from model_0 |
| 3 | 0.0057 | 0.00588 | 0.00632 | 0.00645 | freeze from model_0 |

Read the PO-risk. Nothing else.

