# PO-risk board — covertype (geographic order, class 2 vs rest)

MA stable (below 2× baseline) → all-layer backprop. No online-bootstrap.
Large deviation → from which layer to freeze. Then read MSE_Dict next to PO_Dict.
Tabular PO-risk keeps a separate outcome model and a separate propensity model.
T=1 on the incoming batch. n_ref=10000, n_new=10000. baseline=0.0003472.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 8, n_large = 1
- ma_window = 5, all_layer_backprop = True
- **MA stable → all-layer backprop**

| t | PO-risk | MA | baseline | large | action |
|---:|---:|---:|---:|---|---|
| 0 | 0.000255 | 0.000255 | 0.000347 |  | all trainable |
| 1 | 0.000189 | 0.000222 | 0.000347 |  | all trainable |
| 2 | 0.000523 | 0.000322 | 0.000347 |  | all trainable |
| 3 | 0.00102 | 0.000497 | 0.000347 | yes | all trainable |
| 4 | 0.000267 | 0.000451 | 0.000347 |  | all trainable |
| 5 | 0.000222 | 0.000445 | 0.000347 |  | all trainable |
| 6 | 0.000344 | 0.000476 | 0.000347 |  | all trainable |
| 7 | 0.000159 | 0.000403 | 0.000347 |  | all trainable |

PO_Dict (conditional on freeze-depth):

| t | layer0 | layer1 | layer2 | layer3 | freeze from |
|---:|---:|---:|---:|---:|---|
| 3 | 0.000852 | 0.000834 | 0.000693 | 0.000653 | all trainable |

MSE_Dict (new-batch prediction error, same freeze-depths):

| t | layer0 | layer1 | layer2 | layer3 | MSE-best |
|---:|---:|---:|---:|---:|---|
| 3 | 0.0743 | 0.0671 | 0.0429 | 0.0275 | layer3 |

Read PO_Dict and MSE_Dict. MA stable means every layer can backprop.

