# PO-risk board — electricity (NSW, ordered in time)

MA stable (below 2× baseline) → all-layer backprop. No online-bootstrap.
Large deviation → from which layer to freeze. Then read MSE_Dict next to PO_Dict.
Tabular PO-risk keeps a separate outcome model and a separate propensity model.
T=1 on the incoming batch. n_ref=10000, n_new=5000. baseline=0.0005145.

- hidden_dims = `[64, 32]`, k = 3
- n_batches = 7, n_large = 2
- ma_window = 5, all_layer_backprop = False
- **on large-deviation batches, freeze from model_0**

| t | PO-risk | MA | baseline | large | action |
|---:|---:|---:|---:|---|---|
| 0 | 9.59e-05 | 9.59e-05 | 0.000515 |  | all trainable |
| 1 | 0.000964 | 0.00053 | 0.000515 |  | all trainable |
| 2 | 0.00166 | 0.000906 | 0.000515 | yes | freeze from model_0 |
| 3 | 0.00597 | 0.00217 | 0.000515 | yes | freeze from model_0 |
| 4 | 0.000905 | 0.00192 | 0.000515 |  | all trainable |
| 5 | 9.41e-06 | 0.0019 | 0.000515 |  | all trainable |
| 6 | 6.08e-05 | 0.00172 | 0.000515 |  | all trainable |

PO_Dict (conditional on freeze-depth):

| t | layer0 | layer1 | layer2 | layer3 | freeze from |
|---:|---:|---:|---:|---:|---|
| 2 | 0.00171 | 0.00175 | 0.00175 | 0.00181 | freeze from model_0 |
| 3 | 0.00572 | 0.00589 | 0.00619 | 0.0063 | freeze from model_0 |

MSE_Dict (new-batch prediction error, same freeze-depths):

| t | layer0 | layer1 | layer2 | layer3 | MSE-best |
|---:|---:|---:|---:|---:|---|
| 2 | 0.228 | 0.202 | 0.177 | 0.172 | layer3 |
| 3 | 0.247 | 0.21 | 0.149 | 0.145 | layer3 |

Read PO_Dict and MSE_Dict. MA stable means every layer can backprop.

