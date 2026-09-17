# PO × MSE 对照：评估逻辑

看板的输出是：**该冻结哪一层的 training**。只在 PO 和 MSE 都崩时才开 freeze-depth；`model_i` 只训 top i，其余层停训。

崩 = 因果 MA ≥ 2× 各自的 baseline。不做 online-bootstrap。

## 三条序列

**PO-risk** 问 \(P(Y\mid X)\) 有没有 hop（新 batch \(T=1\) vs \(D_{\mathrm{ref}}\) \(T=0\)）。这是 concept drift 探测器，不是现模型好不好用。

**serving MSE** 问当前 MLP 在新 batch 上还付不付得起房租。评估在 update 之前。

**MMD² of X** 问 \(P(X)\) 有没有动。RBF、median bandwidth 钉在 \(D_{\mathrm{ref}}\) 上，和 stream 同一 \(\sigma\)。baseline 是 \(D_{\mathrm{ref}}\) 对半切。只在「MSE 崩了但 PO 正常」这一格需要它：那一格**不是 concept drift**，先看是不是 X 在变。

| PO-risk | serving MSE | MMD(X) | 读法 | 动作 |
|---|---|---|---|---|
| 安静 | 安静 | — | 没有 hop，模型还在拟合 | **接着 train**（一层都不冻） |
| 崩了 | 没崩 | — | 机制可能动了，误差还在线内 | **再观察**（先不冻） |
| 崩了 | 崩了 | — | hop 可见而且现模型也崩了 | **冻住那一层的 training** |
| 没崩 | 崩了 | 崩了 | 不是 \(P(Y\mid X)\)，是 \(P(X)\) | **X shift**，读 MMD，不按 concept 去冻层 |
| 没崩 | 崩了 | 没崩 | 不是 concept 也不是 X | **tricky**，先不冻 |

冻层只服务「两都崩」：看板告诉你该冻结哪一层的 training。MSE-only 去冻层会把 covariate shift 当成 concept drift。

何时 **update** 仍是业务逻辑。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset all --replay-json
```
