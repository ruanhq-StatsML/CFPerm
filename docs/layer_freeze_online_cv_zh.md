# PO × MSE 对照：评估逻辑

看板的输出是：**该冻结哪一层的 training**。只在 PO 和 MSE 都崩时才开 freeze-depth；`model_i` 只训 top i，其余层停训。

崩 = 因果 MA ≥ 2× 各自的 baseline。不做 online-bootstrap。

## 三条序列

**PO-risk** 问 \(P(Y\mid X)\) 有没有 hop（新 batch \(T=1\) vs \(D_{\mathrm{ref}}\) \(T=0\)）。nuisance 是独立的 RF：\(\mu(Y\mid X)\) 和 \(e(T\mid X)\)，不是 serving MLP，也不是 lstsq。RF 的 PO-risk **不该突然崩**；更可能先崩的是模型的 serving MSE。这是重点。PO 崩了、模型没崩 → **再观察**，先不冻。

**serving MSE** 问当前 MLP 在新 batch 上还付不付得起房租。评估在 update 之前。MSE 先崩、PO 还安静，才是主路径。

**MMD² of X** 问 \(P(X)\) 有没有动。只在「MSE 崩了但 PO 正常」这一格需要它：那一格**不是 concept drift**。

## MMD 口径（钉死）

用的是 **MMD²(\(X_{\mathrm{new}}, X_{\mathrm{ref}}\))**。跟 PO-risk 同一份 \(T=0\) reference batch。RBF、median bandwidth 钉在 \(D_{\mathrm{ref}}\) 上，和 stream 同一 \(\sigma\)。baseline 是 \(D_{\mathrm{ref}}\) 对半切。

另外两种口径**不用**，问题不一样：

| 口径 | 在问什么 | 为什么不是这块看板 |
|---|---|---|
| vs \(D_{\mathrm{ref}}\)（本看板） | 新 batch 相对固定参考分布的 \(P(X)\) | 和 PO-risk 的 \(T=0\) 对齐 |
| 跟之前所有 batch 的 pairwise MMD 均值 | 新 batch 像不像**最近的流** | 慢漂离 \(D_{\mathrm{ref}}\) 时 pairwise 仍可能安静 |
| vs 上个 batch 那一层的 representation | freeze-depth 特征有没有动 | 那是层表征问题，不是 \(X\) 的 covariate shift |

| PO-risk | serving MSE | MMD(\(X_{\mathrm{new}}, X_{\mathrm{ref}}\)) | 读法 | 动作 |
|---|---|---|---|---|
| 安静 | 安静 | — | 没有 hop，模型还在拟合 | **接着 train**（一层都不冻） |
| 崩了 | 没崩 | — | 机制可能动了，误差还在线内；RF 下这一格应少见 | **再观察**（先不冻） |
| 崩了 | 崩了 | — | hop 可见而且现模型也崩了 | **冻住那一层的 training** |
| 没崩 | 崩了 | 崩了 | 不是 \(P(Y\mid X)\)，是 \(P(X)\) | **X shift**，读 vs-ref MMD，不按 concept 去冻层 |
| 没崩 | 崩了 | 没崩 | 不是 concept 也不是 X | **tricky**，先不冻 |

冻层只服务「两都崩」。MSE-only 去冻层会把 covariate shift 当成 concept drift。

何时 **update** 仍是业务逻辑。

## 闭环：OnlineRFPerm 标 WHEN

onset 探针就是冻住参考窗的 RF：

```python
rf = RandomForestRegressor().fit(X_ref, Y_ref)
pred = rf.predict(np.asarray(X_new))
T = mean((Y_new - pred) ** 2) - E_ref
```

last-two hop 是相邻窗的 MSE 比。第一个 hop 就是 **shift-onset**。MMD 同样只对 `X_new` vs `X_ref`。

看板读 WHAT：

| DGP | P(X) | P(Y\|X) | OnlineRFPerm | MMD²(X_new, X_ref) | 看板该看到 |
|---|---|---|---|---|---|
| gradual concept | 固定 | β 慢慢翻 | onset 在 labeled 附近 | 安静 | MSE/PO 动，不是 X shift |
| gradual covariate | μ 慢慢走 | 固定 | onset 在 labeled 附近 | 过线 | MSE 崩、PO 安静 → X shift |

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --justify
```
