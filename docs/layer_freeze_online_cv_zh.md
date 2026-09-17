# 除非有大 deviation，否则一直 trainable

看板直接读 PO-risk。没有大偏差，AnyMLP 各层都 trainable，**全量每一层 back-propagate**。
有大偏差，才做 freeze prototype：从第几层开始冻。
这是 conditional on 当前模型的逻辑。

**何时 update 是业务逻辑**，纯统计 justify 不了。

表格 PO-risk **单独**维护两套 nuisance，不是 serving MLP，也不是一个 lstsq：

- outcome model \(\mu(Y\mid X)\)
- propensity \(e(T\mid X)\)，新 batch 是 \(T=1\)

\[
\varphi=(Y-\mu)(T-e),\quad \widehat\tau(X)\approx\varphi,\quad \mathrm{risk}=\mathrm{mean}(\widehat\tau^2)
\]

和 D_ref 对半切的 baseline 比，stream ≥ 2× baseline 算大 deviation。
有大偏差时，outcome 再看上该 freeze-depth MLP 的预测（conditional on 模型），propensity 仍单独拟合。

raw 点估计在小 `n_new` 上会抖。**不做 online-bootstrap**：要对每个小 batch 重复推 MLP，顶不住。
改做因果 moving average（窗口约覆盖 1000 条流，`n_new=20` → window 50）。

**MA 稳定（不过 2× baseline）→ 可以放心全量每一层 back-propagate。** 这是预期读法。
online-bootstrap justify 不了：每个小 batch 重复推 MLP 顶不住。

大 hop 上 freeze-depth 同时记两套数组，从哪一层开始明显变动就对着 MSE 读：

```
PO_Dict  = {layer0: np.array(), layer1: np.array(), …, layer_k: np.array()}
MSE_Dict = {layer0: np.array(), layer1: np.array(), …, layer_k: np.array()}
```

`layer i` = `model_i`（只训 top i 组）。安静段是 NaN，不在小 `n_new` 上开 k+1 个 clone。

换表 airlines（航班延误，时间序）。n_ref 仍是 10000，n_new 收到 20 只做 MA gate，不做 clone。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset all --replay-json
```

看板：`results/layer_freeze_online_cv/index.html`。
