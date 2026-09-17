# 除非有大 deviation，否则一直 trainable

看板直接读 PO-risk。没有大偏差，AnyMLP 各层都 trainable。
有大偏差，才做 freeze prototype：从第几层开始冻。
这是 conditional on 当前模型的逻辑。

表格 PO-risk **单独**维护两套 nuisance，不是 serving MLP，也不是一个 lstsq：

- outcome model \(\mu(Y\mid X)\)
- propensity \(e(T\mid X)\)，新 batch 是 \(T=1\)

\[
\varphi=(Y-\mu)(T-e),\quad \widehat\tau(X)\approx\varphi,\quad \mathrm{risk}=\mathrm{mean}(\widehat\tau^2)
\]

和 D_ref 对半切的 baseline 比，stream ≥ 2× baseline 算大 deviation。
有大偏差时，outcome 再看上该 freeze-depth MLP 的预测（conditional on 模型），propensity 仍单独拟合。

n_new 不宜过小，否则要 online-bootstrap，计算量太大。electricity 5000，covertype 10000。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset both
```

看板：`results/layer_freeze_online_cv/index.html`。
