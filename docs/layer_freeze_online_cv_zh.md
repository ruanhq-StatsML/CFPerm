# 除非有大 deviation，否则一直 trainable

看板直接读 PO-risk。没有大偏差，AnyMLP 各层都 trainable。
有大偏差，才做 freeze prototype：从第几层开始冻。
这是 conditional on 当前模型的逻辑。

**何时 update 是业务逻辑**，纯统计 justify 不了。缩小 batch 只会让 PO-risk 旗标抖，不是 update 开关。

表格 PO-risk **单独**维护两套 nuisance，不是 serving MLP，也不是一个 lstsq：

- outcome model \(\mu(Y\mid X)\)
- propensity \(e(T\mid X)\)，新 batch 是 \(T=1\)

\[
\varphi=(Y-\mu)(T-e),\quad \widehat\tau(X)\approx\varphi,\quad \mathrm{risk}=\mathrm{mean}(\widehat\tau^2)
\]

和 D_ref 对半切的 baseline 比，stream ≥ 2× baseline 算大 deviation。
有大偏差时，outcome 再看上该 freeze-depth MLP 的预测（conditional on 模型），propensity 仍单独拟合。

n_new 不宜过小：不是“小了才能看到何时 update”，而是点估计会抖，纯统计更不能当 update 开关。何时 update 是业务（政策切、队列换、域换）。PO-risk 只在业务已经要动的时候，读能不能全开、还是从哪层冻。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset both --skip-freeze
```

看板：`results/layer_freeze_online_cv/index.html`。
