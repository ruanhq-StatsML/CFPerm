# 从哪一层开始 update：看板直接读 PO-risk

就这一件事。k 层 AnyMLP 维护 **k+1** 个模型 `model_0 … model_k`。
`model_i` 从最上面第 i 层开始 update（i=0 全冻）。
每个新来的 batch 是 **T=1**，reference 是 **T=0、n_ref=10000**。
CosineAnnealingLR 走完，**直接读** 这个模型的 streaming PO-risk。
i\* = argmin_i PO-risk：从这一层开始 update。没有别的。

\[
\varphi=(Y-\mu)(T-e),\quad \widehat\tau(X)\approx\varphi,\quad \mathrm{risk}=\mathrm{mean}(\widehat\tau^2)
\]

μ 就是该 `model_i` 的预测。

## Batch 不宜过小

PO-risk 是 (ref ∪ new) 上的**点估计**。incoming `n_new` 太小，T=1 那一块很薄，argmin_i 会抖，看板读不稳。

要在小 batch 上读稳，就得对每个 `model_0…model_k` 做 **online-bootstrap**（新 batch 重采样 B 次，每次再算一遍 PO-risk，B 通常上百）。k+1 个模型 × B 次，还叠在每次在线更新后面，计算量太大，看板不干这个。

所以 `n_new` 放到和 `n_ref` 一个量级（几千到一万），**直接读 PO-risk**，不做 bootstrap。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset both
```

看板：`results/layer_freeze_online_cv/index.html`。
electricity 默认 incoming batch 5000；covertype 10000。
