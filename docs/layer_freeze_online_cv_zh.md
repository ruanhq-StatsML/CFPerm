# 哪一层该 freeze：streaming PO-risk 做在线训练的交叉验证

先把 AnyMLP 放在表格流上。k 层就维护 **k+1** 个模型 `model_0 … model_k`。
`model_i` 只训练最上面 i 层（i=0 全冻）。每个新来的 batch 是 **T=1**，
reference 是 **T=0、n_ref=10000**。CosineAnnealingLR / AdamW（或 SGD）在未冻的层上走几步，
然后看 streaming PO-risk：

\[
\varphi=(Y-\mu)(T-e),\quad \widehat\tau(X)\approx\varphi,\quad \mathrm{risk}=\mathrm{mean}(\widehat\tau^2)
\]

μ 是这个 freeze-depth 的 AnyMLP 预测。i\* = argmin_i PO-risk：更新到这一层，下面冻住。
更深会把 reference 的 P(Y|X) 冲掉；更浅跟不上 T=1 的 hop。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset both
```

看板：`results/layer_freeze_online_cv/index.html`。
默认表：electricity（真时间漂移）；更大的 covertype（地理顺序，n≫10000）。
Brier 只是 side gauge，不是 CV。
