# PO × MSE 对照看板

就是个看板。PO-risk 和 serving MSE 两条趋势对着看：

- **两个都正常** → 接着 train，全量每一层 back-propagate
- **PO 崩了、MSE 没崩** → 再观察，先不冻
- **PO 崩了而且 MSE 崩了** → 这个更新策略不太行，冻住；这时才开 freeze-depth clone，对着 `PO_Dict` / `MSE_Dict` 看从哪一层开始动

崩 = 因果 MA ≥ 2× 各自的 baseline。PO-risk 的 baseline 是 D_ref 对半切；MSE 的 baseline 是预训练 MLP 在 D_ref 上的误差。不做 online-bootstrap。

何时 **update** 仍是业务逻辑。看板只给这三档对照。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset all --replay-json
```

看板：`results/layer_freeze_online_cv/index.html`。
