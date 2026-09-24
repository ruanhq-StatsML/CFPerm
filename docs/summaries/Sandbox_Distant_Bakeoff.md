# Sandbox · 离 AGOD 方法论的两条线

> 故意换场景、换方法；**不**接 excess / PO / soft-burn / claim / ToT。  
> 只看古典监督学习与无监督聚类的数字效果。

代码：`sandbox/forecast_bakeoff.py` · `sandbox/prompt_themes.py`  
跑法：`PYTHONPATH=. python3 scripts/run_sandbox_distant.py`

| 线 | 场景 | 方法 | 效果尺子 |
|---|---|---|---|
| A | stream pack 下一步回归 | naive / Ridge / HGB | holdout RMSE·MAE·R² + vs-naive lift |
| B | DiffusionDB prompt | TF-IDF + MiniBatchKMeans | silhouette + top terms |

读法：HGB lift≈0 → 近似随机游走；silhouette 选 k；**不**写成 AGOD 技能/因果主张。
