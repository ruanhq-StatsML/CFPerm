# RSI 转向：PO-risk 适应效率（离开相邻窗 transfer）

> **换探索方向。** 前四轮（null / ECE / block-CI / reservoir）钉在 sample-chunk
> transfer 探针。本转向问另一件事：

> **在 OnlineRFPerm 拒绝批上，ref / probe / refit 谁用更少适应 FLOPs 买到硬样本排序或 sig-MSE 下降？**

不锁广告；数据 = metro / PM25 / stocks / waymo_proxy（已有 `summary.json`）。

相关：[`AGOD_po_ref_vs_refit.md`](../agod/AGOD_po_ref_vs_refit.md) ·
[`RSI_Iteration_Log.md`](./RSI_Iteration_Log.md)

---

## 1. 统计主张

| 量 | 含义 |
|---|---|
| Spearman(PO, truth_hard) | 硬行排序技能（主） |
| sig-only next-MSE | IPTW 下游（次） |
| `relative_flops` | ref≈0；probe=每步重拟合；refit=仅拒绝窗重拟合 |
| `rank_eff` | Δρ vs ref / FLOPs |
| `mse_eff` | (unif−mode) MSE_sig / FLOPs |

**关键点：** 排序赢 ≠ MSE 赢；probe 的 always-on 成本常被低估；refit 拒绝稀疏时可能 **rank_eff 更高**。

## 2. 跑法

```bash
PYTHONPATH=. python3 scripts/run_po_eff_scorecard.py \
  --summary results/agod_po_ref_vs_refit/summary.json \
  --out results/agod_po_eff
```

代码：`agod/po_eff.py` · `scripts/run_po_eff_scorecard.py`

## 3. 后续发散（本方向）

| ID | 点子 | 状态 |
|---|---|---|
| P1 | gate duty × flops → 期望适应预算 | **Iter6 落地** |
| P2 | cbrt vs sqrt IPTW 的 mse_eff | **Iter7 落地**（同 FLOPs，软权重） |
| P3 | image-OOD bench 接入同一 FLOPs 尺 | **关掉**：PO 不适合 image-OOD（负结果）；省 FLOPs |
| P4 | Grad-RFPerm freeze 闭环 MSE–FLOPs 并表 | **Iter9 落地** |
| P5 | Tomorrow demo 双栏图 | **Iter8 落地** |
| P6 | 软权重该不该烧 + 算力账拆分 | **Iter19 落地** |

### P6 读法
加权 FLOPs≈0；α 不是省算力旋钮。默认不烧；beat uniform 才 `BURN_*`；被迫 gate → `SOFTEN_ONLY`（∛）。
