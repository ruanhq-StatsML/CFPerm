# PO-risk 赋能 Post-Training（简要逻辑 + LaTeX）

> 主文（LaTeX）：[`AGOD_PO_Risk_PostTraining.tex`](AGOD_PO_Risk_PostTraining.tex)  
> 代码：`agod/po_risk_train.py` · `agod/po_iptw.py`  
> 对照表：[`AGOD_po_risk_train_compare.md`](AGOD_po_risk_train_compare.md)

---

## 一屏简要逻辑

```text
窗口 t 的 PO / MMD / proto / ΔPO
        │  (外部归因，不反传到 α)
        ▼
   α（模态单纯形）或 w_i（样本权）
        │  只作用到窗口 t+1（防同窗泄漏）
        ▼
  post-training 执行器：
    · Path A  样本：w ∝ √PO（RFPerm reject 后 soft IPTW）
    · Path B  模态：LR_m · step_alloc · BWD freeze · stack prior=α
```

**一句话：** PO-risk 不替换训练 loss，而是把「概念残差漂移」变成 **下一窗训练计划**——谁加 LR、谁多分 step、谁冻反传、哪些难样本抬权。

---

## 为什么叫「赋能」

| 赋能点 | 机制 | 不做什么 |
|---|---|---|
| 钱花在漂移模态 | Softmax(PO)→α→λ_m / s_m | 不在 α 上 BP |
| 省更新 FLOPs | α&lt;θ → BWD freeze（FWD 仍开） | 不声称推理延迟↓ |
| 难样本 | reject 后 √PO IPTW | 不 always-on 硬抬（易伤平静包） |
| 可辩护因果 | t 的传感器只调 t+1 | 不同窗用 Y_{t+1} 回灌权重再报持出 |

---

## 两条路径（详文 §3）

1. **Path A（样本级）** — OnlineRFPerm reject → `w=√PO` → 下一 fit；默认 soft，gated 优于 always-√。  
2. **Path B（模态级）** — 八个 metric version（`equal`…`po_next`）→ `next_step_actuators`。

经验要点（smoke）：Affec 类 $M\ge3$ 不对称包有 FLOPs / 小 Acc 空间（`po_gated` / `po_proto`）；COCO 类已饱和 2-mod 近 noop，别报效率赢。

---

## 编译

```bash
cd docs/agod && pdflatex AGOD_PO_Risk_PostTraining.tex
```
