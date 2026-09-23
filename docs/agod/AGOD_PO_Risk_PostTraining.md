# PO-risk Metrics → 加速后训练（落地 procedure）

> 主文（LaTeX）：[`AGOD_PO_Risk_PostTraining.tex`](AGOD_PO_Risk_PostTraining.tex)  
> 代码：`agod/po_risk_train.py` · `agod/po_iptw.py`  
> 对照表：[`AGOD_po_risk_train_compare.md`](AGOD_po_risk_train_compare.md)

**目标锁定：** 只加速后训练 / 在线适配（更新 FLOPs、step、wall-clock、到 Acc★ 的窗数）。Acc/MSE 是约束，不是主指标。

---

## Metrics 落地 = 传感器 → 执行器 → 速度 KPI

```text
窗口 t 算 PO(/MMD/proto/ΔPO)
        → map → α 或 w_i
        → 执行器(t+1)：freeze / LR / step_alloc / √PO-IPTW
        → 报 FLOPs_rel、steps、wall-clock、T(Acc★)
        → 约束：ΔAcc ≥ −ε 才算落地成功
```

---

## 主 KPI（按这个顺序报）

| KPI | 角色 |
|---|---|
| FLOPs_rel / step 预算 / wall-clock / T(Acc★) | **主：加速** |
| ΔAcc、ΔMSE vs equal | **约束** |
| H(α)、hard-rank | 诊断 |

Acc 涨了但 FLOPs 不降 → **不算**加速落地。

---

## 共享 procedure P0（每个 metric version 都走）

1. **Sense** — 必算 PO_m；可选 MMD / proto / ΔPO  
2. **Map** — `metric_to_alpha(v, …)`（见下表选 v）  
3. **Actuate（加速）** — freeze（主砍 FLOPs）+ step 倾倒 + LR；stack prior 次要  
4. **Train** t+1（FWD 始终开）  
5. **Score** — FLOPs↓ 且 ΔAcc≥−ε → OK；否则调 θ_fr / floor / 换 version  

Path A：RFPerm reject 后 `w∝√PO`（gated）；calm 窗保持 uniform。

---

## 各 metric 怎么落地（focus = 速度旋钮）

| version | 何时开火 | 加速落地 |
|---|---|---|
| `equal` | baseline / Acc 已饱和 | 不定速；定 Acc★ 条 |
| `po_soft` | PO 尖峰、M≥3 | LR+steps→顶 α；冻尾部 |
| `po_minus_cov` | MMD 高、概念平 | **别**给纯 X 漂买 step；冻 cov 重模态 |
| `po_gated` | 噪 vs 漂不清 | Affec 类主 FLOPs 砍法 |
| `po_proto` | 真质心动 + PO | step 集中；约束更好过 |
| `po_delta` | ΔPO 上升 | Acc 掉之前重分预算 → T(Acc★)↓ |
| `po_budget` | 塔必须常开 | 只软 LR；**不冻** |
| `po_next` | 有稳定 hist | 用预报 α 更早 freeze/realloc |

---

## 经验（加速透镜）

- Affec M=5：`po_gated` FLOPs≈0.70 + Acc 近似平 → **可报加速**  
- COCO 等 2-mod 已高 Acc：FLOPs≈1 → **别报加速**；`equal`/`po_budget`  
- 真要加速：把 `step_alloc` 接到真实 optimizer step，不要只调 LR  

```bash
cd docs/agod && pdflatex AGOD_PO_Risk_PostTraining.tex
```
