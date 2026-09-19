# 图谱特征 FSDS × PO-risk 过夜迭代报告

**范围锁定（讲武德）**  
- 输入：TencentGR **图谱特征** \(X\)（已 localize 的 W1/W2 grid）  
- `W` = **时段**（W1 vs W2），**不是** treatment  
- PO-risk = 时段漂移代理 `mean(τ̂²)` → **只做特征排序先验**，交给官方 FSDS  
- **不做**：community / ego-net / ULS·PPR·GraphScan / GNN / ATE 因果宣称  

**PR：** https://github.com/ruanhq-StatsML/CFPerm/pull/72  
**日志：** `results/tencent_gr_fsds_iterate/ITERATION_LOG.md`  
**Cookbook：** `docs/tencent_gr/PO_RISK_FOR_DS.md`

---

## 1. 一句话结论

在 rare-convert、图谱特征宽度 \(d\approx 21\) 的设定下：

> **Period-PO 的 VIMP 作为 FSDS 的排序先验，略优于纯 FSDS baseline；方法差距远小于 seed 方差。**  
> DS 默认路径：`po_help_fsds` @ **k=15**，报 seed mean±std；需要冻名单用 LOO-pos majority。

---

## 2. 协议与评测口径

```
Standardize(W1) → cmean / PO 引导（可选）→ 官方 FSDS
官方 FSDS = Scaler → VarianceThreshold → SelectKBest(F) → HGB / LogReg
```

| 项 | 设定 |
|---|---|
| 选特征 | **只用 W1-train**；W2 **永不**参与 selection |
| 主指标 | W2 temporal HGB AUC（辅：AP；rare-pos 下 AP 极小） |
| 稳健性 | seeds **0/1/2** 的 mean±std |
| 数据 | `W1_localized_grid` / `W2_localized_grid_sample`（~11k / ~12k 行；W1 仅 ~4 正例） |

---

## 3. 核心数字（seeds 0/1/2，k=15 除非注明）

| 方法 | mean W2 HGB | σ | 角色 |
|---|---:|---:|---|
| **P_po_vimp → FSDS** | **0.722** | 0.170 | **默认（最佳均值）** |
| 同上 @ **k=18** | **0.724** | 0.171 | 可选微抬 |
| Z / rare-π 族 | 0.717 | **0.137** | 要稳、略让均值 |
| A baseline FSDS | 0.716 | 0.154 | 对照 |
| LOO-pos majority（冻名单） | 0.708 | **0.114** | 需要一份固定特征表时 |
| CLI `po_help_fsds` smoke | **0.720** | — | 与默认路径一致 |

代表特征故事（PO-VIMP）：`ui_pop_mismatch` 常居首位（约 35% share），随后 covisit / item·user volume 类图特征——可读、不声称因果。

---

## 4. 过夜迭代时间线（iter01→10）

| 轮次 | 主题 | 关键发现 |
|---|---|---|
| 01 | 多步选择初探 | π-stable / cmean+π 略好；**硬 corr@0.92 伤** |
| 02–03 | 官方 FSDS fuse；组合方 | soft-corr / J*-only 不稳定或伤；收敛到 **cmean+π→FSDS** |
| 04 | k / seed / MI | seed σ ≫ method；MI 无增益 |
| 05 | **接入 PO-risk** | P_po **0.722** > Z **0.720** > A **0.716** |
| 06 | DS helper + rare-π / boot-π | `po_help_select`；rare-π 最低 σ；**boot-π 非免费午餐** |
| 07 | α∈{0.3,0.5,0.7} | **完全相同**（d≈21）→ 锁 α=0.5；tight pool=k 伤 |
| 08 | k∈{10,12,15,18}；τ̂² 行滤 | **k&lt;15 伤**；τ̂²-rows ≈ P_po |
| 09 | seed-maj / LOO 稳定性 | **seed-maj2 伤均值**；LOO-pos maj 最低 σ |
| 10 | **收口** | CLI + OVERNIGHT_SUMMARY；停 20min 定时器 |

---

## 5. 锁定给 DS 的用法

### 默认（推荐）

```bash
PYTHONPATH=. python3 scripts/tencent_gr/po_help_fsds.py \
  --w1-grid results/tencent_gr_localize_fsds_time/W1_localized_grid.parquet \
  --w2-grid results/tencent_gr_localize_fsds_time/W2_localized_grid_sample.parquet \
  --select-k 15 --seed 0 \
  --out-dir results/tencent_gr_fsds_iterate/po_help_cli
```

Python：`po_help_select(g_w1, g_w2, cols, k=15)` → 把 `pool` 交给官方 FSDS。

### 可选旋钮

| 需求 | 做法 |
|---|---|
| 略抬均值 | k=18 + 纯 PO-VIMP |
| 降 seed 方差 | rare-pos 封顶 π（`Z_combined_PO_rare`） |
| 冻一份特征名单 | leave-one-positive-out majority（勿用 seed-maj2） |
| 池子宽度 | **k+3 slack**；不要收成恰好 k |

### 明确不要做

- 调 α（本网格无效）  
- 单独 bootstrap-π / 单独 avg-VIMP top-k（不经 FSDS）  
- 硬相关剪枝、J*-only 预筛当主路径  
- 把 PO-risk 说成 treatment effect / ATE  

---

## 6. 交付物索引

| 路径 | 内容 |
|---|---|
| `scripts/tencent_gr/po_risk_fsds.py` | `fit_period_po` / `po_help_select` / τ̂² filter / rare-π helpers |
| `scripts/tencent_gr/po_help_fsds.py` | **DS CLI** |
| `scripts/tencent_gr/run_fsds_multistep_iterate.py` | 过夜变体 harness |
| `docs/tencent_gr/PO_RISK_FOR_DS.md` | 短 cookbook |
| `results/tencent_gr_fsds_iterate/iterNN/FINDINGS.md` | 各轮证据 |
| `tests/test_po_risk_fsds.py` | helper 单测 |

---

## 7. 局限与下一步（暂停点）

- 正例极少 → **seed 噪声主导**；数字只作相对排序，不作业务 KPI。  
- 结论绑在当前 localize 支撑 + \(d\approx 21\)；换支撑应重跑 seed 扫描。  
- 暂停处已够 DS 手递；若续作，优先：**更大支撑 / 更多正例** 或 **业务可读特征卡**，而非再扫 α/boot-π。
