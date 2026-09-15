# 幻觉制度 × 文风/画风漂移 → 业务逻辑与 SQL 价值论证

目标：把 OnlineRFPerm / PO-risk / RF-domain **直接挂到可运营的业务动作与钱**，并用可维护 SQL 持续 justify。  
不是把算法“翻译成故事”，而是：**信号 → 动作落地 → KPI 增量 → 单位经济 → 贡献账**。

SQL 包：`sql/biz_value/`  
可跑种子：`PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py`  
**客服助手主账**：`docs/biz/CS_ASSISTANT_CONTRIBUTION.md`  
种子报告：`results/agod/biz_value_sql/REPORT.md`

---

## 0. 客服助手增量贡献（主落地）

业务只认这三样：少工单、少退款、多机器人承接。钱从单位经济乘出来。

| 对照 | 含义 |
|------|------|
| `fire_acted` | 制度跳变期，客服助手侧动作落地（检索刷新 / 回滚 / Top-k 审计） |
| `fire_ignored` | 同条件不动作 |
| 增量 | acted 相对 ignored 的费率差 × acted 会话量 × 单价 |

种子可复现结果（28 天窗、跳变后 7 acted / 7 ignored）：

| KPI | 数值 |
|-----|------|
| 少工单 | **206** 单 |
| 少退款 | **130** 单 |
| 多承接 | **316** 次 |
| 区间增量贡献 | **¥16,656**（工单¥2,825 + 退款¥5,520 + 承接¥627） |
| 每千会话 | **¥11,897** |
| 30 天跑率 | **¥71,383**/月（demo 日 200 会话） |
| 生产中等（日 1 万会话） | **¥1.92M**/月 |
| 动作拆分 | rollback ¥3,860 / audit ¥2,569 / retrieval ¥2,542 |

SQL：`sql/biz_value/04_cs_assistant_contribution.sql`  
视图：`vw_cs_assist_exec_summary` / `vw_cs_assist_rate_compare` / `vw_cs_assist_daily_ledger`

对外一句：

> 客服助手在幻觉制度跳变的动作日里，相对同条件不动作，少了 113 单工单、69 单退款，多承接 179 次会话，贡献约 ¥8,972；日 1 万会话约 ¥1.92M/月。

---

## 1. 为什么必须 map 到业务（而不是停在 AUROC）

| 算法信号 | 若不上业务 | 会上业务之后 |
|---------|------------|--------------|
| `rfperm_fire`（concept） | “模型漂了”的空告警 | **触发客服助手动作** → 工单/退款/承接增量 |
| `po_risk0` Top-k | 抽象排序 | **审计人力 ROI**（P@10、确认差例） |
| `style_domain_auc` | 误当成质量挂了 | **只动素材配比 / decoding**，避免误重训偏好头 |
| `judge_err_ratio` | 离线表 | **对齐发布决策**（是否放行偏好更新） |

价值论证的正确句子：

> 相对同条件不动作，动作落地少了 \(T\) 单工单、\(R\) 单退款、多承接 \(C\) 次，贡献 ¥\(V\)（见 `vw_cs_assist_exec_summary`）。

而不是：

> 幻觉检测 AUROC=0.9，所以有业务价值。

---

## 2. 幻觉制度 → 业务逻辑（购物助手 / 客服 / RAG）

### 2.1 业务对象

- **Surface**：`shop_assistant` / `cs_bot` / `rag_qa`
- **坏结果**：错误承诺、错商品属性、错政策 → **工单、退款、信任掉点、转化掉点**
- **算法角色**：检测 \(P(Y_{\text{坏}}\mid X)\) 是否进入新制度；**不是**事实核查器

### 2.2 动作路由（必须拆桶）

```text
RFPerm concept fire?
 ├── avg_rag_hit < 0.35  → retrieval_gap
 │                         动作：刷新索引/检索，先别全量 SFT
 └── avg_rag_hit ≥ 0.35  → generation_regime
                           动作：√PO 门控 / 回滚 model_version + po_risk0 Top-k 人工审计
RFPerm quiet 但 ticket_rate 高 → labeler_or_product_other
                           动作：查标注漂移或产品政策，不是模型制度
```

SQL：`vw_halluc_root_bucket` → `vw_halluc_workorders`  
单位经济表：`dim_value_assumption`（`cs_ticket_cost` / `refund_unit_cost`）  
信号表：`fct_shift_signal`；审计队列：`fct_audit_queue`

### 2.3 钱怎么算（单位经济表维护，不写死在查询里）

`dim_value_assumption` 维护：

| metric | 含义 | 例 |
|--------|------|----|
| `cs_ticket_cost` | 一单客服成本 | 25 CNY |
| `refund_unit_cost` | 一单退款期望损失 | 80 CNY |

日成本：

\[
\text{quality\_cost} = n_{\text{tickets}}\cdot c_{\text{cs}} + n_{\text{refunds}}\cdot c_{\text{refund}}
\]

ROI 视图：`vw_halluc_roi_rollup`  
比较 **fire+acted** vs **fire+ignored** 的日均成本差 = 可对外讲的节省代理。

### 2.4 和 HH/HaluEval 原型怎么接

原型 `results/agod/hf_landing/halu_regime_rag.json`：

- `hop_at_cut.fired=1` → 写入 `fct_shift_signal.signal_type='rfperm_fire', axis='concept'`
- 实例 `po_risk0` Top-k → `fct_audit_queue`
- `mean_rag_hit_top10` → 与 `avg_rag_hit` 一起决定 `retrieval_gap` vs `generation_regime`

---

## 3. 文风 / 画风漂移 → 业务逻辑（Feed / 广告 / 商详图）

### 3.1 业务对象

- **Surface**：`feed_caption` / `ad_creative` / `shop_image_gen` / `video_cover`
- **坏结果**：CTR 掉、完播掉、品牌投诉、素材疲劳
- **算法角色**：检测 \(P(X)\) 画像（风格簇、正式度、长度、画风 cluster）是否变；**默认不解释为对齐失败**

### 3.2 动作路由（防误伤）

```text
style_domain_auc ≥ 0.70 且 concept_fire=0  → style_only_portrait
    动作：只调素材配比 / decoding / 风格采样，不重训偏好头

concept_fire=1 且 style_auc < 0.60         → concept_only_quality
    动作：RFPerm-as-Judge 门禁 + √PO；风格采样器先别动

两轴都动                                   → joint_style_and_concept
    动作：拆两张工单，禁止混成“模型坏了”
```

SQL：`vw_style_vs_concept_split` → `vw_style_workorders`

### 3.3 钱怎么算

| metric | 含义 |
|--------|------|
| `value_per_ctr_point` | +0.01 CTR 的 GMV/广告变现代理 |
| `brand_incident_cost` | 一次品牌标记的期望成本 |

描述性（非因果）对比：`vw_style_roi_rollup`  
style_only 日的 brand_cost / CTR vs quiet 日 —— 用于 **证明“风格工单该不该开”**，不是证明某张图导致 GMV。

### 3.4 和 HH 原型怎么接

`hh_judge_style.json`：

- `style_domain_auc` → `fct_shift_signal.signal_type='style_domain_auc'`
- `judge_err_ratio` → `judge_err_ratio`（合并门禁）
- `hop_at_cut.fired` → concept fire（偏好映射轴）

合并门禁 SQL：`vw_alignment_merge_gate`  
→ `BLOCK_MERGE` / `ALLOW_WITH_STYLE_REBALANCE` / `ALLOW`

---

## 4. 一张图：从信号到钱

```text
serve_event (曝光/生成)
    ⊕ shift_signal (RFPerm/PO/style AUC)
    ⊕ audit_queue (po_risk0 Top-k)
    ⊕ value_assumption (单位经济)
            ↓
   root_bucket / axis_split     ← 业务路由（干什么）
            ↓
   value_daily / roi_rollup     ← 价值论证（值多少）
            ↓
   workorders / merge_gate      ← 值班与发布门禁
```

---

## 5. 维护清单（谁改什么）

| 表/视图 | Owner | 刷新 |
|--------|-------|------|
| `dim_value_assumption` | 业务 FP / 数据 | 季度审 |
| `fct_serve_event` | 数据工程 | T+1 |
| `fct_shift_signal` | 算法 on-call | 近实时/小时 |
| `fct_audit_queue` | 质检 | 日内 |
| `vw_*_workorders` | 值班看板 | 自动 |
| `vw_*_roi_rollup` | 周会价值页 | 自动 |

---

## 6. 对外怎么讲（一句话模板）

**客服助手（主）**：  
“幻觉制度跳变期动作落地相对不动作：少 \(T\) 工单、少 \(R\) 退款、多承接 \(C\) 会话，贡献 ¥\(V\)；月跑率见 `vw_cs_assist_exec_summary.monthly_runrate_yen`。”

**文风/画风**：  
“style_domain_auc 升高且无 concept-fire 的日子，只调素材配比/decoding，不重训偏好头；品牌标记成本与 CTR 代理见 `vw_style_value_daily`。”

---

## 7. 跑起来

```bash
pip install duckdb pandas
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
# → docs/biz/CS_ASSISTANT_CONTRIBUTION.md
# → results/agod/biz_value_sql/REPORT.md
# → vw_cs_assist_exec_summary / vw_halluc_roi_rollup / vw_style_roi_rollup
```

### 种子结论（justify 模板）

| 主题 | 对外一句话 |
|------|------------|
| **客服助手** | 少 113 工单 / 69 退款 / 多承接 179，增量 **¥16,656**；日 1 万会话 ≈ **¥1.92M/月** |
| 文风 | style-only 日 CTR 低于 quiet，品牌标记成本高于 quiet → 调素材配比，不重训偏好 |
| 对齐发布 | `judge_err_ratio≥1.5` 或 concept-fire → 暂缓偏好更新上线 |
