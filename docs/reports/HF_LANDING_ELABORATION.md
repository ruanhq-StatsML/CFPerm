# Landing elaborations (文风 / RFPerm-as-Judge / RLHF·DPO / 幻觉 / Router·RAG)

This note expands the directions called out as “ready,” and points at the
two HF subset prototypes shipped under `scripts/agod/hf_landing_protos.py`.

## 1. 文风 / 画风漂移 vs preference 漂移（拆轴）

| 轴 | 统计对象 | 工具 | 动作 |
|----|----------|------|------|
| 文风 / register / 画风 | \(P(X)\) 画像（长度、正式度、标点、嵌入几何） | RF-domain AUC / FSDS covariate VIMP | 刷新风格塔、重采样、改 decoding |
| 偏好 / 正确性映射 | \(P(Y\mid X)\)（chosen/rejected、幻觉标签） | OnlineRFPerm + PO-risk | 门控 \(\sqrt{\mathrm{PO}}\)、抽检 Top-k、必要时 drop-old |

**不要**用 PSI(长度) 解释偏好掉点；**不要**把文风 AUC 当成对齐失败。

HH-RLHF demo：`style_domain_auc`（正式度 tercile）高 ⇒ 文风画像可分；
`judge_err_ratio` + hop@cut fire ⇒ 偏好映射跳了（DPO/RLHF 连更场景）。

## 2. RFPerm-as-Judge（可包成调用）

接口心智模型：

```text
ref_batch  --fit-->  shallow RF judge μ0
stream_batch_t     --OOS err-->  e_now / e_prev  ≥ γ  ?  fire
fire 时:  po_risk0_i = residual vs μ0   → Top-k 审计 / √PO 重加权
quiet 时: w = 1（退火回均匀）
```

落地：

- **离线门禁**：新 DPO 对是否相对真人偏好批次预测性退化 → RFPerm reject/accept。
- **在线连更**：RLHF/DPO 每周合并 → OnlineRFPerm + online FDR。
- **包一层**：现有 `agod.run_rfperm_stream` / `fit_online_probe` / `po_risk0_rows`；
  特征侧 HashingVectorizer（或你们的 sentence embedding）即可，无需先训大裁判模型。

## 3. 幻觉制度检测（不是事实核查）

- **制度火**：幻觉率 / 忠实度标签的 \(P(Y\mid X)\) 跳变 → OnlineRFPerm。
- **排序**：实例 `po_risk0` 只做候选优先级。
- **RAG 轴**：knowledge–answer overlap 低 ⇒ 检索支持缺口（与生成幻觉制度分开工单）。

HaluEval demo：cut 前忠实、cut 后崩溃 ⇒ hop@cut fire；`po_risk0` Top-k 排序。

## 4. 对齐 / 偏好一致性

- 人工 vs 合成偏好：RFPerm-as-Judge 小样本门禁。
- 版本合并 CI：旧 \(T=0\) / 新 \(T=1\)；RF-domain 查画像污染，PO/R 查映射变了。
- 裁判模型校准：裁判分作 \(Y\)，相对黄金集 OnlineRFPerm。

## 5. 多任务 Router + RAG 命中 → 下一步

1. 制度火（OnlineRFPerm）→ 是否该动。
2. FSDS / domain AUC 按 **task block / modality tower** → 动哪一路。
3. RAG-hit / PO 排序 → 抽检与重加权。
4. 门控 √PO 或 drop-old → 下一跳退火。

Router 信号在 HaluEval demo 的 `router_task_shift`：哪类问题质量在 post-cut 更差，
就优先刷新该专家 / 检索器。

## 6. 与「特征库 + 多层异动 + 多模态调度 + PO-OOD」的同一闭环

```text
表征/特征库几何预检
        ↓
OnlineRFPerm 制度检测
        ↓
FSDS / PO 轴诊断（P(X) vs P(Y|X) 代理）
        ↓
分层子集 / 塔 / 任务定位
        ↓
门控适应（√PO / drop-old / KV·LoRA 预算）
        ↓
退火回均匀
```

推荐 / 风控 / 交易：同一闭环换特征库与业务时钟；PO-risk 作 OOD 强度与置信排序，
不是因果根因分数。

## 7. Business SQL (must-have for value justification)

幻觉 / 文风不能停在 AUROC。把信号挂到 CS/退款/CTR/品牌/合并门禁：

- 映射说明：`docs/biz/HALLUC_STYLE_BIZ_SQL_MAP.md`
- SQL 包：`sql/biz_value/`（schema + halluc value + style value + dashboard）
- 种子跑数：`PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py`

对外价值句式：fire+acted vs fire+ignored 的日均质量成本差；style-only 日只开素材工单。

## Run

```bash
PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py --gate 1.25
pytest -q tests/test_hf_landing_protos.py
```
