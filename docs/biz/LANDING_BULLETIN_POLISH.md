# 落地 Bulletin（对外一句 / 可粘贴）

## 打磨后主句（推荐）

> 将分布漂移检测与子集定位应用于 **AI 智能体连续推理路径**、**人工/合成数据质量检测与审核**，以及 **画风/文风漂移监测**；在推理质量衰减时，能把问题归因到 **文本指标与特征层面**，并直接接到刷新、审核与经营动作，具备可落地的生产价值。

更短版：

> 在智能体连续推理、数据质检审核与画风漂移三类场景落地漂移检测：把推理路径衰减定位到文本指标/特征子集，直接驱动审核与刷新动作。

英文对照：

> We apply distribution-shift detection and subset localization to continuous agent inference paths, human/synthetic data quality audit, and style/register drift monitoring—attributing inference-path degradation to text metrics and feature-level subsets with direct production actions.

---

## 相对原稿改了什么

| 原稿 | 问题 | 打磨 |
|------|------|------|
| 「将该检测应用于在…」 | 介词叠床架屋 | 改为「应用于 A、B 以及 C」 |
| 「人工生成数据」 | 易歧义（人标 vs 合成） | 写成「人工/合成数据质量检测与审核」 |
| 「成功将…归因」 | 易被读成因果宣称 | 写成「能把问题归因/定位到…层面」（子集定位，非唯一因果） |
| 「模型推理路径的衰减」 | 偏抽象 | 对齐智能体 **连续推理路径** 质量衰减 |
| 缺落地落点 | 只说到价值 | 补「接到刷新、审核与经营动作」 |

---

## 落地映射（已经在这套原型里对上）

| Bulletin 里的领域 | 本仓库落地 | 归因落到哪 |
|------------------|------------|------------|
| AI 智能体连续推理路径 | OnlineRFPerm 推理流；检索 vs 生成；Router/专家；early detection→delay→¥ | 文本/会话特征、rag_hit、意图桶、路径指标 |
| 人工/合成数据质检与审核 | RFPerm-as-Judge；周更合并门禁；PO Top-k 审核；安全拒答 | 偏好映射 vs 文风轴拆开；实例 Top-k |
| 画风/文风漂移 | style_domain_auc / FSDS 协变量侧；**另账**不进客服主账 | 长度、正式度、register 等文本指标与特征族 |
| 「文本指标和特征层面」 | 图 localization → FSDS / LOGO 维内归因 | 先 2–3 维/塔，再族内特征 |

**不是**：新 LLM 训练法、全空间唯一因果根因、做推荐系统。  
**是**：检测 → 子集定位 → 审核/刷新/对账。

---

## 可跑落地（本仓库已接线）

```bash
PYTHONPATH=. python3 scripts/agod/bulletin_landing_attr.py          # HF cache
PYTHONPATH=. python3 scripts/agod/bulletin_landing_attr.py --synth  # 离线合成
PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py --synth       # 风格+图谱并成特征维度
```

产物：
- `results/agod/bulletin_landing/{summary.json,REPORT.md}`
- `results/agod/feature_dim_attr/{summary.json,REPORT.md}`（并维层）

| 场景 | 检测 | 归因（文本指标 + 特征族） | 直接动作 |
|------|------|---------------------------|----------|
| 智能体连续推理 | OnlineRFPerm 火情 | LOGO/RF-mass：`text_hash` / `rag` / `style` / `path`；路径指标 step_depth·tool_call·retry | 刷检索 / 审工具链 / Top-k / 回滚生成侧 |
| 数据质检审核 | RFPerm-as-Judge + err ratio | 偏好轴 vs 文风轴拆开；PO Top-k | 拒合并 / 限量 / 人工复核 |
| 画风漂移 | style domain AUC | formal / hedge / length 等指标 LOGO | 改模板·decoding（**另账**） |

**并维纪律**：风格归因与图谱归因 **先收成特征维度**（`FEATURE_DIM_UNIFIED_ATTR.md`），再 L1→L2；文风 = Author/style_register 维，不是旁路产品。

---

## 可接的下一句（可选，bulletin 第二点）

> 进一步，在商户·作者·商品·订单等多维业务图上，先做图级 localization 标出问题维度，再在维度内做 FSDS 特征归因，形成可近实时的两级定位闭环。

---

## 索引

- 可跑脚本：`scripts/agod/bulletin_landing_attr.py`
- 特征维度并维：`FEATURE_DIM_UNIFIED_ATTR.md` / `scripts/agod/feature_dim_attr.py`
- **Prototype × 方法论场景清单（Justify）**：`LLM_LANDING_PROTOTYPE_METHOD_LIST.md`
- **OnlineRFPerm→LLM 推理 roadmap**：`ONLINERFPERM_LLM_INFER_ROADMAP.md`
- **画风漂移迭代**：`STYLE_DRIFT_ITER_PROTOTYPE.md` / `scripts/agod/style_drift_iter_proto.py`
- 总图：`LLM_LANDING_USECASES_BIZ.md`
- 审核 / 多塔：`SUBSET_LOCALIZATION_AUDIT_TOWERS.md`
- 周更·安全·图谱→FSDS：`WEEKLY_AUDIT_SAFETY_GRAPH_FSDS.md`
- 商户·作者·商品·订单：`MERCHANT_AUTHOR_PRODUCT_ORDER_GRAPH_ATTR.md`
