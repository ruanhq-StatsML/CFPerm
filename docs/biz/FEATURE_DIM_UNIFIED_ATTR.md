# 风格归因 + 图谱归因 → 先统一成特征维度

> **口径**：风格归因和图谱归因不是两套系统。  
> **先**把两边都收成同一张 **特征维度目录**，**再**跑 L1 localization → L2 FSDS/LOGO。

---

## 1. 为什么先并维

| 分开做 | 并成特征维度 |
|--------|--------------|
| 文风一套指标、图谱一套节点，对不上工单 | 同一层 `D = {d1..dk}`，告警语言统一 |
| 风格旁路、图谱旁路，动作路由分叉 | Top-2/3 **维度** → 维内特征 → 动作 |
| 全特征直接 FSDS 噪 | 先缩维，再维内归因 |

```text
画风指标 (formal/hedge/len…)     ──┐
商户·作者·商品·订单 图节点         ──┼──► 特征维度 D ──► L1 Top维 ──► L2 维内 FSDS/LOGO
智能体 path / rag / text           ──┘
```

---

## 2. 目录怎么收

| 来源 | 特征维度（节点名） | 维内装什么 |
|------|-------------------|------------|
| 画风/文风 | `length` / `punct` / `register`（也可挂在 Author 下） | tok_len、formal、hedge… |
| 业务图 | `merchant` / `author` / `product` / `order` | 各维画像与路径特征 |
| LLM 同构别名 | `rag_retrieval` / `style_register` / `text_payload` / `agent_path` | rag_hit、style10、payload、step/tool/retry |

**关键句**：文风不是旁路产品——它是图上的 **Author / style_register** 维度（或其下的子维）。

---

## 3. 统一流水线（和 bulletin 落地同一句式）

```text
① 检测火了（OnlineRFPerm / RFPerm / domain AUC）
② 样本打到特征维度目录 D（风格指标、图节点、path/rag 都进同一层）
③ L1：维度级 mass / LOGO → Top-2/3 维 S*
④ L2：只在 S* 内 FSDS / 组内 LOCO → 文本指标或特征名
⑤ 动作：只砸 S*（刷检索 / 审路径 / 改 decoding / Top-k 审核）
```

账本纪律不变：`style_register` 的 CTR/品牌路径 **另账**，不进客服主账；偏好维与文风维可同层、分账。

---

## 4. 可跑

```bash
PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py
PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py --synth
PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py --alias graph   # 商户·作者·商品·订单原名
```

产物：`results/agod/feature_dim_attr/{summary.json,REPORT.md}`

Bulletin 总跑板仍见：`scripts/agod/bulletin_landing_attr.py`（三类场景）；本脚本专责 **并维层**。

---

## 5. 索引

- Bulletin：`LANDING_BULLETIN_POLISH.md`
- 图谱两级：`MERCHANT_AUTHOR_PRODUCT_ORDER_GRAPH_ATTR.md`
- **图谱刻画与评估**：`GRAPH_USECASE_EVAL.md`
- 周更·安全·图→FSDS：`WEEKLY_AUDIT_SAFETY_GRAPH_FSDS.md`
