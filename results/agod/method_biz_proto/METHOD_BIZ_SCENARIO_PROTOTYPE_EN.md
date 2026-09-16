# Method × Biz-scenario × Unit-econ Prototype

## One-liner

方法论确认幻觉制度跳变（ratio≈11.7）后，客服助手相对不动作少 206.0 工单 / 130.0 退款 / 多承接 316.0；净贡献 ¥16,590，单价±20% 净带 ¥13,259–¥19,921；优先动作 model_rollback。

## Relevance chain

1. **Method (OnlineRFPerm / PO / RAG / style AUC)** — opens the regime window and routes actions.
2. **Biz scenarios (hop weak/base/strong × traffic × actions)** — turns the window into ticket/refund/containment volumes.
3. **Unit economics (±20% prices)** — turns volumes into a ¥ band under finance assumptions.

## Same

- 都对「制度/画像是否变了」敏感
- 都需要对照臂（quiet / fire_ignored / fire_acted），不能只看绝对值
- 动作路由都依赖轴拆分（concept vs covariate / RAG vs generation）

## Different

- 方法论输出：fired / ratio / AUROC / P@10 / style_auc（诊断）
- 业务情景输出：少工单 / 少退款 / 多承接 / 动作净¥（运营）
- 单位经济输出：同一量下的 ¥ 带（财务假设）
- 方法论不保证因果归因；业务账是对照代理，不是 RCT

## Relevance

OnlineRFPerm 决定「哪段日子算跳变对照窗」；hop 情景把跳变强度变成量；单位经济把量变成可对账的 ¥ 带——三层叠乘，缺一不可。

See CN full report: `docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md`
