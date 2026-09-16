# Business-value SQL package

Maps **hallucination regimes** and **style/creative drift** onto operational KPIs.
Primary landing surface: **客服助手 incremental contribution** (tickets / refunds / contained sessions / ¥).

| File | Role |
|------|------|
| `00_schema.sql` | dims + serve/signal/audit facts + unit economics |
| `01_hallucination_value.sql` | concept-fire → CS/refund cost + act-vs-ignore ROI |
| `02_style_drift_value.sql` | style AUC vs concept-fire split → CTR/brand workorders |
| `03_value_dashboard.sql` | exec dashboard + merge gate + on-call queue |
| `04_cs_assistant_contribution.sql` | **客服助手贡献账**：费率对照、增量¥、月跑率、动作拆分 |
| `05_cs_net_and_weekly.sql` | 扣审计人力净贡献 + 周经营看板 |
| `06_cs_exec_dashboard.sql` | 经营总看板 + 动作推荐 |
| `07_cs_week_attribution_payback.sql` | 周归因贡献账 + 动作回本天数 |
| `08_cs_unit_econ_cumulative.sql` | 单位经济敏感度带 + 累计贡献曲线 |
| `09_cs_opportunity_breakeven.sql` | 不动作留白 + 成本盈亏平衡 |
| `10_cs_action_contain_split.sql` | 按动作拆多承接 |
| `11_cs_week_split_coverage.sql` | 周三项拆分 + 覆盖率外推 |
| `12_cs_payback_contain.sql` | 回本分子含承接 |
| `13_cs_hf_knob_bridge.sql` | HF hop 参数桥接 |
| `14_cs_marginal_day.sql` | 边际动作日贡献 |
| `15_cs_payback_contain_price_stress.sql` | 承接单价 ±50% 仅承接回本承压 |
| `16_cs_detection_delay_profit.sql` | Early detection：delay → 少拿毛/净¥ |
| `17_cs_fully_loaded_capture.sql` | 全成本净（扣审计+动作日）+ 留白捕获 + 周动作净 |
| `18_cs_ops_onepager.sql` | 值班一页纸：Before/After + 全成本 + delay + 优先动作 |
| `19_cs_weekly_fully_loaded.sql` | 周全成本净周报：周毛−审计−动作日 + WoW + 对账 |
| `20_cs_hf_hop_yen_band.sql` | HF hop/rag 重 seed 贡献带：全成本净 + retrieval 混比 |

```bash
pip install duckdb pandas
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
python3 scripts/agod/cs_assist_weekly_ops_brief.py
python3 scripts/agod/cs_assist_hop_sensitivity.py
```

- Narrative map: `docs/biz/HALLUC_STYLE_BIZ_SQL_MAP.md`
- CS contribution: `docs/biz/CS_ASSISTANT_CONTRIBUTION.md`
- Weekly ops brief: `docs/biz/CS_ASSISTANT_WEEKLY_OPS_BRIEF.md`
- Exec dashboard: `docs/biz/CS_ASSISTANT_EXEC_DASHBOARD.md`
