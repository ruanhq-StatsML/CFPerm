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
