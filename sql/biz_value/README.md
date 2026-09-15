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
| `05_cs_net_and_weekly.sql` | 扣审计人力净贡献 + 周经营看板 |

```bash
pip install duckdb pandas
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
# → docs/biz/CS_ASSISTANT_CONTRIBUTION.md
```

Narrative map (CN): `docs/biz/HALLUC_STYLE_BIZ_SQL_MAP.md`  
CS contribution (CN): `docs/biz/CS_ASSISTANT_CONTRIBUTION.md`
