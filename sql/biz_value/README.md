# Business-value SQL package

Maps **hallucination regimes** and **style/creative drift** onto operational KPIs.

| File | Role |
|------|------|
| `00_schema.sql` | dims + serve/signal/audit facts + unit economics |
| `01_hallucination_value.sql` | concept-fire → CS/refund cost + act-vs-ignore ROI |
| `02_style_drift_value.sql` | style AUC vs concept-fire split → CTR/brand workorders |
| `03_value_dashboard.sql` | exec dashboard + merge gate + on-call queue |

```bash
pip install duckdb pandas
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
```

Narrative map (CN): `docs/biz/HALLUC_STYLE_BIZ_SQL_MAP.md`
