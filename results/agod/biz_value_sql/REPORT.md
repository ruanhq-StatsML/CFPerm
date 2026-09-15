# Business-value SQL demo report

Seeded DuckDB run for hallucination + style/creative surfaces.

## Hallucination ROI (`shop_assistant`)

- fire acted days: `7.0`
- fire ignored days: `7.0`
- avg quality cost fire+acted: `710.0`
- avg quality cost fire+ignored: `1211.4285714285713`
- **est daily save act vs ignore: `501.42857142857133` CNY**

Narrative: acting on concept-fire (√PO / rollback / Top-k audit) vs ignoring
reduces CS+refund cost on the same fire regime — this is the value sentence.

## Style / creative ROI (`feed_caption`)

- style-only days: `12.0`
- concept-only days: `0.0`
- joint days: `6.0`
- avg CTR style-only / quiet: `0.059444444444444446` / `0.082`
- avg brand cost style-only / quiet: `637.5` / `90.0`

Narrative: high `style_domain_auc` without concept-fire → open **creative mix**
ticket only; do not retrain preference head.

## Artifacts

- DuckDB: `results/agod/biz_value_sql/biz_value.duckdb`
- JSON: `roi_halluc.json`, `roi_style.json`, `dashboard.json`, `merge_gate_sample.json`
- SQL package: `sql/biz_value/`
- Mapping doc: `docs/biz/HALLUC_STYLE_BIZ_SQL_MAP.md`
