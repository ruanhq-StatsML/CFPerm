# Business-value SQL demo report

## 客服助手增量贡献（主结果）

- tickets avoided: `113.0`
- refunds avoided: `69.0`
- extra sessions contained: `179.0`
- **incremental ¥ realized: `8972.0`**
- ¥ breakdown tickets/refunds/contain: `2825.0` / `5520.0` / `627.0`
- incremental ¥ / 1k sessions: `6408.0`
- monthly run-rate ¥: `38449.0`
- action split: `cs_assist_actions.json`
- traffic scenarios: `cs_assist_traffic.json`
- detail: `CS_ASSISTANT_CONTRIBUTION.md` / `docs/biz/CS_ASSISTANT_CONTRIBUTION.md`

## Hallucination cost rollup (`shop_assistant`)

- est daily save act vs ignore: `1192.142857142857` CNY

## Style / creative (`feed_caption`)

- style-only days: `12.0`
- avg CTR style-only / quiet: `0.05444444444444444` / `0.082`

## Artifacts

- `sql/biz_value/04_cs_assistant_contribution.sql`
- `results/agod/biz_value_sql/CS_ASSISTANT_CONTRIBUTION.md`
- `cs_assist_increment.json` / `cs_assist_actions.json` / `cs_assist_traffic.json`
