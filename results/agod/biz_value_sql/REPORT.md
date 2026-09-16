# Business-value SQL demo report

## 客服助手增量贡献（主结果）

- tickets avoided: `206.0`
- refunds avoided: `130.0`
- extra sessions contained: `316.0`
- **gross incremental ¥: `16656.0`**
- **net incremental ¥ (after audit): `16590.0`**
- ¥ breakdown tickets/refunds/contain: `5150.0` / `10400.0` / `1106.0`
- monthly run-rate ¥ (gross): `71383.0`
- action split: `cs_assist_actions.json`
- traffic / weekly / net: `cs_assist_traffic.json` / `cs_assist_weekly.json` / `cs_assist_net.json`
- detail: `CS_ASSISTANT_CONTRIBUTION.md` / `docs/biz/CS_ASSISTANT_CONTRIBUTION.md`

## Hallucination cost rollup (`shop_assistant`)

- est daily save act vs ignore: `2221.428571428571` CNY

## Style / creative (`feed_caption`)

- style-only days: `12.0`
- avg CTR style-only / quiet: `0.05444444444444444` / `0.082`

## Artifacts

- `sql/biz_value/04_cs_assistant_contribution.sql`
- `sql/biz_value/05_cs_net_and_weekly.sql`
- `results/agod/biz_value_sql/CS_ASSISTANT_CONTRIBUTION.md`
