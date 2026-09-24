# FW-pm25-naive regime watch

**Headline:** FW-pm25-naive watch: 1 naive-routed packs; regime_alerts on naive=0, all_packs=0

| pack | router | MAE | surprise_rate | regime_alerts |
|---|---|---:|---:|---:|
| `metro_interstate` | `hgb` | 590.1 | 0.004 | 0 |
| `beijing_pm25` | `naive_last` | 15.81 | 0.008 | 0 |
| `waymo_proxy` | `ridge` | 0.2097 | 0.124 | 0 |

## Reading

- On naive-routed packs: keep last-value; `regime_alert` ⇒ re-run bakeoff/router.
- Do not promote HGB just because surprises exist — check streak first.

