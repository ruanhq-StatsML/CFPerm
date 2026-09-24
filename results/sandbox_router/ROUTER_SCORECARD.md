# FW-router scorecard

**Headline:** FW-router table={'metro_interstate': 'hgb', 'beijing_pm25': 'naive_last', 'waymo_proxy': 'ridge'}; MAE better than global-HGB on 2/3 packs (mean gain=0.5205)

Routing table: `{'metro_interstate': 'hgb', 'beijing_pm25': 'naive_last', 'waymo_proxy': 'ridge'}`

| pack | routed | MAE routed | MAE global HGB | gain |
|---|---|---:|---:|---:|
| `metro_interstate` | `hgb` | 447.9 | 447.9 | 0 |
| `beijing_pm25` | `naive_last` | 11.4 | 12.95 | 1.55 |
| `waymo_proxy` | `ridge` | 0.1025 | 0.1142 | 0.01177 |

## Reading

- gain>0 ⇒ specialist beats forced HGB on that pack.
- This is the P0 flywheel opportunity: route first, deepen later.

