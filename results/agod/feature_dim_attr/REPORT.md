# 特征维度统一归因 — 风格 + 图谱先并维

> 风格归因与图谱归因 **先整合为特征维度**，再跑同一套 L1→L2。

## 并维口径

| 来源 | 收成的特征维度 |
|------|----------------|
| 画风/文风指标 | length / punct / register（Style 节点下的子维） |
| 商户·作者·商品·订单图 | merchant / author / product / order |
| LLM 同构别名 | rag_retrieval / style_register / text_payload / agent_path |

## Graph → 特征维度

- Top dims: `['style_register', 'rag_retrieval']`
- LOGO share: `{'rag_retrieval': 0.36388518340124837, 'style_register': 0.3650886991605228, 'text_payload': 0.271025666857788, 'agent_path': 0.0}`
- RF mass: `{'rag_retrieval': 0.10732726148939956, 'style_register': 0.20514532624367643, 'text_payload': 0.04977904629528641, 'agent_path': 0.6377483659706374}`
- Actions: ['改模板·decoding（另账，不进客服主账）', '刷检索 / 换索引（别先回滚整模）']
- Within: `{'style_register': [{'rank': 1, 'name': 'char_len', 'vimp': 0.10338999617543976}, {'rank': 2, 'name': 'tok_len', 'vimp': 0.10175533006844184}, {'rank': 3, 'name': 'avg_word', 'vimp': 0.0}, {'rank': 4, 'name': 'qmark', 'vimp': 0.0}, {'rank': 5, 'name': 'bang', 'vimp': 0.0}], 'rag_retrieval': [{'rank': 1, 'name': 'rag_hit', 'vimp': 0.10732726148950689}]}`

## Style → 特征维度（不是旁路产品）

- Top dims: `['register', 'length']`
- LOGO/RF: logo=`{'length': 0.0, 'punct': 0.0, 'register': 0.0}` rf=`{'length': 0.3749999999996251, 'punct': 0.074999999999925, 'register': 0.54999999999945}`
- Actions: ['只刷新维度 register', '只刷新维度 length']

## Audit → 特征维度（偏好维 vs 文风维）

- Top dims: `['style_register', 'preference_text']`
- share: `{'preference_text': 0.0, 'style_register': 0.9999999928714723}` / `{'preference_text': 0.9780452274813715, 'style_register': 0.021954772517628538}`
- Actions: ['改模板·decoding（另账，不进客服主账）', '只刷新维度 preference_text']

## 一句

先并到特征维度，再 localization；风格与图谱是同一层目录上的不同节点，不是两套归因系统。
