# Graph localization（已标超纲 / deferred）

> **决定：本阶段不做。**  
> ULS / PPR-Nibble / GraphScan 听起来可用，但仍算超纲——主协议只保留  
> **图谱特征 → \(X\)** + **cmean \(\delta/D/r/R\)** 缩支撑  
> （[`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md)、
> [`Feature_Entity_Joint_CMean.tex`](./Feature_Entity_Joint_CMean.tex)）。

## 明确不做

| 项 | 状态 |
|---|---|
| Community detection | ❌ 超 scope |
| Ego-centric network | ❌ 超 scope |
| ULS / PPR-Nibble / GraphScan / conductance 收核 | ❌ **本阶段也超纲**（deferred） |
| GNN 归因 | ❌ |

若以后单开实验，可复活分数驱动局部化；**默认闭环不依赖图上收核。**

**Leiden 上游 key（审阅稿，非默认）：**  
见 [`LEIDEN_FOR_DRILL_ELABORATION.md`](./LEIDEN_FOR_DRILL_ELABORATION.md)——社区只作 `drill_key` 候选，**不**用模块度当 Drill 门。

## 默认闭环（仅此）

```text
图谱特征 → X
标准化 → cmean 列 δ / 行 r → 早停得 K* → 一次 FSDS
（多 entity key / F→E 联合按主文档；非因果、非 subgroup）
```
