# 多模态归因（CFPerm）：permute \(T\)，报告的是模态块

组别归因问哪些 **坐标** \(X_j\) 带着组间差。
这一件同一套 permute-\(T\) 循环，单位换成 **模态块**：

\[
X = [\text{text}\mid \text{vision/dense}\mid \text{graph}\mid \text{fusion}]
\]

\(T\) 仍是组别 / 包版本。打乱的是 \(T\)，不是 \(X\)，不是像素，不是原文。
\(I_m=\sum_{j\in m}\mathrm{importance}_j\)。块过跨块阈值才算这个模态带着差。

## 逻辑

```
1. 观察当前包，写下这一跳的 Y
2. X 按模态分块（文本几何 / 视觉或 dense / 图 / 融合）
3. T 是包版本或组别，不是时间 batch 本身
4. permute T → 特征 importance，再在块内求和
5. 跨块阈值；命中的是模态，不是单列 x_*
```

和组别归因、`R/CFPerm_vimp.R` 同一套决定。只是最后报告的粒度是块。

## 盘上的流

| 流 | \(T\) | 块 | \(Y\) |
|---|---|---|---|
| 混合检索 dense hop | 切点前 / 切点后的 dense 包 | query / sparse / dense / fusion | 融合答案成不成立 |
| Graph-RAG community | local 包 / community 包 | query-seed / graph 拓扑 | 支撑节点在包里 |

```bash
PYTHONPATH=. python3 scripts/prototype_multimodal_attribution.py
PYTHONPATH=. python3 scripts/prototype_group_attribution.py
```

读数：`results/manuscript/multimodal_attribution/REPORT.md`。

Planted：三模态里只有 vision 与 \(T\) 交互，top 应是 `vision`。Null：\(Y\) 靠 text/graph、不靠 \(T\)，块应 quiet。

## 两件归因怎么分开

| 对象 | 单位 | 动哪一列 |
|---|---|---|
| 在线 serving 闸 | 时间窗 | `batch` |
| 多组别归因 | 坐标 \(X_j\) | 组别 \(T\) |
| 多模态归因 | 模态块 | 组别 / 包版本 \(T\) |

不要把块命中读成「所以时间闸是这路引起的」。HH chosen、像素、原文都不是 \(Y\) 也不是 \(X\)。
