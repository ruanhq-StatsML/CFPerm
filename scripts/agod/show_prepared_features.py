#!/usr/bin/env python3
"""展示已做好的推荐/对齐特征：在哪、长什么样、怎么读。

特征产物（已跑好）::

  results/agod/hh_online_stream/stream_table.parquet   # tidy: t_idx/q/a/y/feature
  results/agod/hh_online_stream/feature_matrix.npy     # (N, 138) = style10 + hash128
  results/agod/hh_online_stream/stream_table_preview.csv
  results/agod/hh_online_stream/meta.json

原始缓存::

  data/hf_cache/hh_rlhf_helpful_base_2500.jsonl

画风/文风 10 维名字（feature[0:10]）::

  n_tok, n_char, avg_word, qmark, bang, hedge, formal, i_count, newlines, upper_ratio

Usage::

  PYTHONPATH=. python3 scripts/agod/show_prepared_features.py
  PYTHONPATH=. python3 scripts/agod/show_prepared_features.py --n 3
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "results" / "agod" / "hh_online_stream"
DOCS = ROOT / "docs" / "biz"
OUT = ART

STYLE_DIMS = [
    ("n_tok", "词数/200"),
    ("n_char", "字符数/800"),
    ("avg_word", "平均词长/10"),
    ("qmark", "? 密度"),
    ("bang", "! 密度"),
    ("hedge", "maybe/perhaps/... 犹豫词"),
    ("formal", "therefore/however/... 正式词 ← 画风主轴之一"),
    ("i_count", "第一人称 i"),
    ("newlines", "换行密度"),
    ("upper_ratio", "大写比例"),
]


def load_artifacts():
    parquet = ART / "stream_table.parquet"
    npy = ART / "feature_matrix.npy"
    if not parquet.exists() or not npy.exists():
        raise SystemExit(
            f"missing artifacts under {ART}\n"
            "run: PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py"
        )
    df = pd.read_parquet(parquet)
    feat = np.load(npy)
    meta = json.loads((ART / "meta.json").read_text()) if (ART / "meta.json").exists() else {}
    summary = (
        json.loads((ART / "summary.json").read_text())
        if (ART / "summary.json").exists()
        else {}
    )
    return df, feat, meta, summary


def describe_row(df: pd.DataFrame, feat: np.ndarray, i: int) -> dict:
    r = df.iloc[i]
    x = feat[i]
    style = {name: float(x[j]) for j, (name, _) in enumerate(STYLE_DIMS)}
    return {
        "t_idx": int(r["t_idx"]),
        "batch": int(r["batch"]),
        "y": int(r["y"]),
        "y_raw_pref": int(r["y_raw_pref"]),
        "flipped": int(r["flipped"]),
        "question": str(r["question"])[:160],
        "answer": str(r["answer"])[:200],
        "featurizer": str(r["featurizer"]),
        "feature_dim": int(len(x)),
        "style10": style,
        "hash128_l2": float(np.linalg.norm(x[10:])),
        "hash128_head8": [float(v) for v in x[10:18]],
    }


def render_md(examples: list[dict], meta: dict, summary: dict) -> str:
    ex_blocks = []
    for e in examples:
        style_lines = "\n".join(
            f"| `{k}` | {STYLE_DIMS[i][1]} | {v:.4f} |"
            for i, (k, v) in enumerate(e["style10"].items())
        )
        ex_blocks.append(
            f"""### 样例 t_idx={e['t_idx']}（batch={e['batch']}, y={e['y']}）

- **question:** {e['question']}
- **answer:** {e['answer'][:180]}…
- featurizer=`{e['featurizer']}` · dim={e['feature_dim']} · flipped={e['flipped']}

| style 维 | 含义 | 值 |
|----------|------|-----|
{style_lines}

hash128 ‖·‖₂ ≈ {e['hash128_l2']:.3f}；head8 = `{e['hash128_head8']}`
"""
        )

    return f"""# 已做好的特征在哪 + 样例 Prototype

## 路径

| 文件 | 内容 |
|------|------|
| `results/agod/hh_online_stream/stream_table.parquet` | tidy 表：`t_idx, batch, question, answer, y, feature` |
| `results/agod/hh_online_stream/feature_matrix.npy` | 矩阵 **({summary.get('n_rows', '?')}, 138)** = style10 ⊕ hash128 |
| `results/agod/hh_online_stream/stream_table_preview.csv` | 无向量的人眼预览 |
| `results/agod/hh_online_stream/meta.json` | featurizer / cut / n_per |
| `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` | 原始 HH-RLHF 缓存 |

meta：`{json.dumps(meta, ensure_ascii=False)}`

## 特征怎么拼的

```text
HH dialog
  → human_prompt(chosen)  → question
  → assistant_reply(chosen/rejected) → answer
  → style_feature(answer) → 10 维画风/文风
  → HashingVectorizer(question+"\\n"+answer) → 128 维弱语义
  → feature = concat(style10, hash128)     # 共 138
```

**画风检测 intuition：** `formal` / `hedge` / 长度 / 标点 描述的是 \(P(X)\) 说话画像，不是偏好对错。  
高 `style_domain_auc`（见 `results/agod/hf_landing/hh_judge_style.json`）说明文风可分；偏好 hop 另走 OnlineRFPerm + `y`。

## 样例

{chr(10).join(ex_blocks)}

## 怎么读（三行 Python）

```python
import pandas as pd, numpy as np
df = pd.read_parquet("results/agod/hh_online_stream/stream_table.parquet")
X  = np.load("results/agod/hh_online_stream/feature_matrix.npy")  # (2000, 138)
# X[i, :10] = style；X[i, 10:] = hash
```

生产：只用 parquet 的 `question/answer/y`，把 `feature` 换成你们 embedding 即可。

## 还有其他相关产物吗？

| 产物 | 路径 | 用途 |
|------|------|------|
| HH Judge + style_auc | `results/agod/hf_landing/hh_judge_style.json` | 文风 AUC / judge_err_ratio |
| Halu 推理 hop | `results/agod/hf_landing/halu_regime_rag.json` | 幻觉制度 knobs |
| 连续检测 hop 表 | `results/agod/hh_online_stream/hop_history.json` | fire / po Top-k |
| 生成脚本 | `scripts/agod/hh_online_rfperm_stream.py` | data-manipulation |

```bash
PYTHONPATH=. python3 scripts/agod/show_prepared_features.py
```
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2, help="how many example rows to show")
    ap.add_argument(
        "--which",
        choices=["head", "formal", "hedge", "flipped"],
        default="head",
        help="pick examples: head / high-formal / high-hedge / flipped label",
    )
    args = ap.parse_args()

    df, feat, meta, summary = load_artifacts()
    assert feat.shape[0] == len(df)

    if args.which == "head":
        idxs = list(range(min(args.n, len(df))))
    elif args.which == "formal":
        idxs = list(np.argsort(-feat[:, 6])[: args.n])  # formal dim
    elif args.which == "hedge":
        idxs = list(np.argsort(-feat[:, 5])[: args.n])
    else:
        flipped = np.flatnonzero(df["flipped"].to_numpy() == 1)
        idxs = flipped[: args.n].tolist() if len(flipped) else [0]

    examples = [describe_row(df, feat, i) for i in idxs]
    (OUT / "feature_examples.json").write_text(
        json.dumps(
            {
                "paths": {
                    "parquet": str((ART / "stream_table.parquet").relative_to(ROOT)),
                    "feature_npy": str((ART / "feature_matrix.npy").relative_to(ROOT)),
                },
                "shape": list(feat.shape),
                "style_dims": [{"name": n, "meaning": m} for n, m in STYLE_DIMS],
                "examples": examples,
            },
            ensure_ascii=False,
            indent=2,
        )
    )

    md = render_md(examples, meta, summary)
    (OUT / "PREPARED_FEATURES.md").write_text(md)
    (DOCS / "PREPARED_FEATURES.md").write_text(md)
    print(md)
    print(f"wrote {DOCS / 'PREPARED_FEATURES.md'}")
    print(f"wrote {OUT / 'feature_examples.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
