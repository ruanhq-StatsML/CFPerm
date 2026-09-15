#!/usr/bin/env python3
"""LLM 推理变点检测 prototype：两数据集 → 连续时间流 → OnlineRFPerm。

你要的数据 form（一行一个时刻）::

    {
      "t": int,                    # 连续时间下标
      "question": str,
      "answer": str,
      "embedding": np.ndarray,     # shape (p,) 或 (p, 1)
      "score": float,              # 监督/弱标签：幻觉=0/1 或 偏好=0/1
    }

两个 dataset（本地 cache）::

  1) HaluEval QA  — 推理幻觉制度变点（主）
  2) HH-RLHF      — 对齐/偏好映射变点（对照）

方法::

  变点 / 制度跳变  →  OnlineRFPerm（本脚本）
  实例 PO-risk 排序 →  你们实验已论证用 BOCPD 更好（这里不硬塞 OnlineRFPerm 做 PO）

Usage::

  PYTHONPATH=. python3 scripts/agod/llm_stream_changepoint.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import HashingVectorizer

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache"
OUT = ROOT / "results" / "agod" / "llm_changepoint"
DOCS = ROOT / "docs" / "biz"

from agod.online_rfperm import (  # noqa: E402
    error_floor,
    fit_online_probe,
    hop_fires,
    probe_err,
    shift_ratio,
)


# =============================================================================
# 1) 最小解析
# =============================================================================


def human_prompt(dialog: str) -> str:
    m = re.search(r"Human:\s*(.*?)(?:\n\nAssistant:|$)", dialog, flags=re.S)
    return (m.group(1).strip() if m else dialog.strip())[:800]


def assistant_reply(dialog: str) -> str:
    parts = re.split(r"\n\nAssistant:", dialog)
    return parts[-1].strip() if len(parts) >= 2 else dialog.strip()


_HEDGE = ("maybe", "perhaps", "i think", "not sure", "probably", "might")
_FORMAL = ("therefore", "however", "furthermore", "regarding", "consequently")


def style_feature(text: str) -> np.ndarray:
    """10-d register features — cheap proxy until real embeddings land."""
    t = (text or "").strip()
    low = t.lower()
    toks = re.findall(r"[a-zA-Z']+", low)
    chars = max(len(t), 1)
    avg_w = float(np.mean([len(w) for w in toks])) if toks else 0.0
    return np.asarray(
        [
            len(toks) / 200.0,
            chars / 800.0,
            avg_w / 10.0,
            low.count("?") / 5.0,
            low.count("!") / 5.0,
            sum(low.count(h) for h in _HEDGE) / 5.0,
            sum(low.count(h) for h in _FORMAL) / 5.0,
            low.count(" i ") / 10.0,
            t.count("\n") / 10.0,
            (sum(c.isupper() for c in t) / chars),
        ],
        dtype=np.float64,
    )


def embed_qa(
    questions: list[str],
    answers: list[str],
    hash_dim: int = 64,
    *,
    pref_signal: np.ndarray | None = None,
) -> np.ndarray:
    """Demo embedding = style10 + hash(q\\n a) [+ optional pref channel].

    生产换成业务 embedding，form 不变。pref_signal 仅用于 HH 对照里
    制造「偏好映射可学 → cut 后翻转」的清晰变点。
    """
    styles = np.vstack([style_feature(a) for a in answers])
    vec = HashingVectorizer(
        n_features=hash_dim, alternate_sign=False, norm="l2", ngram_range=(1, 2)
    )
    texts = [f"{q}\n{a}" for q, a in zip(questions, answers)]
    h = vec.transform(texts).toarray().astype(np.float64)
    parts = [styles, h]
    if pref_signal is not None:
        parts.append(np.asarray(pref_signal, dtype=np.float64).reshape(-1, 1))
    return np.hstack(parts)


def load_jsonl(path: Path, n: int) -> list[dict]:
    rows = []
    with path.open() as f:
        for line in f:
            rows.append(json.loads(line))
            if len(rows) >= n:
                break
    return rows


# =============================================================================
# 2) Data manipulation → list[{t, question, answer, embedding, score}]
# =============================================================================


def stream_from_halueval(
    rows: list[dict],
    *,
    emb_dim: int = 64,
    n_per: int = 100,
    cut_batch: int = 4,
    seed: int = 0,
) -> list[dict]:
    """推理面：question/answer + score=幻觉标签；cut 后幻觉率抬高制造变点。"""
    _ = seed  # deterministic label remapping; no RNG needed
    questions, answers, scores = [], [], []
    for r in rows:
        q = str(r.get("question") or "")
        a = str(r.get("answer") or "")
        hall = str(r.get("hallucination", "")).strip().lower() in ("yes", "true", "1")
        questions.append(q)
        answers.append(a)
        scores.append(1.0 if hall else 0.0)

    emb = embed_qa(questions, answers, emb_dim)
    n_use = (len(scores) // n_per) * n_per
    questions, answers = questions[:n_use], answers[:n_use]
    scores = np.asarray(scores[:n_use], dtype=np.float64)
    emb = emb[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)

    # 制度变点：cut 前压幻觉、cut 后抬幻觉（模拟推理质量跳变）
    for b in range(int(batch.max()) + 1):
        m = batch == b
        idx = np.flatnonzero(m)
        if b < cut_batch:
            pos = idx[scores[idx] > 0.5]
            scores[pos[max(2, len(pos) // 15) :]] = 0.0
        else:
            neg = idx[scores[idx] < 0.5]
            scores[neg[: int(0.85 * len(neg))]] = 1.0

    stream = []
    for t in range(n_use):
        stream.append(
            {
                "t": int(t),
                "question": questions[t],
                "answer": answers[t],
                "embedding": emb[t].reshape(-1, 1),  # (p, 1)
                "score": float(scores[t]),
                "batch": int(batch[t]),
                "dataset": "halueval_qa",
            }
        )
    return stream


def stream_from_hh(
    rows: list[dict],
    *,
    emb_dim: int = 64,
    n_per: int = 80,
    cut_batch: int = 4,
    seed: int = 1,
) -> list[dict]:
    """对齐面：question/answer + score=偏好；cut 后翻转偏好映射。"""
    rng = np.random.default_rng(seed)
    questions, answers, scores = [], [], []
    for r in rows:
        q = human_prompt(r["chosen"])
        for key, lab in (("chosen", 1.0), ("rejected", 0.0)):
            questions.append(q)
            answers.append(assistant_reply(r[key]))
            scores.append(lab)

    perm = rng.permutation(len(scores))
    questions = [questions[i] for i in perm]
    answers = [answers[i] for i in perm]
    scores = np.asarray(scores, dtype=np.float64)[perm]
    # pref channel = 原始偏好；cut 后只翻 y，不翻 channel → P(Y|X) 在 cut 断裂
    emb = embed_qa(questions, answers, emb_dim, pref_signal=scores)

    n_use = (len(scores) // n_per) * n_per
    questions, answers = questions[:n_use], answers[:n_use]
    scores, emb = scores[:n_use], emb[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)

    # 偏好映射 hop：cut 后翻转标签（embedding 里仍留着旧映射信号）
    after = batch >= cut_batch
    flip = after & (rng.random(n_use) < 0.92)
    scores = scores.copy()
    scores[flip] = 1.0 - scores[flip]

    stream = []
    for t in range(n_use):
        stream.append(
            {
                "t": int(t),
                "question": questions[t],
                "answer": answers[t],
                "embedding": emb[t].reshape(-1, 1),
                "score": float(scores[t]),
                "batch": int(batch[t]),
                "dataset": "hh_rlhf",
            }
        )
    return stream


def stream_to_matrices(stream: list[dict]):
    """list[{t,q,a,emb,score}] → X (n,p), y (n,), batch (n,), t (n,)"""
    X = np.stack([r["embedding"].reshape(-1) for r in stream], axis=0)
    y = np.asarray([r["score"] for r in stream], dtype=np.float64)
    batch = np.asarray([r["batch"] for r in stream], dtype=int)
    t = np.asarray([r["t"] for r in stream], dtype=int)
    return X, y, batch, t


# =============================================================================
# 3) OnlineRFPerm 变点检测（连续 batch）
# =============================================================================


def online_rfperm_changepoints(
    stream: list[dict],
    *,
    gate: float = 1.20,
    seed: int = 0,
) -> list[dict]:
    """连续时间流上的 OnlineRFPerm 变点：batch t-1 → 探针，batch t → e_now / fire。"""
    X, y, batch, t = stream_to_matrices(stream)
    y_cls = (y >= 0.5).astype(int)
    k = int(batch.max()) + 1
    hops: list[dict] = []
    e_prev: float | None = None
    for bt in range(1, k):
        prev_m = batch == (bt - 1)
        cur_m = batch == bt
        if not prev_m.any() or not cur_m.any():
            continue
        probe = fit_online_probe(
            X[prev_m], y_cls[prev_m], seed=int(seed) + bt, task="acc"
        )
        e_now = float(probe_err(probe, X[cur_m], y_cls[cur_m], task="acc"))
        e_fl = error_floor("acc", int(cur_m.sum()))
        fired = hop_fires(e_now, e_prev, gate=gate, e_floor=e_fl)
        ratio = (
            1.0
            if e_prev is None
            else float(shift_ratio(e_now, e_prev, e_floor=e_fl))
        )
        hops.append(
            {
                "batch_t": int(bt),
                "t_start": int(t[cur_m].min()),
                "t_end": int(t[cur_m].max()),
                "fired": bool(fired),
                "ratio": ratio,
                "e_now": e_now,
                "e_prev": None if e_prev is None else float(e_prev),
                "pos_rate": float(np.mean(y_cls[cur_m])),
            }
        )
        e_prev = e_now
    return hops


def save_stream_parquet(stream: list[dict], path: Path) -> None:
    """落盘时 embedding 存 list；另存 npy。"""
    rows = []
    embs = []
    for r in stream:
        e = r["embedding"].reshape(-1)
        embs.append(e)
        rows.append(
            {
                "t": r["t"],
                "question": r["question"],
                "answer": r["answer"],
                "score": r["score"],
                "batch": r["batch"],
                "dataset": r["dataset"],
                "embedding": e.tolist(),
            }
        )
    pd.DataFrame(rows).to_parquet(path, index=False)
    np.save(path.with_suffix(".npy"), np.stack(embs, axis=0))


# =============================================================================
# 4) 报告
# =============================================================================


def render_md(results: dict) -> str:
    blocks = []
    for name, r in results.items():
        hops = r["hops"]
        fires = [h for h in hops if h["fired"]]
        hop_tbl = "\n".join(
            f"| {h['batch_t']} | {h['t_start']}–{h['t_end']} | {h['fired']} | "
            f"{h['ratio']:.3f} | {h['pos_rate']:.2f} |"
            for h in hops[:12]
        )
        ex = r["example"]
        blocks.append(
            f"""### {name}

| 项 | 值 |
|----|-----|
| 源 | `{r['source']}` |
| n / p / n_per | {r['n']} / {r['p']} / {r['n_per']} |
| cut_batch | {r['cut_batch']} |
| fires | {len(fires)} / {len(hops)} |
| first_fire_batch | {(fires[0]['batch_t'] if fires else None)} |

样例时刻 t={ex['t']}::

```json
{{
  "t": {ex['t']},
  "question": {json.dumps(ex['question'][:80], ensure_ascii=False)},
  "answer": {json.dumps(ex['answer'][:80], ensure_ascii=False)},
  "embedding": "np.ndarray shape {ex['embedding_shape']}",
  "score": {ex['score']}
}}
```

| batch | t | fired | ratio | score均值 |
|------:|---|:-----:|------:|----------:|
{hop_tbl}
"""
        )

    return f"""# LLM 推理变点：两数据集 manipulation → OnlineRFPerm

## 数据 form（就这一列时间流）

```python
{{
  "t": int,
  "question": str,
  "answer": str,
  "embedding": np.ndarray,   # shape (p,) or (p, 1)
  "score": float,            # Halu=幻觉0/1；HH=偏好0/1
}}
```

## 方法怎么用

| 问题 | 方法 | 本原型 |
|------|------|--------|
| **大模型推理有没有制度/质量变点？** | **OnlineRFPerm** | ✅ 两数据集已跑 |
| 实例风险排序 / PO-risk | **BOCPD**（你们实验已论证更好） | 不在此硬塞 OnlineRFPerm |

```python
# stream: list[dict] with embedding + score
from agod.online_rfperm import fit_online_probe, probe_err, hop_fires, error_floor
# 每个新 batch：用上一窗拟合探针 → 本窗算 e_now → 与 e_prev 比 → fire?
# fired == 变点/制度火
```

> PO-risk 实例排序请用 **BOCPD**（实验已论证）；本脚本只做 OnlineRFPerm 变点。

## 本跑

{chr(10).join(blocks)}

## 产物

| 文件 | 内容 |
|------|------|
| `results/agod/llm_changepoint/halu_stream.parquet` | Halu 时间流 |
| `results/agod/llm_changepoint/halu_stream.npy` | embedding 矩阵 |
| `results/agod/llm_changepoint/hh_stream.parquet` | HH 时间流 |
| `results/agod/llm_changepoint/*_hops.json` | OnlineRFPerm fire 表 |

```bash
PYTHONPATH=. python3 scripts/agod/llm_stream_changepoint.py
```
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", type=float, default=1.20)
    ap.add_argument("--emb-dim", type=int, default=128)
    ap.add_argument("--halu-n", type=int, default=1200)
    ap.add_argument("--hh-n", type=int, default=1000)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    results = {}

    # --- dataset 1: HaluEval (inference changepoint) ---
    halu_path = CACHE / "halueval_qa_3000.jsonl"
    halu_rows = load_jsonl(halu_path, args.halu_n)
    halu_stream = stream_from_halueval(halu_rows, emb_dim=args.emb_dim, n_per=100, cut_batch=4)
    save_stream_parquet(halu_stream, OUT / "halu_stream.parquet")
    halu_hops = online_rfperm_changepoints(halu_stream, gate=args.gate, seed=0)
    (OUT / "halu_hops.json").write_text(json.dumps(halu_hops, indent=2))
    ex0 = halu_stream[0]
    results["HaluEval（推理幻觉变点）"] = {
        "source": str(halu_path.relative_to(ROOT)),
        "n": len(halu_stream),
        "p": int(ex0["embedding"].shape[0]),
        "n_per": 100,
        "cut_batch": 4,
        "hops": halu_hops,
        "example": {
            "t": ex0["t"],
            "question": ex0["question"],
            "answer": ex0["answer"],
            "score": ex0["score"],
            "embedding_shape": list(ex0["embedding"].shape),
        },
    }

    # --- dataset 2: HH-RLHF (preference map changepoint) ---
    hh_path = CACHE / "hh_rlhf_helpful_base_2500.jsonl"
    hh_rows = load_jsonl(hh_path, args.hh_n)
    hh_stream = stream_from_hh(hh_rows, emb_dim=args.emb_dim, n_per=80, cut_batch=4)
    save_stream_parquet(hh_stream, OUT / "hh_stream.parquet")
    hh_hops = online_rfperm_changepoints(hh_stream, gate=args.gate, seed=1)
    (OUT / "hh_hops.json").write_text(json.dumps(hh_hops, indent=2))
    ex1 = hh_stream[0]
    results["HH-RLHF（偏好映射变点）"] = {
        "source": str(hh_path.relative_to(ROOT)),
        "n": len(hh_stream),
        "p": int(ex1["embedding"].shape[0]),
        "n_per": 80,
        "cut_batch": 4,
        "hops": hh_hops,
        "example": {
            "t": ex1["t"],
            "question": ex1["question"],
            "answer": ex1["answer"],
            "score": ex1["score"],
            "embedding_shape": list(ex1["embedding"].shape),
        },
    }

    summary = {
        "form": ["t", "question", "answer", "embedding(p,1)", "score"],
        "method_changepoint": "OnlineRFPerm",
        "method_po_risk_recommended": "BOCPD (per prior experiments; not OnlineRFPerm)",
        "gate": args.gate,
        "halu_fires": sum(1 for h in halu_hops if h["fired"]),
        "hh_fires": sum(1 for h in hh_hops if h["fired"]),
        "halu_first_fire": next((h["batch_t"] for h in halu_hops if h["fired"]), None),
        "hh_first_fire": next((h["batch_t"] for h in hh_hops if h["fired"]), None),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    md = render_md(results)
    (OUT / "LLM_STREAM_CHANGEPOINT.md").write_text(md)
    (DOCS / "LLM_STREAM_CHANGEPOINT.md").write_text(md)
    print(md)
    print(f"wrote {OUT}")
    print(f"wrote {DOCS / 'LLM_STREAM_CHANGEPOINT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
