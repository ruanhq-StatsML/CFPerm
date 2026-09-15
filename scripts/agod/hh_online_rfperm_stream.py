#!/usr/bin/env python3
"""HH-RLHF → tidy stream table → OnlineRFPerm 连续实时检测 prototype.

你要的数据 form（先整理，再检测）::

    t_idx | batch | question | answer | y | feature[d]

  - t_idx   : 时间序（流里第几条；可用真实 timestamp 替换）
  - batch   : 实时窗（每 n_per 条一个 batch = 一次探针 hop）
  - question: Human prompt
  - answer  : Assistant reply（chosen / rejected）
  - y       : 偏好标签（chosen=1, rejected=0；cut 后可模拟映射 hop）
  - feature : 向量（默认 style+hash；可换成你们 embedding）

包（就 call 这个）::

    from agod.online_rfperm import (
        fit_online_probe, probe_err, hop_fires, po_risk0_rows, error_floor, shift_ratio
    )
    from agod.po_iptw import po_iptw_weights

数据::

    data/hf_cache/hh_rlhf_helpful_base_2500.jsonl   # Anthropic/hh-rlhf helpful-base

Usage::

    PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py
    PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py --featurizer style_hash
    PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py --featurizer style_only
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
OUT = ROOT / "results" / "agod" / "hh_online_stream"
DOCS = ROOT / "docs" / "biz"

from agod.online_rfperm import (  # noqa: E402
    error_floor,
    fit_online_probe,
    hop_fires,
    po_risk0_rows,
    probe_err,
    shift_ratio,
)
from agod.po_iptw import po_iptw_weights  # noqa: E402


# =============================================================================
# 1) 数据 form：HH → time-index / feature / question / answer
# =============================================================================

_HEDGE = ("maybe", "perhaps", "i think", "not sure", "probably", "might")
_FORMAL = ("therefore", "however", "furthermore", "regarding", "consequently")


def ensure_hh(n: int = 2500) -> Path:
    path = CACHE / "hh_rlhf_helpful_base_2500.jsonl"
    if path.exists():
        return path
    CACHE.mkdir(parents=True, exist_ok=True)
    from datasets import load_dataset

    ds = load_dataset(
        "Anthropic/hh-rlhf", data_dir="helpful-base", split="train", streaming=True
    )
    with path.open("w") as f:
        for i, row in enumerate(ds):
            f.write(
                json.dumps(
                    {"chosen": row["chosen"], "rejected": row["rejected"]},
                    ensure_ascii=False,
                )
                + "\n"
            )
            if i + 1 >= n:
                break
    return path


def human_prompt(dialog: str) -> str:
    m = re.search(r"Human:\s*(.*?)(?:\n\nAssistant:|$)", dialog, flags=re.S)
    return (m.group(1).strip() if m else dialog.strip())[:800]


def assistant_reply(dialog: str) -> str:
    parts = re.split(r"\n\nAssistant:", dialog)
    return parts[-1].strip() if len(parts) >= 2 else dialog.strip()


def style_feature(text: str) -> np.ndarray:
    """可解释文风 / register 特征（不依赖 LLM embedding）。

    Intuition：对齐漂移经常先表现为「说话语气」变了；这 10 维便宜、稳定、
    可审计。不够语义时再叠 embedding。
    """
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
        dtype=float,
    )


def hash_feature(texts: list[str], n_features: int = 128) -> np.ndarray:
    """HashingVectorizer — 仅作 **无 GPU / 无模型** 的语义弱代理。

    Intuition（为什么还留着）：
      - 碰撞哈希把 n-gram 投到固定维，streaming 友好、不存词表
      - 对「字面偏好翻转」这种 synthetic hop 够敏感，能把 OnlineRFPerm 跑通
    为什么不太靠谱：
      - 无语义几何；同义改写几乎随机
      - 碰撞噪声大；不适合当最终生产特征
    生产替换：sentence-transformers / 你们 serving embedding / 裁判隐层
    """
    vec = HashingVectorizer(
        n_features=n_features,
        alternate_sign=False,
        norm="l2",
        ngram_range=(1, 2),
    )
    return vec.transform(texts).toarray().astype(float)


def build_tidy_table(
    rows: list[dict],
    *,
    n_per: int = 80,
    cut_batch: int = 4,
    featurizer: str = "style_hash",
    hash_dim: int = 128,
    seed: int = 0,
) -> pd.DataFrame:
    """整理成你要的 form，再切 batch 做实时窗。"""
    rng = np.random.default_rng(seed)
    questions, answers, labels, styles = [], [], [], []
    for r in rows:
        q = human_prompt(r["chosen"])
        for key, lab in (("chosen", 1), ("rejected", 0)):
            a = assistant_reply(r[key])
            questions.append(q)
            answers.append(a)
            labels.append(lab)
            styles.append(style_feature(a))

    # 打乱对，避免 chosen/rejected 交错伪时间结构
    perm = rng.permutation(len(labels))
    questions = [questions[i] for i in perm]
    answers = [answers[i] for i in perm]
    labels = np.asarray(labels, dtype=int)[perm]
    styles = np.vstack(styles)[perm]

    # 截成整 batch
    n_use = (len(labels) // n_per) * n_per
    questions, answers = questions[:n_use], answers[:n_use]
    labels, styles = labels[:n_use], styles[:n_use]
    t_idx = np.arange(n_use, dtype=int)
    batch = np.repeat(np.arange(n_use // n_per), n_per)

    # 特征
    if featurizer == "style_only":
        feat = styles
        feat_name = "style10"
    elif featurizer == "hash_only":
        feat = hash_feature([f"{q}\n{a}" for q, a in zip(questions, answers)], hash_dim)
        feat_name = f"hash{hash_dim}"
    else:  # style_hash — demo 默认：可解释 + 弱语义
        h = hash_feature([f"{q}\n{a}" for q, a in zip(questions, answers)], hash_dim)
        feat = np.hstack([styles, h])
        feat_name = f"style10+hash{hash_dim}"

    # 模拟「偏好映射 hop」：cut 后 92% 标签翻转（DPO 连更失败模式）
    # 真实线上：不要人工翻转；y 用人工偏好 / 裁判分 / 纠错标
    y = labels.copy()
    after = batch >= cut_batch
    flip = after & (rng.random(n_use) < 0.92)
    y[flip] = 1 - y[flip]

    df = pd.DataFrame(
        {
            "t_idx": t_idx,
            "batch": batch,
            "question": questions,
            "answer": answers,
            "y": y,
            "y_raw_pref": labels,  # 翻转前的原始 chosen/rejected
            "flipped": flip.astype(int),
            "featurizer": feat_name,
        }
    )
    # feature 存成 list 方便 parquet；矩阵另存 npy
    df["feature"] = [feat[i].tolist() for i in range(n_use)]
    meta = {
        "cut_batch": int(cut_batch),
        "n_per": int(n_per),
        "featurizer": feat_name,
        "feature_dim": int(feat.shape[1]),
    }
    return df, feat, meta


# =============================================================================
# 2) 连续实时检测：按 batch 走 OnlineRFPerm（不是离线一次 AUROC）
# =============================================================================


def continuous_online_rfperm(
    feat: np.ndarray,
    y: np.ndarray,
    batch: np.ndarray,
    t_idx: np.ndarray,
    *,
    gate: float = 1.25,
    topk: int = 10,
    seed: int = 0,
) -> list[dict]:
    """真正的「连续实时」循环：每个新 batch 到齐 → 探针 → 火？→ PO → 动作。

    时钟：
      batch t-1 = T=0 探针训练集
      batch t   = 当前窗（算 e_now / po_risk0）
      e_prev    = 上一跳的 e_now（连续 OOS 比）

    quiet → w=1
    fire  → w=√po_risk0；审计 Top-k
    """
    feat = np.asarray(feat, dtype=float)
    y = np.asarray(y, dtype=int).ravel()
    batch = np.asarray(batch, dtype=int).ravel()
    t_idx = np.asarray(t_idx, dtype=int).ravel()
    k = int(batch.max()) + 1

    history: list[dict] = []
    e_prev: float | None = None

    for t in range(1, k):
        prev_m = batch == (t - 1)
        cur_m = batch == t
        if not prev_m.any() or not cur_m.any():
            continue

        X0, y0 = feat[prev_m], y[prev_m]
        X1, y1 = feat[cur_m], y[cur_m]
        t_idx_cur = t_idx[cur_m]

        # --- call agod ---
        probe = fit_online_probe(X0, y0, seed=seed + t, task="acc")
        e_now = probe_err(probe, X1, y1, task="acc")
        e_fl = error_floor("acc", int(cur_m.sum()))
        fired = hop_fires(e_now, e_prev, gate=gate, e_floor=e_fl)
        ratio = 1.0 if e_prev is None else shift_ratio(e_now, e_prev, e_floor=e_fl)
        po = po_risk0_rows(probe, X1, y1, task="acc")
        w = po_iptw_weights(po, mode="sqrt") if fired else np.ones_like(po)

        order = np.argsort(-po)
        top = order[: min(topk, len(order))]
        audit_t_idx = t_idx_cur[top].tolist()
        p_at_k = float(np.mean(y1[top])) if len(top) else float("nan")

        # 评估 / 刻画（下一步）
        rec = {
            "batch_t": int(t),
            "t_idx_start": int(t_idx_cur.min()),
            "t_idx_end": int(t_idx_cur.max()),
            "n": int(cur_m.sum()),
            "fired": bool(fired),
            "ratio": float(ratio),
            "e_now": float(e_now),
            "e_prev": None if e_prev is None else float(e_prev),
            "pos_rate": float(np.mean(y1)),
            "mean_po": float(np.mean(po)),
            "mean_w": float(np.mean(w)),
            "precision_at_k": p_at_k,
            "audit_t_idx_topk": audit_t_idx,
            "action": (
                "REJECT_or_LIMIT_merge + audit_topk(po_risk0) + optional_sqrtPO"
                if fired
                else "ACCEPT_merge + w=1"
            ),
        }
        history.append(rec)
        e_prev = e_now  # 连续：下一跳用本跳误差

    return history


# =============================================================================
# 3) 报告
# =============================================================================


def render_md(
    df: pd.DataFrame,
    hist: list[dict],
    meta: dict,
    *,
    data_path: str,
) -> str:
    cut = int(meta["cut_batch"])
    fires = [h for h in hist if h["fired"]]
    first = fires[0] if fires else None
    at_cut = next((h for h in hist if h["batch_t"] == cut), hist[-1] if hist else {})

    hist_rows = "\n".join(
        f"| {h['batch_t']} | {h['t_idx_start']}–{h['t_idx_end']} | "
        f"{h['fired']} | {h['ratio']:.3f} | {h['e_now']:.3f} | "
        f"{h['precision_at_k']:.2f} | `{h['action'][:28]}…` |"
        for h in hist
    )

    return f"""# HH-RLHF → tidy 表 → OnlineRFPerm 连续实时检测

## 数据是什么

| 项 | 值 |
|----|-----|
| 源 | **Anthropic/hh-rlhf** `helpful-base`（对齐 / 偏好） |
| 缓存 | `{data_path}` |
| 整理后 | `results/agod/hh_online_stream/stream_table.parquet` |
| 行数 | {len(df)}（chosen+rejected 展开） |
| featurizer | `{meta.get('featurizer')}` |
| n_per / cut_batch | {meta.get('n_per')} / {cut} |

### 你要的 form（样例）

| t_idx | batch | question | answer | y | feature |
|------:|------:|----------|--------|--:|---------|
| {df.iloc[0]['t_idx']} | {df.iloc[0]['batch']} | {(df.iloc[0]['question'] or '')[:40]}… | {(df.iloc[0]['answer'] or '')[:40]}… | {df.iloc[0]['y']} | dim={len(df.iloc[0]['feature'])} |
| {df.iloc[1]['t_idx']} | {df.iloc[1]['batch']} | {(df.iloc[1]['question'] or '')[:40]}… | {(df.iloc[1]['answer'] or '')[:40]}… | {df.iloc[1]['y']} | … |

真实线上：`t_idx` 换成 serving / 标注时间戳；`y` 换成人工偏好或裁判分；**不要**做 demo 里的标签翻转。

## 怎么连续实时检测（核心）

```text
新 batch t 到齐
    → μ0 = fit_online_probe(batch t-1)     # 上一窗当制度探针
    → e_now = err(μ0, batch t)
    → fire ⇔ e_now / e_prev ≥ γ 且 e_prev≥floor
    → po_i  = po_risk0_rows(μ0, batch t)
    → fire? w=√po , audit Top-k(po)
      quiet? w=1 退火
    → e_prev ← e_now                       # 连续：下一跳接着比
```

包调用（就这几行）::

```python
from agod.online_rfperm import fit_online_probe, probe_err, hop_fires, po_risk0_rows
from agod.po_iptw import po_iptw_weights

probe = fit_online_probe(X_prev, y_prev, task="acc")
e_now = probe_err(probe, X_cur, y_cur, task="acc")
fire  = hop_fires(e_now, e_prev, gate=1.25)
po    = po_risk0_rows(probe, X_cur, y_cur, task="acc")
w     = po_iptw_weights(po, mode="sqrt") if fire else np.ones_like(po)
```

## 本跑 hop 表（连续）

| batch | t_idx | fired | ratio | e_now | P@k | action |
|------:|-------|:-----:|------:|------:|----:|--------|
{hist_rows}

- 首次 fire：batch=`{(first or {}).get('batch_t')}`，ratio≈`{(first or {}).get('ratio')}`
- cut={cut} 处：fired=`{at_cut.get('fired')}`，ratio≈`{at_cut.get('ratio')}`，P@k=`{at_cut.get('precision_at_k')}`
- 总 fires：{len(fires)} / {len(hist)} hops

## HashingVectorizer 靠谱吗？Intuition + 替代

| 方案 | Intuition | 何时用 | 风险 |
|------|-----------|--------|------|
| **style 10 维** | 语气/长度/正式度 = P(X) 画像，便宜可审计 | 文风轴、快速探针 | 不解语义偏好 |
| **HashingVectorizer** | 固定维 n-gram 投影，streaming、无词表 | **仅 demo / 无 embedding 时** | 无语义、碰撞噪声；**生产不推荐当唯一特征** |
| **style + hash**（本 demo 默认） | 可解释轴 + 弱字面信号，保证 hop 可复现 | 把 OnlineRFPerm 闭环跑通 | 仍非语义 |
| **sentence-transformers / 业务 embedding** | 语义几何接近裁判/模型内部表征 | **生产默认** | 要模型与版本钉扎 |
| **LLM hidden / reward head** | 与对齐目标同空间 | 已有 reward model 时 | 贵；要防泄漏 |

**结论：** HashingVectorizer **不够靠谱当最终特征**；它只是「没有 GPU/没有 embedding 服务时，仍能把时间序 + OnlineRFPerm 闭环演示出来」的弱代理。你接生产时把 `feature` 列换成你们的 embedding 即可，**检测逻辑一行都不用改**。

## 下一步评估 / 刻画

1. **制度**：fire 率、ratio 分布、cut 对齐（本表）
2. **排序**：po_risk0 P@k / AUROC（相对 y）
3. **动作**：fire→拒合并/审计；quiet→放行
4. **接到¥**：同一 fire 窗 → 客服 act-vs-ignore 账（另见 biz SQL）

## 怎么跑

```bash
PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py
PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py --featurizer style_only
```
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--featurizer",
        choices=["style_hash", "style_only", "hash_only"],
        default="style_hash",
    )
    ap.add_argument("--n-pairs", type=int, default=1000)
    ap.add_argument("--n-per", type=int, default=80)
    ap.add_argument("--cut-batch", type=int, default=4)
    ap.add_argument("--gate", type=float, default=1.20)
    ap.add_argument("--topk", type=int, default=10)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    path = ensure_hh()
    rows = []
    with path.open() as f:
        for i, line in enumerate(f):
            rows.append(json.loads(line))
            if i + 1 >= args.n_pairs:
                break

    df, feat, meta = build_tidy_table(
        rows,
        n_per=args.n_per,
        cut_batch=args.cut_batch,
        featurizer=args.featurizer,
    )

    # 落盘：tidy 表 + 特征矩阵
    out_parquet = OUT / "stream_table.parquet"
    df.to_parquet(out_parquet, index=False)
    np.save(OUT / "feature_matrix.npy", feat)
    df.drop(columns=["feature"]).head(40).to_csv(
        OUT / "stream_table_preview.csv", index=False
    )
    (OUT / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2))

    hist = continuous_online_rfperm(
        feat,
        df["y"].to_numpy(),
        df["batch"].to_numpy(),
        df["t_idx"].to_numpy(),
        gate=args.gate,
        topk=args.topk,
    )
    (OUT / "hop_history.json").write_text(
        json.dumps(hist, ensure_ascii=False, indent=2)
    )

    summary = {
        "data": str(path.relative_to(ROOT)),
        "featurizer": meta["featurizer"],
        "n_rows": len(df),
        "n_batches": int(df["batch"].max()) + 1,
        "n_per": int(meta["n_per"]),
        "cut_batch": int(meta["cut_batch"]),
        "gate": args.gate,
        "n_fires": int(sum(1 for h in hist if h["fired"])),
        "first_fire_batch": next((h["batch_t"] for h in hist if h["fired"]), None),
        "columns": ["t_idx", "batch", "question", "answer", "y", "feature"],
        "package": "agod.online_rfperm + agod.po_iptw",
        "artifacts": {
            "parquet": str(out_parquet.relative_to(ROOT)),
            "feature_npy": "results/agod/hh_online_stream/feature_matrix.npy",
            "hop_history": "results/agod/hh_online_stream/hop_history.json",
        },
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2)
    )

    md = render_md(df, hist, meta, data_path=str(path.relative_to(ROOT)))
    (OUT / "HH_ONLINE_RFPERM_STREAM.md").write_text(md)
    (DOCS / "HH_ONLINE_RFPERM_STREAM.md").write_text(md)
    print(md)
    print(f"wrote {out_parquet}")
    print(f"wrote {DOCS / 'HH_ONLINE_RFPERM_STREAM.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
