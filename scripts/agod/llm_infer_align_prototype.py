#!/usr/bin/env python3
"""LLM 推理 / 对齐 prototype — 直接 call ``agod`` package + HF 子集数据。

=============================================================================
你要 call 的包（本仓库）
=============================================================================
  from agod import (
      run_rfperm_stream,      # OnlineRFPerm 流：火不火
      fit_online_probe,       # 浅层 RF 制度探针 μ0
      po_risk0_rows,          # 实例 PO-risk（审谁）
  )
  from agod.online_rfperm import hop_fires, quantify_last_two, probe_err
  from agod.po_iptw import po_iptw_weights   # w = √po_risk0（门控重加权）
  from agod.po_refit import Stream

=============================================================================
数据在哪
=============================================================================
  data/hf_cache/hh_rlhf_helpful_base_2500.jsonl   # 对齐 / Judge / 文风
  data/hf_cache/halueval_qa_3000.jsonl            # 推理幻觉制度 / RAG
  （没有就自动从 HF 拉一小截缓存）

=============================================================================
两条落地路径
=============================================================================
  推理：HaluEval → OnlineRFPerm fire → po_risk0 Top-k 审计 / √PO 升配
  对齐：HH-RLHF  → RFPerm-as-Judge → fire 拒合并 / Top-k 复核 / √PO

Usage::

  PYTHONPATH=. python3 scripts/agod/llm_infer_align_prototype.py
  PYTHONPATH=. python3 scripts/agod/llm_infer_align_prototype.py --only infer
  PYTHONPATH=. python3 scripts/agod/llm_infer_align_prototype.py --only align
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.metrics import roc_auc_score

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache"
OUT = ROOT / "results" / "agod" / "llm_infer_align"
DOCS = ROOT / "docs" / "biz"

# ---------------------------------------------------------------------------
# Package you call
# ---------------------------------------------------------------------------
from agod.online_rfperm import (  # noqa: E402
    fit_online_probe,
    hop_fires,
    po_risk0_rows,
    probe_err,
    quantify_last_two,
    run_rfperm_stream,
    shift_ratio,
    error_floor,
)
from agod.po_iptw import po_iptw_weights  # noqa: E402
from agod.po_refit import Stream  # noqa: E402


# =============================================================================
# PO 核心逻辑（你要贴的那段）
# =============================================================================
#
# 1) 制度探针：上一批 fit μ0（浅层 RF）
# 2) 连续 OOS：e_now / e_prev ≥ γ  → fire
# 3) fire 时：对当前批算 po_risk0_i，审 Top-k；T=1 权重 w=√po_risk0
# 4) quiet 时：w=1（退火）
#
# 伪代码::
#
#   μ0 = fit_online_probe(X_prev, y_prev, task="acc")
#   e_now  = probe_err(μ0, X_cur, y_cur, task="acc")
#   fire   = hop_fires(e_now, e_prev, gate=1.25)
#   po     = po_risk0_rows(μ0, X_cur, y_cur, task="acc")
#   if fire:
#       w = po_iptw_weights(po, mode="sqrt")   # √PO，mean≈1
#       audit_idx = argsort(-po)[:k]
#   else:
#       w = ones_like(po)
#


def po_gate_step(
    X_prev: np.ndarray,
    y_prev: np.ndarray,
    X_cur: np.ndarray,
    y_cur: np.ndarray,
    *,
    e_prev: float | None,
    gate: float = 1.25,
    task: str = "acc",
    topk: int = 10,
    seed: int = 0,
) -> dict:
    """单步：探针 → 是否 fire → po_risk0 → √PO 权重 / Top-k 审计索引。"""
    probe = fit_online_probe(X_prev, y_prev, seed=seed, task=task)
    e_now = probe_err(probe, X_cur, y_cur, task=task)
    e_fl = error_floor(task, len(y_cur))
    fired = hop_fires(e_now, e_prev, gate=gate, e_floor=e_fl)
    ratio = (
        1.0
        if e_prev is None
        else shift_ratio(e_now, e_prev, e_floor=e_fl)
    )
    po = po_risk0_rows(probe, X_cur, y_cur, task=task)
    if fired:
        w = po_iptw_weights(po, mode="sqrt")
    else:
        w = np.ones_like(po, dtype=float)
    order = np.argsort(-po)
    audit_idx = order[: min(int(topk), len(order))].tolist()
    return {
        "fired": bool(fired),
        "ratio": float(ratio),
        "e_now": float(e_now),
        "e_prev": None if e_prev is None else float(e_prev),
        "po_risk0": np.asarray(po, dtype=float),
        "w": np.asarray(w, dtype=float),
        "audit_topk_local": audit_idx,
        "mean_po": float(np.mean(po)),
        "mean_w": float(np.mean(w)),
        "precision_at_k": float(
            np.mean(np.asarray(y_cur).ravel()[audit_idx].astype(int))
        )
        if len(audit_idx)
        else float("nan"),
    }


# =============================================================================
# 数据 IO
# =============================================================================


def load_jsonl(path: Path, n: int | None = None) -> list[dict]:
    rows: list[dict] = []
    with path.open() as f:
        for line in f:
            rows.append(json.loads(line))
            if n is not None and len(rows) >= n:
                break
    return rows


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


def ensure_halu(n: int = 3000) -> Path:
    path = CACHE / "halueval_qa_3000.jsonl"
    if path.exists():
        return path
    CACHE.mkdir(parents=True, exist_ok=True)
    from datasets import load_dataset

    ds = load_dataset("pminervini/HaluEval", "qa_samples", split="data", streaming=True)
    with path.open("w") as f:
        for i, row in enumerate(ds):
            f.write(
                json.dumps(
                    {
                        "knowledge": row.get("knowledge", ""),
                        "question": row.get("question", ""),
                        "answer": row.get("answer", ""),
                        "hallucination": row.get("hallucination", ""),
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )
            if i + 1 >= n:
                break
    return path


def hash_texts(texts: list[str], n_features: int = 192) -> np.ndarray:
    vec = HashingVectorizer(
        n_features=n_features, alternate_sign=False, norm="l2", ngram_range=(1, 2)
    )
    return vec.transform(texts).toarray().astype(float)


def pack_batches(X, y, n_per: int):
    n_use = (len(y) // n_per) * n_per
    X, y = X[:n_use], y[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)
    return X, y, batch


def auroc(y, s) -> float:
    y = np.asarray(y).astype(int)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, s))


def precision_at_k(y, s, k=10) -> float:
    y, s = np.asarray(y).astype(int), np.asarray(s, dtype=float)
    k = min(int(k), len(y))
    if k <= 0:
        return float("nan")
    return float(np.mean(y[np.argsort(-s)[:k]]))


# =============================================================================
# 推理 prototype：HaluEval 幻觉制度
# =============================================================================


def build_infer_stream(rows: list[dict], *, n_per=100, cut_batch=4, seed=0):
    """前半忠实、后半幻觉暴涨 → OnlineRFPerm 应 fire。"""
    rng = np.random.default_rng(seed)
    texts, y, rag = [], [], []
    for i, r in enumerate(rows):
        know, q, ans = r.get("knowledge", ""), r.get("question", ""), r.get("answer", "")
        hall = str(r.get("hallucination", "")).strip().lower() in ("yes", "true", "1")
        texts.append(f"{q}\n{ans}")
        y.append(1 if hall else 0)
        kset = set(re.findall(r"[a-z0-9]+", know.lower()))
        aset = set(re.findall(r"[a-z0-9]+", ans.lower()))
        rag.append(len(kset & aset) / max(len(aset), 1))
    X = hash_texts(texts)
    y = np.asarray(y, dtype=int)
    rag = np.asarray(rag, dtype=float)
    X, y, batch = pack_batches(X, y, n_per)
    rag = rag[: len(y)]
    # force a regime hop after cut_batch for demo clarity
    n_b = int(batch.max()) + 1
    for b in range(n_b):
        m = batch == b
        idx = np.flatnonzero(m)
        if b < cut_batch:
            # mostly faithful
            flip = idx[y[idx] == 1]
            drop = flip[: max(0, len(flip) - max(2, len(idx) // 20))]
            y[drop] = 0
        else:
            # mostly hallucinated
            keep = idx[y[idx] == 0]
            y[keep[: int(0.85 * len(keep))]] = 1
            # scramble rag a bit post-cut
            rag[idx] = np.clip(rag[idx] * rng.uniform(0.2, 0.7, size=idx.size), 0, 1)
    return {"X": X, "y": y, "batch": batch, "rag": rag, "cut_batch": cut_batch}


def run_infer_proto(*, gate=1.25, n_rows=1200) -> dict:
    path = ensure_halu()
    rows = load_jsonl(path, n=n_rows)
    stream = build_infer_stream(rows)
    X, y, batch = stream["X"], stream["y"], stream["batch"]
    rag = stream["rag"]

    # --- call package stream ---
    packed = run_rfperm_stream(
        Stream(X=X, y=y, batch=batch), gate=gate, learner="rf", seed=0, detail=True
    )
    hops = packed.get("history") or packed.get("hops") or []
    # find first fire at/after cut
    cut = stream["cut_batch"]
    fire_hop = next((h for h in hops if h.get("fired") and h.get("t", 0) >= cut - 1), None)
    if fire_hop is None:
        fire_hop = next((h for h in hops if h.get("fired")), hops[-1] if hops else {})

    # --- manual one-step PO gate (same math, explicit) ---
    t = int(fire_hop.get("t", cut))
    prev, cur = batch == (t - 1), batch == t
    step = po_gate_step(
        X[prev], y[prev], X[cur], y[cur], e_prev=None, gate=gate, topk=10, seed=t
    )
    # re-run with real e_prev from stream if available
    if t >= 2:
        pprev = batch == (t - 2)
        if pprev.any() and prev.any():
            probe0 = fit_online_probe(X[pprev], y[pprev], seed=t - 1, task="acc")
            e_prev = probe_err(probe0, X[prev], y[prev], task="acc")
            step = po_gate_step(
                X[prev],
                y[prev],
                X[cur],
                y[cur],
                e_prev=e_prev,
                gate=gate,
                topk=10,
                seed=t,
            )

    # ranking on T=1 of fire hop
    ranking = {}
    if fire_hop.get("po_risk0") is not None:
        idx = np.asarray(fire_hop["index"], dtype=int)
        treated = np.asarray(fire_hop["treated"], dtype=int).astype(bool)
        po = np.asarray(fire_hop["po_risk0"], dtype=float)
        t1 = idx[treated] if treated.any() else idx
        po_t1 = po[treated] if treated.any() else po
        y_t1 = y[t1]
        ranking = {
            "n_t1": int(len(y_t1)),
            "halluc_rate": float(np.mean(y_t1)),
            "auroc_po_risk0": auroc(y_t1, po_t1),
            "precision_at_10": precision_at_k(y_t1, po_t1, 10),
            "mean_rag_hit_top10": float(
                np.mean(rag[t1[np.argsort(-po_t1)[:10]]])
            )
            if len(t1)
            else float("nan"),
        }

    action = (
        "retrieval_refresh"
        if (ranking.get("mean_rag_hit_top10") or 1) < 0.35
        else "model_rollback_or_audit_topk"
    )

    return {
        "surface": "inference_hallucination",
        "package": "agod.online_rfperm / agod.po_iptw",
        "data": str(path.relative_to(ROOT)),
        "n": int(len(y)),
        "gate": gate,
        "stream_fires": int(sum(1 for h in hops if h.get("fired"))),
        "hop_at_cut": {
            "t": fire_hop.get("t"),
            "fired": fire_hop.get("fired"),
            "ratio": fire_hop.get("ratio"),
            "ranking": ranking,
        },
        "po_gate_step": {
            "fired": step["fired"],
            "ratio": step["ratio"],
            "mean_po": step["mean_po"],
            "mean_w": step["mean_w"],
            "precision_at_k": step["precision_at_k"],
            "audit_topk_local": step["audit_topk_local"],
        },
        "recommended_action": action,
        "how_to_use_in_serving": {
            "on_fire": "升配检索 / 人工审计 Top-k(po_risk0) / 可选 √PO 重加权下一跳训练",
            "on_quiet": "标准路径，w=1",
        },
    }


# =============================================================================
# 对齐 prototype：HH-RLHF RFPerm-as-Judge
# =============================================================================


def assistant_reply(dialog: str) -> str:
    parts = re.split(r"\n\nAssistant:", dialog)
    return parts[-1].strip() if len(parts) >= 2 else dialog.strip()


def build_align_stream(rows: list[dict], *, n_per=80, cut_batch=4, seed=0):
    """chosen=好 / rejected=坏；cut 后偏好映射翻转 → Judge 应 fire。"""
    rng = np.random.default_rng(seed)
    texts, y = [], []
    for r in rows:
        for key, lab in (("chosen", 1), ("rejected", 0)):
            texts.append(assistant_reply(r[key]))
            y.append(lab)
    X = hash_texts(texts)
    y = np.asarray(y, dtype=int)
    X, y, batch = pack_batches(X, y, n_per)
    # after cut: flip labels (preference map hopped)
    for b in range(int(batch.max()) + 1):
        if b >= cut_batch:
            m = batch == b
            y[m] = 1 - y[m]
            # add noise
            flip = np.flatnonzero(m)
            y[flip[rng.random(len(flip)) < 0.1]] ^= 1
    return {"X": X, "y": y, "batch": batch, "cut_batch": cut_batch}


def run_align_proto(*, gate=1.25, n_pairs=800) -> dict:
    path = ensure_hh()
    rows = load_jsonl(path, n=n_pairs)
    stream = build_align_stream(rows)
    X, y, batch = stream["X"], stream["y"], stream["batch"]

    packed = run_rfperm_stream(
        Stream(X=X, y=y, batch=batch), gate=gate, learner="rf", seed=1, detail=True
    )
    hops = packed.get("history") or packed.get("hops") or []
    cut = stream["cut_batch"]
    fire_hop = next((h for h in hops if h.get("fired") and h.get("t", 0) >= cut - 1), None)
    if fire_hop is None:
        fire_hop = next((h for h in hops if h.get("fired")), hops[-1] if hops else {})

    t = int(fire_hop.get("t", cut))
    prev, cur = batch == (t - 1), batch == t
    e_prev = None
    if t >= 2 and (batch == (t - 2)).any() and prev.any():
        p0 = fit_online_probe(X[batch == (t - 2)], y[batch == (t - 2)], seed=t, task="acc")
        e_prev = probe_err(p0, X[prev], y[prev], task="acc")
    step = po_gate_step(
        X[prev], y[prev], X[cur], y[cur], e_prev=e_prev, gate=gate, topk=10, seed=t
    )

    ranking = {}
    if fire_hop.get("po_risk0") is not None:
        idx = np.asarray(fire_hop["index"], dtype=int)
        treated = np.asarray(fire_hop["treated"], dtype=int).astype(bool)
        po = np.asarray(fire_hop["po_risk0"], dtype=float)
        t1 = idx[treated] if treated.any() else idx
        po_t1 = po[treated] if treated.any() else po
        ranking = {
            "n_t1": int(len(t1)),
            "auroc_po_risk0": auroc(y[t1], po_t1),
            "precision_at_10": precision_at_k(y[t1], po_t1, 10),
        }

    # merge gate decision
    merge = "REJECT_or_LIMIT" if fire_hop.get("fired") else "ACCEPT"
    return {
        "surface": "alignment_rfperm_as_judge",
        "package": "agod.online_rfperm / agod.po_iptw",
        "data": str(path.relative_to(ROOT)),
        "n": int(len(y)),
        "gate": gate,
        "stream_fires": int(sum(1 for h in hops if h.get("fired"))),
        "hop_at_cut": {
            "t": fire_hop.get("t"),
            "fired": fire_hop.get("fired"),
            "ratio": fire_hop.get("ratio"),
            "ranking": ranking,
        },
        "po_gate_step": {
            "fired": step["fired"],
            "ratio": step["ratio"],
            "mean_po": step["mean_po"],
            "mean_w": step["mean_w"],
            "precision_at_k": step["precision_at_k"],
            "audit_topk_local": step["audit_topk_local"],
        },
        "merge_gate": merge,
        "how_to_use_in_dpo_rlhf": {
            "on_fire": "拒合并或限量合并；po_risk0 Top-k 人工复核偏好对；可选 √PO 重加权",
            "on_quiet": "放行合并，w=1",
        },
    }


# =============================================================================
# Report
# =============================================================================


def render_md(infer: dict | None, align: dict | None) -> str:
    def blk(title: str, p: dict) -> str:
        h = p.get("hop_at_cut") or {}
        r = h.get("ranking") or {}
        s = p.get("po_gate_step") or {}
        return f"""### {title}

| 项 | 值 |
|----|-----|
| 数据 | `{p.get('data')}` |
| package | `{p.get('package')}` |
| n | {p.get('n')} |
| stream fires | {p.get('stream_fires')} |
| hop fired / ratio | `{h.get('fired')}` / `{h.get('ratio')}` |
| po_risk0 P@10 / AUROC | `{r.get('precision_at_10')}` / `{r.get('auroc_po_risk0')}` |
| 单步 gate fired / mean_w | `{s.get('fired')}` / `{s.get('mean_w')}` |
| 动作 / 门禁 | `{p.get('recommended_action') or p.get('merge_gate')}` |
"""

    parts = [
        "# LLM 推理 / 对齐 Prototype（call ``agod``）",
        "",
        "## 你 call 哪个包",
        "",
        "```python",
        "from agod.online_rfperm import (",
        "    fit_online_probe, hop_fires, po_risk0_rows,",
        "    probe_err, run_rfperm_stream,",
        ")",
        "from agod.po_iptw import po_iptw_weights",
        "from agod.po_refit import Stream",
        "```",
        "",
        "## 数据在哪",
        "",
        "| 用途 | 路径 | 源 |",
        "|------|------|----|",
        "| 对齐 / Judge | `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` | Anthropic/hh-rlhf |",
        "| 推理 / 幻觉 | `data/hf_cache/halueval_qa_3000.jsonl` | pminervini/HaluEval |",
        "",
        "## PO 核心逻辑（贴这段）",
        "",
        "```python",
        "probe = fit_online_probe(X_prev, y_prev, task='acc')",
        "e_now = probe_err(probe, X_cur, y_cur, task='acc')",
        "fire  = hop_fires(e_now, e_prev, gate=1.25)",
        "po    = po_risk0_rows(probe, X_cur, y_cur, task='acc')",
        "w     = po_iptw_weights(po, mode='sqrt') if fire else np.ones_like(po)",
        "audit = np.argsort(-po)[:10]   # Top-k 人工 / 升配",
        "```",
        "",
        "## 本跑结果",
        "",
    ]
    if infer:
        parts.append(blk("推理（幻觉制度）", infer))
        parts.append(
            f"- serving：`{(infer.get('how_to_use_in_serving') or {})}`\n"
        )
    if align:
        parts.append(blk("对齐（RFPerm-as-Judge）", align))
        parts.append(
            f"- DPO/RLHF：`{(align.get('how_to_use_in_dpo_rlhf') or {})}`\n"
        )
    parts += [
        "## 怎么跑",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/agod/llm_infer_align_prototype.py",
        "# → results/agod/llm_infer_align/",
        "# → docs/biz/LLM_INFER_ALIGN_PROTOTYPE.md",
        "```",
        "",
        "## 还有别的吗？",
        "",
        "| 还要什么 | 包 / 文件 |",
        "|----------|-----------|",
        "| 流式一整段 hop 历史 | `run_rfperm_stream(Stream(...))` |",
        "| 门控 √PO 权重 | `po_iptw_weights(..., mode='sqrt')` |",
        "| 更重的 refit 流 | `agod.po_refit.run_refit_stream` |",
        "| 接到客服¥账 | `scripts/agod/run_biz_value_sql_demo.py` |",
        "| 方法×情景异同 | `docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md` |",
        "",
        "边界：不是事实核查器；不是因果归因；fire 只打开对照窗 / 门禁。",
        "",
    ]
    return "\n".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["infer", "align", "both"], default="both")
    ap.add_argument("--gate", type=float, default=1.25)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    infer = align = None
    if args.only in ("infer", "both"):
        infer = run_infer_proto(gate=args.gate)
        (OUT / "infer.json").write_text(
            json.dumps(infer, ensure_ascii=False, indent=2, default=str)
        )
    if args.only in ("align", "both"):
        align = run_align_proto(gate=args.gate)
        (OUT / "align.json").write_text(
            json.dumps(align, ensure_ascii=False, indent=2, default=str)
        )

    summary = {"infer": infer, "align": align, "package": "agod", "gate": args.gate}
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str)
    )
    md = render_md(infer, align)
    (OUT / "LLM_INFER_ALIGN_PROTOTYPE.md").write_text(md)
    (DOCS / "LLM_INFER_ALIGN_PROTOTYPE.md").write_text(md)
    print(md)
    print(f"wrote {DOCS / 'LLM_INFER_ALIGN_PROTOTYPE.md'}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
