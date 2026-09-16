#!/usr/bin/env python3
"""HF landing prototypes — RFPerm-as-Judge, style hop, hallucination regime.

Cached subsets under ``data/hf_cache/``:

1. ``Anthropic/hh-rlhf`` helpful-base — preference / DPO-style continuous
   change / 文风(style) hop / RFPerm-as-Judge.
2. ``pminervini/HaluEval`` qa_samples — hallucination regime fire +
   ``po_risk0`` ranking + RAG-hit proxy + multi-task router signal.

Stance: regime detection, ranking, gated PO weights. Not fact-checking,
not unique causal attribution.

Usage::

    PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py
    PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py --only judge
    PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py --only halluc
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score

# Prefer site-packages ``datasets`` over workspace/datasets namespace shadow.
sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache"
OUT = ROOT / "results" / "agod" / "hf_landing"

from agod.online_rfperm import (  # noqa: E402
    fit_online_probe,
    po_risk0_rows,
    probe_err,
    run_rfperm_stream,
)
from agod.po_refit import Stream  # noqa: E402


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------

_HEDGE = ("maybe", "perhaps", "i think", "not sure", "probably", "might")
_FORMAL = ("therefore", "however", "furthermore", "regarding", "consequently")


def assistant_reply(dialog: str) -> str:
    parts = re.split(r"\n\nAssistant:", dialog)
    return parts[-1].strip() if len(parts) >= 2 else dialog.strip()


def human_prompt(dialog: str) -> str:
    m = re.search(r"Human:\s*(.*?)(?:\n\nAssistant:|$)", dialog, flags=re.S)
    return (m.group(1).strip() if m else dialog.strip())[:800]


def style_vector(text: str) -> np.ndarray:
    """Lightweight 文风 / register features (no LLM call)."""
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


def hash_text_matrix(
    texts: list[str], *, n_features: int = 256, seed: int = 0
) -> np.ndarray:
    _ = seed
    vec = HashingVectorizer(
        n_features=n_features,
        alternate_sign=False,
        norm="l2",
        ngram_range=(1, 2),
    )
    return vec.transform(texts).toarray().astype(float)


def pack_batches(
    X: np.ndarray, y: np.ndarray, n_per: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_use = (len(y) // n_per) * n_per
    X = X[:n_use]
    y = y[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)
    return X, y, batch


def auroc_safe(y: np.ndarray, s: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, s))


def precision_at_k(y: np.ndarray, s: np.ndarray, k: int = 10) -> float:
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    k = min(int(k), len(y))
    if k <= 0:
        return float("nan")
    idx = np.argsort(-s)[:k]
    return float(np.mean(y[idx]))


def domain_auc(X_a: np.ndarray, X_b: np.ndarray, seed: int = 0) -> float:
    X = np.vstack([X_a, X_b])
    y = np.concatenate([np.zeros(len(X_a)), np.ones(len(X_b))]).astype(int)
    if len(np.unique(y)) < 2:
        return float("nan")
    clf = RandomForestClassifier(
        n_estimators=40,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )
    scores = cross_val_score(clf, X, y, cv=3, scoring="roc_auc")
    return float(np.mean(scores))


def token_set(s: str) -> set[str]:
    return set(re.findall(r"[a-z0-9]+", (s or "").lower()))


def rag_hit_proxy(knowledge: str, answer: str) -> float:
    k, a = token_set(knowledge), token_set(answer)
    if not a:
        return 0.0
    return float(len(k & a) / max(len(a), 1))


def rank_hop(hop: dict | None, y_global: np.ndarray) -> dict | None:
    if not hop or "po_risk0" not in hop:
        return None
    idx = np.asarray(hop["index"], dtype=int)
    treated = np.asarray(hop["treated"], dtype=int).astype(bool)
    po = np.asarray(hop["po_risk0"], dtype=float)
    t1_idx = idx[treated] if treated.any() else idx
    po_t1 = po[treated] if treated.any() else po
    y_t1 = y_global[t1_idx]
    return {
        "n_t1": int(len(y_t1)),
        "pos_rate": float(np.mean(y_t1)),
        "auroc_po_risk0": auroc_safe(y_t1, po_t1),
        "precision_at_10": precision_at_k(y_t1, po_t1, 10),
        "mean_po_t1": float(np.mean(po_t1)),
    }


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if obj is None:
        return None
    return obj


# ---------------------------------------------------------------------------
# Demo 1 — HH-RLHF: RFPerm-as-Judge + preference hop + style hop
# ---------------------------------------------------------------------------


def build_hh_stream(
    rows: list[dict],
    *,
    n_per: int = 80,
    cut_batch: int = 4,
    n_features: int = 192,
    seed: int = 0,
) -> dict:
    """Expand chosen/rejected; casual→formal order; hard preference flip after cut."""
    rng = np.random.default_rng(seed)
    prompts, replies, labels, styles = [], [], [], []
    for r in rows:
        p = human_prompt(r["chosen"])
        for text, lab in (
            (assistant_reply(r["chosen"]), 1),
            (assistant_reply(r["rejected"]), 0),
        ):
            prompts.append(p)
            replies.append(text)
            labels.append(lab)
            styles.append(style_vector(text))

    styles = np.vstack(styles)
    labels = np.asarray(labels, dtype=int)
    # Shuffle so the preference hop is not confounded with gradual 文风 order.
    # Style portrait shift is measured separately via domain AUC on formality tails.
    perm = rng.permutation(len(labels))
    prompts = [prompts[i] for i in perm]
    replies = [replies[i] for i in perm]
    labels = labels[perm]
    styles = styles[perm]

    texts = [f"{p}\n\n{a}" for p, a in zip(prompts, replies)]
    X_txt = hash_text_matrix(texts, n_features=n_features, seed=seed)
    X = np.hstack([X_txt, styles])
    y = labels.copy()
    X, y, batch = pack_batches(X, y, n_per)
    styles = styles[: len(y)]

    # Hard preference-map hop (DPO/RLHF continuous-update failure mode).
    after = batch >= cut_batch
    flip = after & (rng.random(len(y)) < 0.92)
    y = y.copy()
    y[flip] = 1 - y[flip]

    return {
        "X": X,
        "y": y,
        "batch": batch,
        "styles": styles,
        "n_text": n_features,
        "cut_batch": cut_batch,
        "n_per": n_per,
        "flip_rate_after": float(flip.mean()) if after.any() else 0.0,
        "formality_early": float(np.mean(styles[batch < cut_batch][:, 6]))
        if np.any(batch < cut_batch)
        else float("nan"),
        "formality_late": float(np.mean(styles[batch >= cut_batch][:, 6]))
        if np.any(batch >= cut_batch)
        else float("nan"),
    }


def run_judge_demo(pack: dict, *, gate: float = 1.35, seed: int = 0) -> dict:
    X, y, batch = pack["X"], pack["y"], pack["batch"]
    cut = int(pack["cut_batch"])
    stream = Stream(X=X, y=y, batch=batch, name="hh_judge", task="acc")
    rec = run_rfperm_stream(stream, gate=gate, learner="rf", seed=seed, detail=True)

    pre, post = batch < cut, batch >= cut
    judge = fit_online_probe(X[pre], y[pre], seed=seed, task="acc")
    e_pre = float(probe_err(judge, X[pre], y[pre], task="acc"))
    e_post = float(probe_err(judge, X[post], y[post], task="acc"))
    po_post = po_risk0_rows(
        judge, X[post], y[post], batch_po=max(e_post - e_pre, 0.0), task="acc"
    )

    n_text = int(pack["n_text"])
    # Stream pre/post may be IID in X after shuffle — use formality tails as
    # the explicit 文风 portrait contrast (low vs high register terciles).
    formality = pack["styles"][:, 6] + 0.5 * pack["styles"][:, 1]
    q_lo, q_hi = np.quantile(formality, [0.33, 0.67])
    lo = formality <= q_lo
    hi = formality >= q_hi
    style_auc = domain_auc(pack["styles"][lo], pack["styles"][hi], seed=seed)
    text_auc = domain_auc(X[pre][:, :n_text], X[post][:, :n_text], seed=seed + 1)

    hist = rec["history"]
    fires = [h for h in hist if h.get("fired")]
    hop_cut = next((h for h in hist if h["t"] == cut), None)
    hop_quiet = next((h for h in hist if h["t"] == max(cut - 2, 1)), None)

    return {
        "scenario": "hh_rlhf_rfperm_as_judge_style",
        "dataset": "Anthropic/hh-rlhf@helpful-base",
        "n": int(len(y)),
        "n_batches": int(batch.max() + 1),
        "n_per": int(pack["n_per"]),
        "cut_batch": cut,
        "gate": gate,
        "flip_rate_after": pack["flip_rate_after"],
        "formality_early": pack["formality_early"],
        "formality_late": pack["formality_late"],
        "style_domain_auc": style_auc,
        "text_domain_auc": text_auc,
        "judge_err_pre": e_pre,
        "judge_err_post": e_post,
        "judge_err_ratio": float(e_post / max(e_pre, 1e-6)),
        "judge_post_po_auroc": auroc_safe(y[post], po_post),
        "fire_rate": float(rec.get("fire_rate", 0.0)),
        "first_fire_t": int(fires[0]["t"]) if fires else None,
        "fires": [
            {
                "t": h["t"],
                "ratio": h["ratio"],
                "mean_r0": h.get("mean_r0"),
                "mean_r1": h.get("mean_r1"),
            }
            for h in fires
        ],
        "hop_at_cut": {
            "t": hop_cut["t"],
            "fired": hop_cut["fired"],
            "ratio": hop_cut["ratio"],
            "ranking": rank_hop(hop_cut, y),
        }
        if hop_cut
        else None,
        "hop_quiet": {
            "t": hop_quiet["t"],
            "fired": hop_quiet["fired"],
            "ratio": hop_quiet["ratio"],
            "ranking": rank_hop(hop_quiet, y),
        }
        if hop_quiet
        else None,
        "reading": {
            "style_axis": "style_domain_auc high => register/文风 portrait moved (P(X))",
            "preference_axis": "judge_err_ratio / OnlineRFPerm fire => preference map P(Y|X) moved",
            "action": "quiet->uniform; fire->audit Top-k by po_risk0; optional sqrt(PO) on T=1",
            "not_a_claim": "not a moral judge; not causal feature attribution",
        },
    }


# ---------------------------------------------------------------------------
# Demo 2 — HaluEval: hallucination regime + RAG + router
# ---------------------------------------------------------------------------


def build_halu_stream(
    rows: list[dict],
    *,
    n_per: int = 100,
    cut_batch: int = 4,
    n_features: int = 192,
    seed: int = 0,
) -> dict:
    """Clean batches first, hallucinated mass after cut — sharp regime hop."""
    rng = np.random.default_rng(seed)
    y_raw = np.asarray(
        [
            1 if str(r.get("hallucination", "")).lower().startswith("y") else 0
            for r in rows
        ],
        dtype=int,
    )
    # Put faithful answers first, hallucinations last (stable within class).
    order = np.argsort(y_raw, kind="mergesort")
    rows = [rows[i] for i in order]
    y_raw = y_raw[order]

    texts = [
        f"Q: {r.get('question', '')}\nK: {r.get('knowledge', '')}\nA: {r.get('answer', '')}"
        for r in rows
    ]
    X_txt = hash_text_matrix(texts, n_features=n_features, seed=seed)
    rag = np.asarray(
        [rag_hit_proxy(r.get("knowledge", ""), r.get("answer", "")) for r in rows],
        dtype=float,
    ).reshape(-1, 1)
    styles = np.vstack([style_vector(r.get("answer", "")) for r in rows])
    X = np.hstack([X_txt, rag, styles])

    X, y_raw, batch = pack_batches(X, y_raw, n_per)
    rag = rag[: len(y_raw)].ravel()
    rows = rows[: len(y_raw)]

    # Hard regime with non-vacuous pre-cut denominator:
    # sprinkle a few positives before cut so e_prev is not ~0 (OnlineRFPerm
    # refuses vacuous ratios). Then jump to mostly-halluc after cut.
    y = y_raw.copy()
    pre = batch < cut_batch
    post = batch >= cut_batch
    y[pre] = 0
    sprinkle = pre & (rng.random(len(y)) < 0.08)
    y[sprinkle] = 1
    convert = post & (y_raw == 0) & (rng.random(len(y)) < 0.85)
    y[post] = np.where(y_raw[post] == 1, 1, 0)
    y[convert] = 1

    qlen = np.asarray([len(token_set(r.get("question", ""))) for r in rows], dtype=float)
    task_id = (qlen > np.median(qlen)).astype(int)

    return {
        "X": X,
        "y": y,
        "batch": batch,
        "rag": rag,
        "task_id": task_id,
        "cut_batch": cut_batch,
        "n_per": n_per,
        "rate_before": float(y[pre].mean()) if pre.any() else float("nan"),
        "rate_after": float(y[post].mean()) if post.any() else float("nan"),
    }


def run_halluc_demo(pack: dict, *, gate: float = 1.35, seed: int = 0) -> dict:
    X, y, batch = pack["X"], pack["y"], pack["batch"]
    cut = int(pack["cut_batch"])
    stream = Stream(X=X, y=y, batch=batch, name="halu_qa", task="acc")
    rec = run_rfperm_stream(stream, gate=gate, learner="rf", seed=seed, detail=True)

    hist = rec["history"]
    fires = [h for h in hist if h.get("fired")]
    hop_cut = next((h for h in hist if h["t"] == cut), None)
    hop_quiet = next((h for h in hist if h["t"] == max(cut - 2, 1)), None)

    def _rank(hop: dict | None) -> dict | None:
        base = rank_hop(hop, y)
        if not base or not hop:
            return base
        idx = np.asarray(hop["index"], dtype=int)
        treated = np.asarray(hop["treated"], dtype=int).astype(bool)
        po = np.asarray(hop["po_risk0"], dtype=float)
        t1_idx = idx[treated] if treated.any() else idx
        po_t1 = po[treated] if treated.any() else po
        rag_t1 = pack["rag"][t1_idx]
        top = np.argsort(-po_t1)[: min(10, len(po_t1))]
        out = dict(base)
        out["halluc_rate"] = out.pop("pos_rate")
        out["mean_rag_hit_top10"] = (
            float(np.mean(rag_t1[top])) if len(top) else float("nan")
        )
        return out

    pre, post = batch < cut, batch >= cut
    tid = pack["task_id"]
    task_shift = {
        "task0_rate_pre": float(np.mean(tid[pre] == 0)) if pre.any() else float("nan"),
        "task0_rate_post": float(np.mean(tid[post] == 0)) if post.any() else float("nan"),
        "halluc_by_task_post": {
            "task0": float(np.mean(y[post & (tid == 0)]))
            if np.any(post & (tid == 0))
            else float("nan"),
            "task1": float(np.mean(y[post & (tid == 1)]))
            if np.any(post & (tid == 1))
            else float("nan"),
        },
        "domain_auc_pre_post": domain_auc(X[pre], X[post], seed=seed),
    }

    probe = fit_online_probe(X[pre], y[pre], seed=seed, task="acc")
    po_post = po_risk0_rows(probe, X[post], y[post], task="acc")

    return {
        "scenario": "halueval_hallucination_regime_rag",
        "dataset": "pminervini/HaluEval@qa_samples",
        "n": int(len(y)),
        "n_batches": int(batch.max() + 1),
        "n_per": int(pack["n_per"]),
        "cut_batch": cut,
        "gate": gate,
        "rate_before": pack["rate_before"],
        "rate_after": pack["rate_after"],
        "fire_rate": float(rec.get("fire_rate", 0.0)),
        "first_fire_t": int(fires[0]["t"]) if fires else None,
        "fires": [
            {
                "t": h["t"],
                "ratio": h["ratio"],
                "mean_r0": h.get("mean_r0"),
                "mean_r1": h.get("mean_r1"),
            }
            for h in fires
        ],
        "hop_at_cut": {
            "t": hop_cut["t"],
            "fired": hop_cut["fired"],
            "ratio": hop_cut["ratio"],
            "ranking": _rank(hop_cut),
        }
        if hop_cut
        else None,
        "hop_quiet": {
            "t": hop_quiet["t"],
            "fired": hop_quiet["fired"],
            "ratio": hop_quiet["ratio"],
            "ranking": _rank(hop_quiet),
        }
        if hop_quiet
        else None,
        "instance_probe_auroc_post": auroc_safe(y[post], po_post),
        "router_task_shift": task_shift,
        "reading": {
            "regime": "OnlineRFPerm fire at cut => hallucination label law hopped",
            "ranking": "po_risk0 sorts audit candidates under shift (not a fact checker)",
            "rag": "low rag_hit among Top-k => retrieval-support gap, separate from generation hop",
            "router": "task_shift guides which expert/route to refresh next",
        },
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def write_markdown(summary: dict, path: Path) -> None:
    j = summary["judge"]
    h = summary["halluc"]
    hop_j = j.get("hop_at_cut") or {}
    hop_h = h.get("hop_at_cut") or {}
    lines = [
        "# HF landing prototypes — run report",
        "",
        "Packages: local `agod` OnlineRFPerm + PO-risk; sklearn HashingVectorizer RF.",
        "Datasets: `Anthropic/hh-rlhf` (helpful-base subset), `pminervini/HaluEval` (qa_samples).",
        "",
        "## 1) RFPerm-as-Judge + 文风 / preference hop (HH-RLHF)",
        "",
        f"- judge err pre→post: `{j['judge_err_pre']:.3f}` → `{j['judge_err_post']:.3f}` "
        f"(ratio `{j['judge_err_ratio']:.2f}`)",
        f"- style domain AUC (P(X) register): `{j['style_domain_auc']:.3f}`",
        f"- text domain AUC: `{j['text_domain_auc']:.3f}`",
        f"- fire_rate `{j['fire_rate']:.3f}`, first_fire_t `{j['first_fire_t']}`",
        f"- hop@cut fired=`{hop_j.get('fired')}` ratio=`{hop_j.get('ratio')}`",
        "",
        "Read: high style AUC ⇒ 文风画像变了；judge ratio / fire ⇒ 偏好映射变了（DPO/RLHF 连更）。",
        "动作：fire 后按 `po_risk0` Top-k 抽检偏好对；安静回均匀。",
        "",
        "## 2) Hallucination regime + RAG-hit + router (HaluEval)",
        "",
        f"- halluc rate before/after: `{h['rate_before']:.3f}` / `{h['rate_after']:.3f}`",
        f"- fire_rate `{h['fire_rate']:.3f}`, first_fire_t `{h['first_fire_t']}`",
        f"- instance probe AUROC(post): `{h['instance_probe_auroc_post']}`",
        f"- hop@cut fired=`{hop_h.get('fired')}` ranking=`{hop_h.get('ranking')}`",
        f"- router_task_shift: `{h['router_task_shift']}`",
        "",
        "Read: regime fire = 幻觉标签制度跳变；`po_risk0` 只排序；RAG-hit 低 ⇒ 检索支持缺口；",
        "task_shift ⇒ 下一步刷新哪个 router 专家。",
        "",
        "## Stance",
        "",
        "Detection / localization / ranking / gated reweight — not fact-checking, not unique causal attribution.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["all", "judge", "halluc"], default="all")
    ap.add_argument("--gate", type=float, default=1.35)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    summary: dict = {"gate": args.gate, "seed": args.seed}

    if args.only in ("all", "judge"):
        hh_rows = load_jsonl(ensure_hh(), n=2000)
        pack = build_hh_stream(hh_rows, n_per=80, cut_batch=4, seed=args.seed)
        judge = run_judge_demo(pack, gate=args.gate, seed=args.seed)
        summary["judge"] = judge
        (OUT / "hh_judge_style.json").write_text(
            json.dumps(jsonable(judge), indent=2), encoding="utf-8"
        )
        print(
            f"[judge] err {judge['judge_err_pre']:.3f}→{judge['judge_err_post']:.3f} "
            f"ratio={judge['judge_err_ratio']:.2f} style_auc={judge['style_domain_auc']:.3f} "
            f"fire_rate={judge['fire_rate']:.3f} first_fire={judge['first_fire_t']}"
        )

    if args.only in ("all", "halluc"):
        hu_rows = load_jsonl(ensure_halu(), n=2400)
        pack_h = build_halu_stream(hu_rows, n_per=100, cut_batch=4, seed=args.seed)
        halluc = run_halluc_demo(pack_h, gate=args.gate, seed=args.seed)
        summary["halluc"] = halluc
        (OUT / "halu_regime_rag.json").write_text(
            json.dumps(jsonable(halluc), indent=2), encoding="utf-8"
        )
        print(
            f"[halluc] rate {halluc['rate_before']:.3f}→{halluc['rate_after']:.3f} "
            f"fire_rate={halluc['fire_rate']:.3f} first_fire={halluc['first_fire_t']} "
            f"auroc={halluc['instance_probe_auroc_post']}"
        )

    if args.only == "all":
        (OUT / "summary.json").write_text(
            json.dumps(jsonable(summary), indent=2), encoding="utf-8"
        )
        write_markdown(summary, OUT / "REPORT.md")
        print(f"wrote {OUT / 'REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
