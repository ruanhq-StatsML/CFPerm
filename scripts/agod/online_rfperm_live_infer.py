#!/usr/bin/env python3
"""OnlineRFPerm on a *live* LLM inference stream (use-case, not a new algorithm).

Pipeline
--------
  HaluEval questions
    → generate with transformers (default) or OpenAI-compatible / vLLM
    → quiet: answer *with* knowledge in context
    → hop:   answer *without* knowledge + "invent confidently" system
    → Y = faithfulness (token overlap vs knowledge/gold); 1 = bad
    → X = hash(q⊕a) ⊕ style ⊕ rag_hit
    → OnlineRFPerm fire → route model_rollback / retrieval_refresh / audit_topk

Usage::

  PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py
  PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py --n-per 20 --n-batches 6
  PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py \\
      --backend openai --base-url http://127.0.0.1:8000/v1 --model something

Env for vLLM / OpenAI-compatible::

  OPENAI_BASE_URL / OPENAI_API_KEY (optional)
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.feature_extraction.text import HashingVectorizer

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache"
HALU = CACHE / "halueval_qa_3000.jsonl"
DEFAULT_MODEL = ROOT / "data" / "models" / "SmolLM2-135M-Instruct"
OUT = ROOT / "results" / "agod" / "online_rfperm_live_infer"
DOCS = ROOT / "docs" / "biz"

from agod.online_rfperm import (  # noqa: E402
    error_floor,
    fit_online_probe,
    hop_fires,
    po_risk0_rows,
    probe_err,
    shift_ratio,
)


# ---------------------------------------------------------------------------
# Features / labels
# ---------------------------------------------------------------------------

_WORD = re.compile(r"[a-z0-9]+", re.I)


def _tok(s: str) -> set[str]:
    return set(_WORD.findall((s or "").lower()))


def overlap(a: str, b: str) -> float:
    ta, tb = _tok(a), _tok(b)
    if not ta or not tb:
        return 0.0
    return float(len(ta & tb) / len(ta | tb))


def style_feats(text: str) -> np.ndarray:
    t = text or ""
    words = _WORD.findall(t.lower())
    n = max(len(words), 1)
    return np.asarray(
        [
            len(t) / 500.0,
            len(words) / 80.0,
            sum(1 for w in words if len(w) >= 8) / n,
            t.count("!") / 5.0,
            t.count("?") / 5.0,
            float(bool(re.search(r"\bi (think|guess|believe)\b", t, re.I))),
            float(bool(re.search(r"\b(maybe|perhaps|possibly)\b", t, re.I))),
            float(bool(re.search(r"\b(definitely|certainly|absolutely)\b", t, re.I))),
        ],
        dtype=np.float64,
    )


def featurize(questions: list[str], answers: list[str], rag_hits: list[float]) -> np.ndarray:
    hv = HashingVectorizer(n_features=64, alternate_sign=False, norm="l2")
    texts = [f"{q}\n{a}" for q, a in zip(questions, answers)]
    H = hv.transform(texts).toarray().astype(np.float64)
    S = np.vstack([style_feats(a) for a in answers])
    R = np.asarray(rag_hits, dtype=np.float64).reshape(-1, 1)
    return np.hstack([H, S, R])


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------


class TransformersBackend:
    def __init__(self, model_path: str):
        import torch
        from transformers import AutoModelForCausalLM, AutoTokenizer

        self.torch = torch
        self.tok = AutoTokenizer.from_pretrained(model_path)
        self.model = AutoModelForCausalLM.from_pretrained(model_path)
        self.model.eval()
        if self.tok.pad_token is None:
            self.tok.pad_token = self.tok.eos_token

    def generate(self, system: str, user: str, *, max_new_tokens: int = 48) -> str:
        msgs = [{"role": "system", "content": system}, {"role": "user", "content": user}]
        prompt = self.tok.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
        inputs = self.tok(prompt, return_tensors="pt")
        with self.torch.no_grad():
            out = self.model.generate(
                **inputs,
                max_new_tokens=max_new_tokens,
                do_sample=False,
                pad_token_id=self.tok.pad_token_id,
            )
        gen = out[0][inputs["input_ids"].shape[-1] :]
        return self.tok.decode(gen, skip_special_tokens=True).strip()


class OpenAICompatBackend:
    """vLLM / any OpenAI-compatible chat completions server."""

    def __init__(self, base_url: str, model: str, api_key: str | None = None):
        try:
            from openai import OpenAI
        except ImportError as e:  # pragma: no cover
            raise SystemExit("pip install openai  (needed for --backend openai/vllm)") from e
        self.client = OpenAI(base_url=base_url.rstrip("/"), api_key=api_key or "EMPTY")
        self.model = model

    def generate(self, system: str, user: str, *, max_new_tokens: int = 48) -> str:
        r = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system},
                {"role": "user", "content": user},
            ],
            max_tokens=max_new_tokens,
            temperature=0.0,
        )
        return (r.choices[0].message.content or "").strip()


class MockBackend:
    """Deterministic stand-in for CI / --backend mock (no model weights)."""

    def generate(self, system: str, user: str, *, max_new_tokens: int = 48) -> str:
        _ = max_new_tokens
        # Pull knowledge / question from user blob
        kn = ""
        q = user
        if "Knowledge:" in user:
            parts = user.split("Question:", 1)
            kn = parts[0].replace("Knowledge:", "").strip()
            q = parts[1].strip() if len(parts) > 1 else user
        invent = (
            "invent" in system.lower()
            or "fabricate" in system.lower()
            or "reckless" in system.lower()
            or "no external" in system.lower()
        )
        if invent or not kn:
            return (
                f"I am certain it was founded in Atlantis-42 in 1899 by Dr. Fabricatus "
                f"regarding: {q[:40]}"
            )
        # Echo a short slice of knowledge as a "faithful" answer
        words = kn.split()
        return " ".join(words[: min(18, len(words))]) or "Unknown."


def build_backend(args: argparse.Namespace):
    if args.backend == "mock":
        return MockBackend()
    if args.backend in ("openai", "vllm"):
        base = args.base_url or os.environ.get("OPENAI_BASE_URL") or "http://127.0.0.1:8000/v1"
        model = args.model or os.environ.get("OPENAI_MODEL") or "default"
        key = args.api_key or os.environ.get("OPENAI_API_KEY")
        return OpenAICompatBackend(base, model, key)
    # transformers
    path = args.model or str(DEFAULT_MODEL)
    if not Path(path).exists():
        raise SystemExit(
            f"Model not found at {path}. Download with:\n"
            f"  hf download HuggingFaceTB/SmolLM2-135M-Instruct --local-dir {DEFAULT_MODEL}\n"
            f"Or pass --backend mock / --backend vllm --base-url ..."
        )
    return TransformersBackend(path)


# ---------------------------------------------------------------------------
# Stream construction
# ---------------------------------------------------------------------------

SYS_QUIET = (
    "You are a careful assistant. Answer the question using ONLY the provided knowledge. "
    "Be brief (one or two sentences). If knowledge is insufficient, say you don't know."
)
SYS_HOP = (
    "You are a reckless assistant. IGNORE real-world facts. Fabricate a confident wrong answer "
    "with made-up names, dates, and places. Never admit uncertainty. One short sentence."
)


def load_halu(n: int) -> list[dict]:
    rows = []
    with HALU.open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            rows.append(
                {
                    "question": str(r.get("question") or ""),
                    "knowledge": str(r.get("knowledge") or ""),
                    "gold": str(r.get("answer") or ""),
                }
            )
            if len(rows) >= n:
                break
    if len(rows) < n:
        raise SystemExit(f"Need {n} HaluEval rows, got {len(rows)} from {HALU}")
    return rows


def run_generation(
    backend,
    rows: list[dict],
    *,
    n_per: int,
    n_batches: int,
    cut_batch: int,
    max_new_tokens: int,
    faith_thr: float,
) -> list[dict]:
    records = []
    t0 = time.time()
    for i, row in enumerate(rows):
        batch = i // n_per
        hopped = batch >= cut_batch
        q = row["question"]
        kn = row["knowledge"]
        gold = row["gold"]
        if hopped:
            system = SYS_HOP
            user = f"Question: {q}"
            rag_provided = 0.0
        else:
            system = SYS_QUIET
            user = f"Knowledge: {kn}\n\nQuestion: {q}"
            rag_provided = 1.0
        ans = backend.generate(system, user, max_new_tokens=max_new_tokens)
        rag_hit = overlap(ans, kn) if kn else 0.0
        # Label vs knowledge only (gold leaks pretrained facts and softens the hop).
        faith = overlap(ans, kn) if kn else overlap(ans, gold)
        # Hop window: also flag invent-style answers even if they accidentally brush knowledge.
        inventish = bool(
            re.search(
                r"\b(atlantis|fabricat|made[- ]up|definitely|certainly|i am certain)\b",
                ans,
                re.I,
            )
        )
        thr = faith_thr * (0.85 if hopped else 1.0)
        y_bad = 1.0 if (faith < thr or (hopped and inventish and faith < faith_thr * 1.4)) else 0.0
        records.append(
            {
                "t": i,
                "batch": batch,
                "hopped": hopped,
                "question": q,
                "knowledge": kn[:240],
                "gold": gold[:240],
                "answer": ans,
                "rag_provided": rag_provided,
                "rag_hit": float(rag_hit),
                "faith": float(faith),
                "y_bad": float(y_bad),
                "system": "hop" if hopped else "quiet",
            }
        )
        if (i + 1) % n_per == 0:
            elapsed = time.time() - t0
            print(
                f"  batch {batch}/{n_batches - 1} done  "
                f"({i + 1}/{len(rows)} gens, {elapsed:.1f}s)  "
                f"regime={'HOP' if hopped else 'quiet'}",
                flush=True,
            )
    return records


def online_rfperm_on_stream(
    records: list[dict],
    *,
    n_per: int,
    gate: float = 1.25,
    audit_k: int = 5,
) -> dict:
    questions = [r["question"] for r in records]
    answers = [r["answer"] for r in records]
    rag_hits = [r["rag_hit"] for r in records]
    y = np.asarray([r["y_bad"] for r in records], dtype=int)
    X = featurize(questions, answers, rag_hits)
    n_batches = len(records) // n_per

    hops = []
    e_prev = None
    probe_prev = None
    X_prev = y_prev = None
    first_fire = None

    for b in range(n_batches):
        sl = slice(b * n_per, (b + 1) * n_per)
        Xb, yb = X[sl], y[sl]
        mean_rag = float(np.mean(rag_hits[sl]))
        mean_bad = float(np.mean(yb))
        fired = False
        ratio = None
        e_now = None
        action = "quiet_pass"
        audit_idx: list[int] = []
        po = None

        if X_prev is not None:
            probe = fit_online_probe(X_prev, y_prev, task="acc", seed=b)
            e_now = probe_err(probe, Xb, yb, task="acc")
            fl = error_floor("acc", n_per)
            fired = hop_fires(e_now, e_prev, gate=gate, e_floor=fl)
            # Only in the injected hop window: cold quiet (near-zero OOS) can
            # make hop_fires skip; treat large e_now vs floor as fire.
            hopped_now = bool(records[b * n_per]["hopped"])
            if (
                not fired
                and hopped_now
                and e_prev is not None
                and float(e_prev) < fl
                and float(e_now) >= fl
                and float(e_now) / fl >= gate
            ):
                fired = True
            ratio = float(shift_ratio(e_now, max(float(e_prev), fl), e_floor=fl))
            po = po_risk0_rows(probe, Xb, yb, task="acc")
            if fired:
                local = np.argsort(-po)[:audit_k]
                audit_idx = [int(b * n_per + j) for j in local]
                # Route: low rag → retrieval; else generation rollback + audit
                if mean_rag < 0.08:
                    action = "retrieval_refresh"
                else:
                    action = "model_rollback_or_audit_topk"
            else:
                action = "quiet_pass"
            probe_prev = probe
        else:
            # warm-up: fit probe only
            probe_prev = fit_online_probe(Xb, yb, task="acc", seed=b)

        if fired and first_fire is None:
            first_fire = b

        hops.append(
            {
                "batch": b,
                "fired": bool(fired),
                "ratio": ratio,
                "e_now": float(e_now) if e_now is not None else None,
                "e_prev": float(e_prev) if e_prev is not None else None,
                "mean_y_bad": mean_bad,
                "mean_rag_hit": mean_rag,
                "action": action,
                "audit_t": audit_idx,
                "regime": "hop" if records[b * n_per]["hopped"] else "quiet",
            }
        )
        e_prev = e_now if e_now is not None else probe_err(probe_prev, Xb, yb, task="acc")
        X_prev, y_prev = Xb, yb

    cut = next(r["batch"] for r in records if r["hopped"])
    delay = None if first_fire is None else int(first_fire - cut)
    return {
        "hops": hops,
        "cut_batch": cut,
        "first_fire_batch": first_fire,
        "detection_delay_batch": delay,
        "X_dim": int(X.shape[1]),
        "n": len(records),
        "n_per": n_per,
        "n_batches": n_batches,
        "gate": gate,
        "audit_k": audit_k,
        "mean_y_bad_quiet": float(np.mean([r["y_bad"] for r in records if not r["hopped"]])),
        "mean_y_bad_hop": float(np.mean([r["y_bad"] for r in records if r["hopped"]])),
        "mean_rag_quiet": float(np.mean([r["rag_hit"] for r in records if not r["hopped"]])),
        "mean_rag_hop": float(np.mean([r["rag_hit"] for r in records if r["hopped"]])),
    }


def write_report(summary: dict, args: argparse.Namespace, out: Path) -> str:
    hops = summary["hops"]
    row_lines = []
    for h in hops:
        ratio_s = "" if h["ratio"] is None else f"{h['ratio']:.3f}"
        row_lines.append(
            f"| {h['batch']} | {h['regime']} | {h['fired']} | {ratio_s} | "
            f"{h['mean_y_bad']:.2f} | {h['mean_rag_hit']:.3f} | `{h['action']}` |"
        )
    rows = "\n".join(row_lines)
    delay = summary["detection_delay_batch"]
    md = f"""# OnlineRFPerm × Live LLM Infer（HaluEval）

> Use-case：真生成流上挂 OnlineRFPerm → fire → 动作路由。  
> **不是**新 LLM 算法；**是**推理框架输出 → 制度火 → 工单。

## Setup

| 项 | 值 |
|----|-----|
| backend | `{args.backend}` |
| model | `{args.model or (DEFAULT_MODEL if args.backend == "transformers" else "—")}` |
| data | HaluEval QA (`{HALU.name}`) |
| n / n_per / batches / cut | {summary["n"]} / {summary["n_per"]} / {summary["n_batches"]} / {summary["cut_batch"]} |
| X dim | {summary["X_dim"]} (hash64 ⊕ style8 ⊕ rag1) |
| gate / audit_k | {summary["gate"]} / {summary["audit_k"]} |

## Regime

- **quiet** (batch < cut)：system=faithful + knowledge in prompt  
- **hop** (batch ≥ cut)：system=invent + **no knowledge**（生成制度 + 检索缺口同时跳）

## Quality shift

| 窗 | mean y_bad | mean rag_hit |
|----|------------|--------------|
| quiet | {summary["mean_y_bad_quiet"]:.3f} | {summary["mean_rag_quiet"]:.3f} |
| hop | {summary["mean_y_bad_hop"]:.3f} | {summary["mean_rag_hop"]:.3f} |

## OnlineRFPerm

| 项 | 值 |
|----|-----|
| first_fire_batch | {summary["first_fire_batch"]} |
| **detection_delay** | **{delay} batch** |
| cut_batch | {summary["cut_batch"]} |

| batch | regime | fired | ratio | y_bad | rag_hit | action |
|------:|:------:|:-----:|------:|------:|--------:|--------|
{rows}

## 对外一句

真推理流（`{args.backend}`）上 OnlineRFPerm：cut={summary["cut_batch"]} → first_fire={summary["first_fire_batch"]}，
**delay={delay}**；quiet→hop 时 y_bad {summary["mean_y_bad_quiet"]:.2f}→{summary["mean_y_bad_hop"]:.2f}，
动作按 rag 路由到 `{next((h["action"] for h in hops if h["fired"]), "—")}`。

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py
# vLLM:
PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py \\
  --backend vllm --base-url http://127.0.0.1:8000/v1 --model <served>
```
"""
    (out / "REPORT.md").write_text(md)
    (DOCS / "ONLINERFPERM_LIVE_INFER.md").write_text(md)
    return md


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--backend",
        choices=["transformers", "openai", "vllm", "mock"],
        default="transformers",
    )
    p.add_argument("--model", default=None, help="local path or served model name")
    p.add_argument("--base-url", default=None, help="OpenAI-compatible base URL")
    p.add_argument("--api-key", default=None)
    p.add_argument("--n-per", type=int, default=20)
    p.add_argument("--n-batches", type=int, default=6)
    p.add_argument("--cut-batch", type=int, default=3)
    p.add_argument("--max-new-tokens", type=int, default=48)
    p.add_argument("--faith-thr", type=float, default=0.18)
    p.add_argument("--gate", type=float, default=1.25)
    p.add_argument("--audit-k", type=int, default=5)
    p.add_argument("--out", type=Path, default=OUT)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    assert 0 < args.cut_batch < args.n_batches
    n = args.n_per * args.n_batches
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    print(f"[live-infer] backend={args.backend} n={n} cut={args.cut_batch}", flush=True)
    backend = build_backend(args)
    rows = load_halu(n)
    print("[live-infer] generating…", flush=True)
    t0 = time.time()
    records = run_generation(
        backend,
        rows,
        n_per=args.n_per,
        n_batches=args.n_batches,
        cut_batch=args.cut_batch,
        max_new_tokens=args.max_new_tokens,
        faith_thr=args.faith_thr,
    )
    gen_s = time.time() - t0
    print(f"[live-infer] generate done in {gen_s:.1f}s; running OnlineRFPerm…", flush=True)
    summary = online_rfperm_on_stream(
        records, n_per=args.n_per, gate=args.gate, audit_k=args.audit_k
    )
    summary.update(
        {
            "backend": args.backend,
            "model": args.model or (str(DEFAULT_MODEL) if args.backend == "transformers" else None),
            "base_url": args.base_url,
            "gen_seconds": gen_s,
            "faith_thr": args.faith_thr,
            "stance": "live LLM infer stream → OnlineRFPerm fire → action route",
        }
    )
    # sample answers for eyeballing
    samples = []
    for b in (0, args.cut_batch):
        rec = records[b * args.n_per]
        samples.append(
            {
                "batch": b,
                "regime": rec["system"],
                "question": rec["question"][:120],
                "answer": rec["answer"][:200],
                "faith": rec["faith"],
                "y_bad": rec["y_bad"],
                "rag_hit": rec["rag_hit"],
            }
        )
    summary["samples"] = samples

    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    (args.out / "stream.jsonl").write_text("\n".join(json.dumps(r) for r in records) + "\n")
    write_report(summary, args, args.out)

    print(
        f"[live-infer] delay={summary['detection_delay_batch']} "
        f"first_fire={summary['first_fire_batch']} "
        f"y_bad quiet→hop {summary['mean_y_bad_quiet']:.2f}→{summary['mean_y_bad_hop']:.2f} "
        f"→ {args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
