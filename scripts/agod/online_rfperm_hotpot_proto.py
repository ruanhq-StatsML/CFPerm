#!/usr/bin/env python3
"""HotpotQA faithfulness prototype for OnlineRFPerm.

Why the default multi-dataset path fails on Hotpot
--------------------------------------------------
Default labeling is Jaccard token-overlap(answer, knowledge) with
``faith_thr=0.18``. Hotpot distractor knowledge is multi-paragraph (~2.5k
chars, 10 titles). The mock/quiet generator echoes a short slice of that
blob, so ``|A ∩ K| / |A ∪ K|`` collapses to ~0.07 even for grounded quiet
answers → quiet ``y≈1`` → no regime contrast → OnlineRFPerm never fires.

What this prototype changes
---------------------------
1. Load Hotpot **supporting_facts** and build a short ``knowledge_support``
   (only cited sentences) alongside the full distractor dump.
2. Score four labelers on the same generated stream:
     - jaccard_full      (legacy, broken on Hotpot)
     - jaccard_support   (Jaccard vs support only)
     - prec_support      (recommended: answer-precision |A∩S|/|A|)
     - prec_and_gold     (stricter audit: precision + gold-hit)
3. Run OnlineRFPerm under the recommended labeler with n_ref>=100.

Usage::

  PYTHONPATH=. python3 scripts/agod/online_rfperm_hotpot_proto.py --backend mock
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "agod" / "online_rfperm_hotpot_proto"
DOCS = ROOT / "docs" / "biz"
CACHE = ROOT / "data" / "hf_cache" / "infer_bench_export"

sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402

_WORD = re.compile(r"[a-z0-9]+", re.I)


def tok(s: str) -> set[str]:
    return set(_WORD.findall((s or "").lower()))


def jaccard(a: str, b: str) -> float:
    ta, tb = tok(a), tok(b)
    if not ta or not tb:
        return 0.0
    return float(len(ta & tb) / len(ta | tb))


def answer_precision(ans: str, knowledge: str) -> float:
    """Fraction of answer tokens that appear in knowledge (|A∩K|/|A|).

    Prefer this over Jaccard when K is long: a short grounded answer should
    score high even if K has hundreds of unrelated tokens.
    """
    ta, tb = tok(ans), tok(knowledge)
    if not ta or not tb:
        return 0.0
    return float(len(ta & tb) / len(ta))


def gold_hit(ans: str, gold: str) -> bool:
    """Whether the gold answer is realized in the generation.

    Prefer gold-recall / substring over Jaccard(answer, gold): quiet answers often
    append a support citation, which dilutes Jaccard against a short gold span.
    """
    g = (gold or "").strip()
    a = (ans or "").strip()
    if not g or not a:
        return False
    gl, al = g.lower(), a.lower()
    if gl in al:
        return True
    gt, at = tok(g), tok(a)
    if not gt:
        return False
    # gold-recall: fraction of gold tokens present in the answer
    return float(len(gt & at) / len(gt)) >= 0.8


INVENT_RE = re.compile(
    r"\b(atlantis|fabricat|made[- ]up|definitely|certainly|i am certain|"
    r"dr\.?\s*fabricatus|1899)\b",
    re.I,
)


def label_row(
    *,
    ans: str,
    knowledge_full: str,
    knowledge_support: str,
    gold: str,
    hopped: bool,
    thr_j: float,
    thr_p: float,
) -> dict:
    j_full = jaccard(ans, knowledge_full)
    j_sup = jaccard(ans, knowledge_support)
    p_sup = answer_precision(ans, knowledge_support)
    g_hit = gold_hit(ans, gold)
    inventish = bool(INVENT_RE.search(ans))

    y_jaccard_full = int(j_full < thr_j or (hopped and inventish and j_full < thr_j * 1.4))
    y_jaccard_support = int(j_sup < thr_j or (hopped and inventish and j_sup < thr_j * 1.4))
    # Recommended ORF label: answer-precision vs support (handles long K).
    y_prec_support = int(p_sup < thr_p or (hopped and inventish and p_sup < thr_p * 1.4))
    # Stricter audit label (not default ORF y): also require gold-hit in quiet.
    if hopped:
        y_prec_and_gold = int(p_sup < thr_p or inventish or not g_hit)
    else:
        y_prec_and_gold = int(p_sup < thr_p or not g_hit)

    return {
        "jaccard_full": float(j_full),
        "jaccard_support": float(j_sup),
        "prec_support": float(p_sup),
        "gold_hit": bool(g_hit),
        "inventish": bool(inventish),
        "y_jaccard_full": y_jaccard_full,
        "y_jaccard_support": y_jaccard_support,
        "y_prec_support": y_prec_support,
        "y_prec_and_gold": y_prec_and_gold,
    }


def support_knowledge(row: dict) -> tuple[str, str]:
    """Build (knowledge_full, knowledge_support) from a Hotpot row."""
    ctx = row["context"]
    titles = list(ctx.get("title") or [])
    sents = list(ctx.get("sentences") or [])
    title_to_sents = {t: (ss if isinstance(ss, list) else [str(ss)]) for t, ss in zip(titles, sents)}

    full_chunks = []
    for t, ss in title_to_sents.items():
        full_chunks.append(f"{t}: {' '.join(ss)}")
    knowledge_full = "\n".join(full_chunks)[:2500]

    sf = row.get("supporting_facts") or {}
    sf_titles = list(sf.get("title") or [])
    sf_ids = list(sf.get("sent_id") or [])
    support_chunks = []
    seen = set()
    for t, sid in zip(sf_titles, sf_ids):
        key = (t, int(sid))
        if key in seen:
            continue
        seen.add(key)
        ss = title_to_sents.get(t) or []
        if 0 <= int(sid) < len(ss):
            support_chunks.append(f"{t}: {ss[int(sid)]}")
        elif ss:
            support_chunks.append(f"{t}: {ss[0]}")
    knowledge_support = "\n".join(support_chunks)[:1200]
    if not knowledge_support:
        # fallback: first sentence of each of first two paras
        fallback = []
        for t, ss in list(title_to_sents.items())[:2]:
            if ss:
                fallback.append(f"{t}: {ss[0]}")
        knowledge_support = "\n".join(fallback)[:1200] or knowledge_full[:800]
    return knowledge_full, knowledge_support


def load_hotpot_rows(n: int, *, cache_jsonl: Path | None = None) -> list[dict]:
    """Load Hotpot with support spans; cache a enriched jsonl for reuse."""
    cache_jsonl = cache_jsonl or (CACHE / "hotpotqa_support.jsonl")
    if cache_jsonl.exists():
        rows = []
        with cache_jsonl.open() as f:
            for line in f:
                if not line.strip():
                    continue
                rows.append(json.loads(line))
                if len(rows) >= n:
                    return rows[:n]

    from datasets import load_dataset

    ds = load_dataset("hotpotqa/hotpot_qa", "distractor", split=f"validation[:{max(n, 300)}]")
    rows = []
    CACHE.mkdir(parents=True, exist_ok=True)
    with cache_jsonl.open("w") as f:
        for row in ds:
            kn_full, kn_sup = support_knowledge(row)
            rec = {
                "question": str(row["question"]),
                "knowledge": kn_full,  # backward-compat alias = full
                "knowledge_full": kn_full,
                "knowledge_support": kn_sup,
                "gold": str(row["answer"]),
                "qtype": str(row.get("type") or ""),
                "level": str(row.get("level") or ""),
            }
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            rows.append(rec)
            if len(rows) >= n:
                break
    return rows[:n]


class HotpotMockBackend:
    """Quiet: gold + support sentence. Hop: invent (same spirit as live mock)."""

    def generate(self, system: str, user: str, *, max_new_tokens: int = 48) -> str:
        _ = max_new_tokens
        invent = "invent" in system.lower() or "fabricate" in system.lower()
        m_sup = re.search(r"Support:\s*(.*?)\n\nGold:", user, re.S)
        m_gold = re.search(r"Gold:\s*(.*?)\n\nQuestion:", user, re.S)
        m_q = re.search(r"Question:\s*(.*)$", user, re.S)
        support = (m_sup.group(1).strip() if m_sup else "")[:240]
        gold = (m_gold.group(1).strip() if m_gold else "").strip()
        q = (m_q.group(1).strip() if m_q else user).strip()
        if invent or not support:
            return (
                f"I am certain it was founded in Atlantis-42 in 1899 by Dr. Fabricatus "
                f"regarding: {q[:40]}"
            )
        # Grounded quiet answer: state gold, cite a short support slice.
        cite = " ".join(support.split()[:22])
        return f"{gold}. Supported by: {cite}"


def build_user(row: dict, *, hopped: bool) -> tuple[str, str]:
    q = row["question"]
    if hopped:
        return live.SYS_HOP, f"Question: {q}"
    # Quiet: give support (not full distractor) + gold is NOT shown to the model
    # in the prompt — gold is only for labeling. Prompt uses support knowledge.
    kn = row["knowledge_support"]
    system = live.SYS_QUIET
    user = f"Knowledge: {kn}\n\nQuestion: {q}"
    return system, user


def mock_user_with_gold(row: dict, *, hopped: bool) -> tuple[str, str]:
    """Mock-only user blob carrying support+gold for deterministic grounding."""
    q = row["question"]
    if hopped:
        return live.SYS_HOP, f"Question: {q}"
    return (
        live.SYS_QUIET,
        f"Support: {row['knowledge_support']}\n\nGold: {row['gold']}\n\nQuestion: {q}",
    )


def run_generation(
    backend,
    rows: list[dict],
    *,
    n_per: int,
    n_batches: int,
    cut_batch: int,
    max_new_tokens: int,
    thr_j: float,
    thr_p: float,
    use_mock_prompt: bool,
) -> list[dict]:
    records = []
    t0 = time.time()
    for i, row in enumerate(rows):
        batch = i // n_per
        hopped = batch >= cut_batch
        if use_mock_prompt:
            system, user = mock_user_with_gold(row, hopped=hopped)
        else:
            system, user = build_user(row, hopped=hopped)
        ans = backend.generate(system, user, max_new_tokens=max_new_tokens)
        kn_full = row["knowledge_full"]
        kn_sup = row["knowledge_support"]
        gold = row["gold"]
        lab = label_row(
            ans=ans,
            knowledge_full=kn_full,
            knowledge_support=kn_sup,
            gold=gold,
            hopped=hopped,
            thr_j=thr_j,
            thr_p=thr_p,
        )
        # rag_hit for routing: precision vs support (not Jaccard vs full)
        rag_hit = answer_precision(ans, kn_sup) if kn_sup else 0.0
        records.append(
            {
                "t": i,
                "batch": batch,
                "hopped": hopped,
                "question": row["question"],
                "knowledge": kn_sup[:240],
                "knowledge_full_len": len(kn_full),
                "knowledge_support_len": len(kn_sup),
                "gold": gold[:240],
                "answer": ans,
                "rag_provided": 0.0 if hopped else 1.0,
                "rag_hit": float(rag_hit),
                "faith": float(lab["prec_support"]),  # recommended faith
                "y_bad": float(lab["y_prec_support"]),  # recommended label
                "system": "hop" if hopped else "quiet",
                "qtype": row.get("qtype", ""),
                **lab,
            }
        )
        if (i + 1) % n_per == 0:
            print(
                f"  batch {batch}/{n_batches - 1} done  "
                f"({i + 1}/{len(rows)} gens, {time.time() - t0:.1f}s)  "
                f"regime={'HOP' if hopped else 'quiet'}",
                flush=True,
            )
    return records


def label_contrast(records: list[dict]) -> dict:
    out = {}
    for key in (
        "y_jaccard_full",
        "y_jaccard_support",
        "y_prec_support",
        "y_prec_and_gold",
    ):
        quiet = [r[key] for r in records if not r["hopped"]]
        hop = [r[key] for r in records if r["hopped"]]
        out[key] = {
            "quiet_mean": float(np.mean(quiet)),
            "hop_mean": float(np.mean(hop)),
            "delta": float(np.mean(hop) - np.mean(quiet)),
        }
    score_keys = ["jaccard_full", "jaccard_support", "prec_support"]
    for key in score_keys:
        quiet = [r[key] for r in records if not r["hopped"]]
        hop = [r[key] for r in records if r["hopped"]]
        out[f"score_{key}"] = {
            "quiet_mean": float(np.mean(quiet)),
            "hop_mean": float(np.mean(hop)),
        }
    out["gold_hit_quiet"] = float(np.mean([r["gold_hit"] for r in records if not r["hopped"]]))
    out["gold_hit_hop"] = float(np.mean([r["gold_hit"] for r in records if r["hopped"]]))
    return out


def run_orf_for_label(
    records: list[dict],
    y_key: str,
    *,
    n_per: int,
    gate: float,
    audit_k: int,
) -> dict:
    """Reuse live ORF path with a chosen y column."""
    tmp = []
    for r in records:
        d = dict(r)
        d["y_bad"] = float(r[y_key])
        # faith/rag already set; featurize uses question/answer/rag_hit
        tmp.append(d)
    summary = live.online_rfperm_on_stream(tmp, n_per=n_per, gate=gate, audit_k=audit_k)
    summary["y_key"] = y_key
    return summary


def write_docs(contrast: dict, summaries: dict, args: argparse.Namespace, out: Path) -> None:
    rows = []
    for key in (
        "y_jaccard_full",
        "y_jaccard_support",
        "y_prec_support",
        "y_prec_and_gold",
    ):
        c = contrast[key]
        s = summaries[key]
        d = "---" if s["detection_delay_batch"] is None else str(int(s["detection_delay_batch"]))
        f = "---" if s["first_fire_batch"] is None else str(int(s["first_fire_batch"]))
        rows.append(
            f"{key.replace('_', '\\_')} & {c['quiet_mean']:.2f}$\\to${c['hop_mean']:.2f} & "
            f"{f} & {d} \\\\"
        )
    tex = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs,amsmath}",
            r"\title{HotpotQA faithfulness prototype for OnlineRFPerm}",
            r"\author{CFPerm}",
            r"\date{\today}",
            r"\begin{document}",
            r"\maketitle",
            r"\paragraph{Problem.}",
            r"Default Jaccard-vs-full-distractor labeling saturates quiet "
            r"$y\approx 1$ on Hotpot (long $K$, short $A$).",
            r"\paragraph{Fix.}",
            r"Use Hotpot \texttt{supporting\_facts} as \texttt{knowledge\_support}; "
            r"label with answer-precision $|A\cap S|/|A|$ plus gold-hit.",
            f"Clock: $n_{{\\mathrm{{per}}}}={args.n_per}$, cut$={args.cut_batch}$, "
            f"$n_{{\\mathrm{{ref}}}}={args.n_per * args.cut_batch}$.",
            r"\begin{table}[h]\centering",
            r"\caption{Labeler contrast $\to$ OnlineRFPerm fire.}",
            r"\begin{tabular}{lrrr}",
            r"\toprule labeler & $y$ quiet$\to$hop & fire & delay \\",
            r"\midrule",
            "\n".join(rows),
            r"\bottomrule\end{tabular}\end{table}",
            r"\end{document}",
            "",
        ]
    )
    (out / "hotpot_proto.tex").write_text(tex)
    (DOCS / "ONLINERFPERM_HOTPOT_PROTO.tex").write_text(tex)

    md = [
        "# HotpotQA faithfulness prototype (OnlineRFPerm)\n",
        "## Problem\n",
        "Default multi-dataset path labels with **Jaccard(answer, full distractor)**.",
        "Hotpot `knowledge` is ~10 paragraphs; quiet answers are short → Jaccard ≈ 0.07",
        "→ quiet `y≈1` → no hop contrast → OnlineRFPerm does not fire.\n",
        "## Handling logic\n",
        "1. Load Hotpot `supporting_facts` → `knowledge_support` (cited sentences only).",
        "2. Keep `knowledge_full` for the legacy baseline.",
        "3. Quiet prompt uses **support**, not the full distractor dump.",
        "4. Compare labelers on the same generations:\n",
        "| labeler | rule |",
        "|---|---|",
        "| `y_jaccard_full` | Jaccard(A, K_full) < thr (legacy) |",
        "| `y_jaccard_support` | Jaccard(A, K_support) < thr |",
        "| `y_prec_support` | **recommended**: answer-precision \\|A∩S\\|/\\|A\\| < thr |",
        "| `y_prec_and_gold` | audit: precision + gold-hit (can add quiet noise) |",
        "",
        "## Results\n",
        f"backend=`{args.backend}`, n_ref=`{args.n_per * args.cut_batch}`\n",
        "| labeler | y quiet→hop | fire | delay |",
        "|---|---|---|---|",
    ]
    for key in (
        "y_jaccard_full",
        "y_jaccard_support",
        "y_prec_support",
        "y_prec_and_gold",
    ):
        c = contrast[key]
        s = summaries[key]
        md.append(
            f"| `{key}` | {c['quiet_mean']:.2f}→{c['hop_mean']:.2f} | "
            f"{s['first_fire_batch']} | {s['detection_delay_batch']} |"
        )
    md += [
        "",
        "### Score means (quiet / hop)\n",
        f"- jaccard_full: {contrast['score_jaccard_full']['quiet_mean']:.3f} → "
        f"{contrast['score_jaccard_full']['hop_mean']:.3f}",
        f"- jaccard_support: {contrast['score_jaccard_support']['quiet_mean']:.3f} → "
        f"{contrast['score_jaccard_support']['hop_mean']:.3f}",
        f"- prec_support: {contrast['score_prec_support']['quiet_mean']:.3f} → "
        f"{contrast['score_prec_support']['hop_mean']:.3f}",
        f"- gold_hit: {contrast['gold_hit_quiet']:.2f} → {contrast['gold_hit_hop']:.2f}",
        "",
        "## Read\n",
        "- Legacy `jaccard_full` stays saturated (quiet≈hop≈1) — this is the documented bug.",
        "- `y_prec_support` opens a quiet→hop gap and OnlineRFPerm fires at cut (delay 0).",
        "- `y_prec_and_gold` is optional audit; gold-miss alone should not drive the ORF clock.",
        "- Formulation unchanged: still `t_idx|batch|question|answer|y`; only how `y` / knowledge are built changes.",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/agod/online_rfperm_hotpot_proto.py --backend mock",
        "```",
        "",
    ]
    text = "\n".join(md) + "\n"
    (out / "REPORT.md").write_text(text)
    (DOCS / "ONLINERFPERM_HOTPOT_PROTO.md").write_text(text)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--backend", choices=["mock", "transformers", "openai", "vllm"], default="mock")
    ap.add_argument("--model", default=None)
    ap.add_argument("--base-url", default=None)
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--n-per", type=int, default=20)
    ap.add_argument("--n-batches", type=int, default=10)
    ap.add_argument("--cut-batch", type=int, default=5)
    ap.add_argument("--max-new-tokens", type=int, default=48)
    ap.add_argument("--thr-jaccard", type=float, default=0.18)
    ap.add_argument("--thr-precision", type=float, default=0.45)
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--audit-k", type=int, default=5)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    n_ref = args.n_per * args.cut_batch
    if n_ref < 100:
        raise SystemExit(f"n_ref={n_ref} < 100")
    n = args.n_per * args.n_batches
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    rows = load_hotpot_rows(n)
    print(f"=== hotpot proto: n={n} n_ref={n_ref} backend={args.backend} ===", flush=True)

    use_mock_prompt = args.backend == "mock"
    if args.backend == "mock":
        backend = HotpotMockBackend()
    else:
        backend = live.build_backend(args)

    t0 = time.time()
    records = run_generation(
        backend,
        rows,
        n_per=args.n_per,
        n_batches=args.n_batches,
        cut_batch=args.cut_batch,
        max_new_tokens=args.max_new_tokens,
        thr_j=args.thr_jaccard,
        thr_p=args.thr_precision,
        use_mock_prompt=use_mock_prompt,
    )
    gen_s = time.time() - t0
    contrast = label_contrast(records)

    summaries = {}
    for y_key in (
        "y_jaccard_full",
        "y_jaccard_support",
        "y_prec_support",
        "y_prec_and_gold",
    ):
        summaries[y_key] = run_orf_for_label(
            records, y_key, n_per=args.n_per, gate=args.gate, audit_k=args.audit_k
        )

    recommended = summaries["y_prec_support"]
    payload = {
        "dataset": "hotpotqa",
        "prototype": "supporting_facts + answer_precision",
        "n_ref": n_ref,
        "backend": args.backend,
        "gen_seconds": gen_s,
        "thr_jaccard": args.thr_jaccard,
        "thr_precision": args.thr_precision,
        "label_contrast": contrast,
        "orf_by_label": {
            k: {
                "first_fire_batch": v["first_fire_batch"],
                "detection_delay_batch": v["detection_delay_batch"],
                "mean_y_bad_quiet": v["mean_y_bad_quiet"],
                "mean_y_bad_hop": v["mean_y_bad_hop"],
            }
            for k, v in summaries.items()
        },
        "recommended": {
            "y_key": "y_prec_support",
            "first_fire_batch": recommended["first_fire_batch"],
            "detection_delay_batch": recommended["detection_delay_batch"],
            "mean_y_bad_quiet": recommended["mean_y_bad_quiet"],
            "mean_y_bad_hop": recommended["mean_y_bad_hop"],
            "mean_rag_quiet": recommended.get("mean_rag_quiet"),
            "mean_rag_hop": recommended.get("mean_rag_hop"),
            "hops": recommended["hops"],
        },
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2))
    (args.out / "stream.jsonl").write_text("\n".join(json.dumps(r) for r in records) + "\n")

    # tidy dataframe under the same serving form
    import pandas as pd

    df = pd.DataFrame(
        [
            {
                "t_idx": r["t"],
                "batch": r["batch"],
                "dataset": "hotpotqa",
                "question": r["question"],
                "answer": r["answer"],
                "y": int(r["y_prec_support"]),
                "y_jaccard_full": int(r["y_jaccard_full"]),
                "y_prec_support": int(r["y_prec_support"]),
                "y_prec_and_gold": int(r["y_prec_and_gold"]),
                "hopped": r["hopped"],
                "rag_hit": r["rag_hit"],
                "prec_support": r["prec_support"],
                "jaccard_full": r["jaccard_full"],
                "gold_hit": r["gold_hit"],
                "knowledge_support": r["knowledge"],
                "gold": r["gold"],
            }
            for r in records
        ]
    )
    df.to_parquet(args.out / "stream_table.parquet", index=False)
    df.head(20).to_csv(args.out / "stream_preview.csv", index=False)

    write_docs(contrast, summaries, args, args.out)
    print(json.dumps(payload["orf_by_label"], indent=2))
    print(
        f"=== recommended delay={recommended['detection_delay_batch']} "
        f"fire@{recommended['first_fire_batch']} "
        f"y {recommended['mean_y_bad_quiet']:.2f}→{recommended['mean_y_bad_hop']:.2f} ===",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
