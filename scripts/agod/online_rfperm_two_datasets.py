#!/usr/bin/env python3
"""Two LLM-inference datasets × CPU generate × OnlineRFPerm → LaTeX.

Datasets
--------
  1) HaluEval QA   — knowledge-grounded hallucination / faithfulness
  2) SQuAD         — context-grounded extractive QA (classic infer)

Same blunt use-case on both: quiet (context in prompt) → hop (no context + invent)
→ OnlineRFPerm fire → action route. CPU transformers only.

Usage::

  PYTHONPATH=. python3 scripts/agod/online_rfperm_two_datasets.py
  PYTHONPATH=. python3 scripts/agod/online_rfperm_two_datasets.py --backend mock
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache"
HALU = CACHE / "halueval_qa_3000.jsonl"
SQUAD_CACHE = CACHE / "squad_val_500.jsonl"
OUT = ROOT / "results" / "agod" / "online_rfperm_two_datasets"
DOCS = ROOT / "docs" / "biz"

# Reuse the live-infer plumbing
sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402


def ensure_squad_cache(n: int = 500) -> Path:
    if SQUAD_CACHE.exists():
        return SQUAD_CACHE
    from datasets import load_dataset

    ds = load_dataset("rajpurkar/squad", split=f"validation[:{n}]")
    SQUAD_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with SQUAD_CACHE.open("w") as f:
        for row in ds:
            ans = row["answers"]["text"][0] if row["answers"]["text"] else ""
            f.write(
                json.dumps(
                    {
                        "question": row["question"],
                        "knowledge": row["context"],
                        "gold": ans,
                        "id": row["id"],
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )
    return SQUAD_CACHE


def load_jsonl_qa(path: Path, n: int) -> list[dict]:
    rows = []
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            rows.append(
                {
                    "question": str(r.get("question") or ""),
                    "knowledge": str(r.get("knowledge") or r.get("context") or ""),
                    "gold": str(r.get("gold") or r.get("answer") or ""),
                }
            )
            if len(rows) >= n:
                break
    if len(rows) < n:
        raise SystemExit(f"Need {n} rows from {path}, got {len(rows)}")
    return rows


def run_one(name: str, rows: list[dict], backend, args: argparse.Namespace) -> dict:
    print(f"\n=== {name}: generating n={len(rows)} ===", flush=True)
    t0 = time.time()
    records = live.run_generation(
        backend,
        rows,
        n_per=args.n_per,
        n_batches=args.n_batches,
        cut_batch=args.cut_batch,
        max_new_tokens=args.max_new_tokens,
        faith_thr=args.faith_thr,
    )
    gen_s = time.time() - t0
    summary = live.online_rfperm_on_stream(
        records, n_per=args.n_per, gate=args.gate, audit_k=args.audit_k
    )
    summary["dataset"] = name
    summary["gen_seconds"] = gen_s
    summary["backend"] = args.backend
    first = next((h for h in summary["hops"] if h["fired"]), None)
    summary["first_action"] = first["action"] if first else "—"
    # keep a few samples
    summary["samples"] = [
        {
            "batch": r["batch"],
            "regime": r["system"],
            "question": r["question"][:100],
            "answer": r["answer"][:160],
            "y_bad": r["y_bad"],
            "rag_hit": r["rag_hit"],
        }
        for r in (records[0], records[args.cut_batch * args.n_per])
    ]
    ds_out = args.out / name
    ds_out.mkdir(parents=True, exist_ok=True)
    (ds_out / "summary.json").write_text(json.dumps(summary, indent=2))
    (ds_out / "stream.jsonl").write_text("\n".join(json.dumps(r) for r in records) + "\n")
    print(
        f"=== {name}: delay={summary['detection_delay_batch']} "
        f"fire@{summary['first_fire_batch']} "
        f"y_bad {summary['mean_y_bad_quiet']:.2f}→{summary['mean_y_bad_hop']:.2f} ===",
        flush=True,
    )
    return summary


def latex_escape(s: str) -> str:
    return (
        str(s)
        .replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("_", "\\_")
        .replace("#", "\\#")
        .replace("{", "\\{")
        .replace("}", "\\}")
    )


def write_latex(summaries: list[dict], args: argparse.Namespace, out: Path) -> str:
    rows = []
    for s in summaries:
        delay = s["detection_delay_batch"]
        delay_s = "---" if delay is None else str(int(delay))
        fire = "---" if s["first_fire_batch"] is None else str(int(s["first_fire_batch"]))
        rows.append(
            f"{latex_escape(s['dataset'])} & {s['n']} & {s['cut_batch']} & {fire} & "
            f"{delay_s} & {s['mean_y_bad_quiet']:.2f}$\\to${s['mean_y_bad_hop']:.2f} & "
            f"{s['mean_rag_quiet']:.2f}$\\to${s['mean_rag_hop']:.2f} & "
            f"\\texttt{{{latex_escape(s['first_action'])}}} \\\\"
        )
    body_rows = "\n".join(rows)
    hop_blocks = []
    for s in summaries:
        lines = [
            f"\\paragraph{{{latex_escape(s['dataset'])}.}}"
            f" cut$={s['cut_batch']}$, first fire$={s['first_fire_batch']}$, "
            f"delay$={s['detection_delay_batch']}$."
        ]
        lines.append("\\begin{center}\\small")
        lines.append("\\begin{tabular}{rllrrrl}")
        lines.append(
            "\\toprule batch & regime & fired & ratio & "
            "$y_{\\mathrm{bad}}$ & rag & action \\\\ \\midrule"
        )
        for h in s["hops"]:
            ratio = "---" if h["ratio"] is None else f"{h['ratio']:.2f}"
            lines.append(
                f"{h['batch']} & {h['regime']} & {h['fired']} & {ratio} & "
                f"{h['mean_y_bad']:.2f} & {h['mean_rag_hit']:.2f} & "
                f"\\texttt{{{latex_escape(h['action'])}}} \\\\"
            )
        lines.append("\\bottomrule\\end{tabular}\\end{center}")
        hop_blocks.append("\n".join(lines))
    hops_tex = "\n\n".join(hop_blocks)
    backend = latex_escape(args.backend)
    halu_name = latex_escape(HALU.name)
    squad_name = latex_escape(SQUAD_CACHE.name)

    tex = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs,hyperref,amsmath}",
            r"\title{OnlineRFPerm on Live LLM Inference\\(Two Datasets, CPU)}",
            r"\author{CFPerm use-case note}",
            r"\date{\today}",
            r"\begin{document}",
            r"\maketitle",
            "",
            r"\paragraph{Setup.}",
            (
                f"Blunt serving use-case only: generate with a small CPU model "
                f"(\\texttt{{{backend}}}/\\texttt{{SmolLM2-135M-Instruct}}), "
                r"label faithfulness vs.\ context, run OnlineRFPerm, route "
                r"\texttt{model\_rollback\_or\_audit\_topk} / \texttt{retrieval\_refresh}. "
                f"Quiet batches keep context in the prompt; hop batches drop context and "
                f"force invent-style decoding. Gate $\\gamma={args.gate}$, "
                f"$n_{{\\mathrm{{per}}}}={args.n_per}$, $B={args.n_batches}$, "
                f"cut$={args.cut_batch}$."
            ),
            "",
            r"\paragraph{Datasets.}",
            (
                f"(i)~\\textbf{{HaluEval}} QA with knowledge snippets "
                f"(\\texttt{{{halu_name}}}); "
                f"(ii)~\\textbf{{SQuAD}} validation contexts "
                f"(\\texttt{{{squad_name}}})."
            ),
            "",
            r"\begin{table}[h]",
            r"\centering",
            r"\caption{OnlineRFPerm on CPU live infer (two datasets).}",
            r"\begin{tabular}{lrrrrlll}",
            r"\toprule",
            (
                r"dataset & $n$ & cut & fire & delay & "
                r"$y_{\mathrm{bad}}$ & rag\_hit & first action \\"
            ),
            r"\midrule",
            body_rows,
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
            hops_tex,
            "",
            r"\paragraph{Takeaway.}",
            (
                r"Same wiring on both sets: \emph{live generate $\to$ OnlineRFPerm fire "
                r"$\to$ action}. No new LLM method---only a timestamped quality hop on a "
                r"real inference stream."
            ),
            "",
            r"\end{document}",
            "",
        ]
    )
    (out / "two_datasets.tex").write_text(tex)
    (DOCS / "ONLINERFPERM_TWO_DATASETS.tex").write_text(tex)
    md = (
        "# OnlineRFPerm × two LLM infer datasets (CPU)\n\n"
        "LaTeX: `docs/biz/ONLINERFPERM_TWO_DATASETS.tex` "
        "· `results/agod/online_rfperm_two_datasets/two_datasets.tex`\n\n"
        f"| dataset | delay | fire | y_bad | action |\n|---|---|---|---|---|\n"
        + "\n".join(
            f"| {s['dataset']} | {s['detection_delay_batch']} | {s['first_fire_batch']} | "
            f"{s['mean_y_bad_quiet']:.2f}→{s['mean_y_bad_hop']:.2f} | `{s['first_action']}` |"
            for s in summaries
        )
        + "\n\n```bash\nPYTHONPATH=. python3 scripts/agod/online_rfperm_two_datasets.py\n```\n"
    )
    (out / "REPORT.md").write_text(md)
    (DOCS / "ONLINERFPERM_TWO_DATASETS.md").write_text(md)
    return tex


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", choices=["transformers", "mock", "openai", "vllm"], default="transformers")
    p.add_argument("--model", default=None)
    p.add_argument("--base-url", default=None)
    p.add_argument("--api-key", default=None)
    p.add_argument("--n-per", type=int, default=12)
    p.add_argument("--n-batches", type=int, default=6)
    p.add_argument("--cut-batch", type=int, default=3)
    p.add_argument("--max-new-tokens", type=int, default=40)
    p.add_argument(
        "--faith-thr",
        type=float,
        default=0.45,
        help="answer-precision threshold for binary Y",
    )
    p.add_argument("--gate", type=float, default=1.25)
    p.add_argument("--audit-k", type=int, default=5)
    p.add_argument("--out", type=Path, default=OUT)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    assert 0 < args.cut_batch < args.n_batches
    n = args.n_per * args.n_batches
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    ensure_squad_cache(500)
    backend = live.build_backend(args)

    specs = [
        ("halueval", load_jsonl_qa(HALU, n)),
        ("squad", load_jsonl_qa(SQUAD_CACHE, n)),
    ]
    summaries = [run_one(name, rows, backend, args) for name, rows in specs]
    payload = {
        "stance": "two LLM infer datasets, CPU generate, OnlineRFPerm fire→route",
        "backend": args.backend,
        "n_per": args.n_per,
        "n_batches": args.n_batches,
        "cut_batch": args.cut_batch,
        "datasets": summaries,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2))
    tex = write_latex(summaries, args, args.out)
    print("\n===== LaTeX =====\n")
    print(tex)
    print(f"\n[wrote] {args.out / 'two_datasets.tex'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
