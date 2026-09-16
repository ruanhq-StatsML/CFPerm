#!/usr/bin/env python3
"""OnlineRFPerm live-infer on multiple QA datasets with n_ref >= 100.

Default clock: n_per=20, cut_batch=5 → reference = 100; n_batches=10 → trail=100.

  PYTHONPATH=. python3 scripts/agod/prepare_infer_datasets.py
  PYTHONPATH=. python3 scripts/agod/online_rfperm_multi_datasets.py --backend mock
  PYTHONPATH=. python3 scripts/agod/online_rfperm_multi_datasets.py --backend transformers
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]
ROOT = Path(__file__).resolve().parents[2]
EXPORT = ROOT / "data" / "hf_cache" / "infer_bench_export"
OUT = ROOT / "results" / "agod" / "online_rfperm_multi_datasets"
DOCS = ROOT / "docs" / "biz"

sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402

DATASETS = ("halueval", "squad", "hotpotqa", "truthfulqa")


def load_rows(name: str, n: int) -> list[dict]:
    path = EXPORT / f"{name}.jsonl"
    if not path.exists():
        raise SystemExit(f"missing {path}; run prepare_infer_datasets.py first")
    rows = []
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            rows.append(
                {
                    "question": str(r["question"]),
                    "knowledge": str(r.get("knowledge") or ""),
                    "gold": str(r.get("gold") or ""),
                }
            )
            if len(rows) >= n:
                break
    if len(rows) < n:
        raise SystemExit(f"{name}: need {n}, got {len(rows)}")
    return rows


def run_one(name: str, backend, args: argparse.Namespace) -> dict:
    n = args.n_per * args.n_batches
    assert args.n_per * args.cut_batch >= 100, "n_ref = n_per*cut_batch must be >= 100"
    rows = load_rows(name, n)
    print(f"\n=== {name}: n={n} n_ref={args.n_per * args.cut_batch} ===", flush=True)
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
    summary.update(
        {
            "dataset": name,
            "n_ref": args.n_per * args.cut_batch,
            "gen_seconds": gen_s,
            "backend": args.backend,
        }
    )
    ds_out = args.out / name
    ds_out.mkdir(parents=True, exist_ok=True)
    (ds_out / "summary.json").write_text(json.dumps(summary, indent=2))
    (ds_out / "stream.jsonl").write_text("\n".join(json.dumps(r) for r in records) + "\n")
    # tidy DataFrame export (shared form)
    try:
        import pandas as pd

        df_dir = ROOT / "results" / "agod" / "infer_dataframes"
        df_dir.mkdir(parents=True, exist_ok=True)
        rows_df = []
        for r in records:
            rows_df.append(
                {
                    "t_idx": int(r["t"]),
                    "batch": int(r["batch"]),
                    "dataset": name,
                    "question": r["question"],
                    "answer": r["answer"],
                    "y": int(r["y_bad"]),
                    "hopped": bool(r["hopped"]),
                    "system": r["system"],
                    "rag_hit": float(r["rag_hit"]),
                    "faith": float(r["faith"]),
                    "faith_metric": r.get("faith_metric", "answer_precision"),
                    "knowledge": r.get("knowledge", ""),
                    "gold": r.get("gold", ""),
                }
            )
        df = pd.DataFrame(rows_df)
        df.to_parquet(df_dir / f"{name}_stream_table.parquet", index=False)
        df.head(20).to_csv(df_dir / f"{name}_stream_preview.csv", index=False)
        src = pd.DataFrame(
            [
                {
                    "t_idx": i,
                    "dataset": name,
                    "question": row["question"],
                    "knowledge": row["knowledge"],
                    "gold": row["gold"],
                }
                for i, row in enumerate(rows)
            ]
        )
        src.to_parquet(df_dir / f"{name}_source.parquet", index=False)
        src.head(20).to_csv(df_dir / f"{name}_source_preview.csv", index=False)
    except Exception as e:  # pragma: no cover
        print(f"[warn] dataframe export failed for {name}: {e}", flush=True)
    print(
        f"=== {name}: delay={summary['detection_delay_batch']} "
        f"fire@{summary['first_fire_batch']} n_ref={summary['n_ref']} "
        f"y {summary['mean_y_bad_quiet']:.2f}→{summary['mean_y_bad_hop']:.2f} "
        f"faith_metric=answer_precision ===",
        flush=True,
    )
    return summary


def write_latex(summaries: list[dict], args: argparse.Namespace, out: Path) -> str:
    rows = []
    for s in summaries:
        d = "---" if s["detection_delay_batch"] is None else str(int(s["detection_delay_batch"]))
        f = "---" if s["first_fire_batch"] is None else str(int(s["first_fire_batch"]))
        rows.append(
            f"{s['dataset']} & {s['n']} & {s['n_ref']} & {s['cut_batch']} & {f} & {d} & "
            f"{s['mean_y_bad_quiet']:.2f}$\\to${s['mean_y_bad_hop']:.2f} \\\\"
        )
    body = "\n".join(rows)
    tex = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs,amsmath}",
            r"\title{OnlineRFPerm multi-dataset LLM infer\\($n_{\mathrm{ref}}\ge 100$)}",
            r"\author{CFPerm}",
            r"\date{\today}",
            r"\begin{document}",
            r"\maketitle",
            r"\paragraph{Datasets.}",
            r"HaluEval, SQuAD, HotpotQA, TruthfulQA. "
            f"Clock: $n_{{\\mathrm{{per}}}}={args.n_per}$, cut$={args.cut_batch}$ "
            f"$\\Rightarrow n_{{\\mathrm{{ref}}}}={args.n_per * args.cut_batch}$.",
            r"\begin{table}[h]\centering",
            r"\caption{OnlineRFPerm on four infer datasets.}",
            r"\begin{tabular}{lrrrrrl}",
            r"\toprule dataset & $n$ & $n_{\mathrm{ref}}$ & cut & fire & delay & $y_{\mathrm{bad}}$ \\",
            r"\midrule",
            body,
            r"\bottomrule\end{tabular}\end{table}",
            r"\end{document}",
            "",
        ]
    )
    (out / "multi_datasets.tex").write_text(tex)
    (DOCS / "ONLINERFPERM_MULTI_DATASETS.tex").write_text(tex)
    return tex


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--backend", choices=["transformers", "mock", "openai", "vllm"], default="mock")
    ap.add_argument("--model", default=None)
    ap.add_argument("--base-url", default=None)
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--n-per", type=int, default=20)
    ap.add_argument("--n-batches", type=int, default=10)
    ap.add_argument("--cut-batch", type=int, default=5)
    ap.add_argument("--max-new-tokens", type=int, default=40)
    ap.add_argument(
        "--faith-thr",
        type=float,
        default=0.45,
        help="answer-precision threshold for binary Y (default 0.45)",
    )
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--audit-k", type=int, default=5)
    ap.add_argument("--datasets", nargs="+", default=list(DATASETS))
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    n_ref = args.n_per * args.cut_batch
    if n_ref < 100:
        raise SystemExit(f"n_ref={n_ref} < 100; raise --n-per/--cut-batch")
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    backend = live.build_backend(args)
    summaries = [run_one(name, backend, args) for name in args.datasets]
    (args.out / "summary.json").write_text(
        json.dumps({"n_ref": n_ref, "backend": args.backend, "datasets": summaries}, indent=2)
    )
    tex = write_latex(summaries, args, args.out)
    md = [
        "# Multi-dataset OnlineRFPerm (n_ref >= 100)\n",
        f"n_ref={n_ref}, backend={args.backend}\n",
        "| dataset | n | n_ref | fire | delay | y_bad |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for s in summaries:
        md.append(
            f"| {s['dataset']} | {s['n']} | {s['n_ref']} | {s['first_fire_batch']} | "
            f"{s['detection_delay_batch']} | {s['mean_y_bad_quiet']:.2f}→{s['mean_y_bad_hop']:.2f} |"
        )
    md.append("\nDatasets: `data/hf_cache/infer_bench_export/`\n")
    (args.out / "REPORT.md").write_text("\n".join(md) + "\n")
    (DOCS / "ONLINERFPERM_MULTI_DATASETS.md").write_text("\n".join(md) + "\n")
    print(tex)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
