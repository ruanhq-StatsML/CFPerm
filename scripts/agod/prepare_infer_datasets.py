#!/usr/bin/env python3
"""Materialize a few LLM-infer QA datasets for OnlineRFPerm (n_ref >= 100).

Writes unified {question, knowledge, gold} jsonl under data/hf_cache/infer_bench_export/.

Datasets
--------
  1) halueval   — local HaluEval QA
  2) squad      — SQuAD validation contexts
  3) hotpotqa   — HotpotQA distractor (context paragraphs)
  4) truthfulqa — TruthfulQA generation (correct answers as knowledge)

Usage::

  PYTHONPATH=. python3 scripts/agod/prepare_infer_datasets.py
  PYTHONPATH=. python3 scripts/agod/prepare_infer_datasets.py --n 400
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache"
OUT = CACHE / "infer_bench_export"
HALU = CACHE / "halueval_qa_3000.jsonl"
SQUAD = CACHE / "squad_val_500.jsonl"


def _write(name: str, rows: list[dict], out: Path) -> Path:
    path = out / f"{name}.jsonl"
    with path.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    meta = {
        "dataset": name,
        "n": len(rows),
        "fields": ["question", "knowledge", "gold"],
        "path": str(path),
        "note": "quiet=knowledge in prompt; hop=drop knowledge + invent",
    }
    (out / f"{name}_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[ok] {name}: n={len(rows)} → {path}", flush=True)
    return path


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
    return rows


def load_squad(n: int) -> list[dict]:
    if not SQUAD.exists():
        from datasets import load_dataset

        ds = load_dataset("rajpurkar/squad", split=f"validation[:{max(n, 500)}]")
        SQUAD.parent.mkdir(parents=True, exist_ok=True)
        with SQUAD.open("w") as f:
            for row in ds:
                ans = row["answers"]["text"][0] if row["answers"]["text"] else ""
                f.write(
                    json.dumps(
                        {
                            "question": row["question"],
                            "knowledge": row["context"],
                            "gold": ans,
                        },
                        ensure_ascii=False,
                    )
                    + "\n"
                )
    rows = []
    with SQUAD.open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            rows.append(
                {
                    "question": str(r["question"]),
                    "knowledge": str(r.get("knowledge") or r.get("context") or ""),
                    "gold": str(r.get("gold") or r.get("answer") or ""),
                }
            )
            if len(rows) >= n:
                break
    return rows


def load_hotpot(n: int) -> list[dict]:
    from datasets import load_dataset

    ds = load_dataset("hotpotqa/hotpot_qa", "distractor", split=f"validation[:{n}]")
    rows = []
    for row in ds:
        # context: {'title': [...], 'sentences': [[...], ...]}
        ctx = row["context"]
        titles = ctx.get("title") or []
        sents = ctx.get("sentences") or []
        chunks = []
        for t, ss in zip(titles, sents):
            para = " ".join(ss) if isinstance(ss, list) else str(ss)
            chunks.append(f"{t}: {para}")
        knowledge = "\n".join(chunks)[:2500]
        rows.append(
            {
                "question": str(row["question"]),
                "knowledge": knowledge,
                "gold": str(row["answer"]),
            }
        )
    return rows[:n]


def load_truthfulqa(n: int) -> list[dict]:
    from datasets import load_dataset

    ds = load_dataset("truthfulqa/truthful_qa", "generation", split="validation")
    rows = []
    for row in ds:
        corrects = row.get("correct_answers") or []
        if isinstance(corrects, str):
            corrects = [corrects]
        knowledge = " | ".join(str(c) for c in corrects[:5])
        rows.append(
            {
                "question": str(row["question"]),
                "knowledge": knowledge or str(row.get("best_answer") or ""),
                "gold": str(row.get("best_answer") or (corrects[0] if corrects else "")),
            }
        )
        if len(rows) >= n:
            break
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=300, help="rows per dataset (>= n_ref+trail)")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    assert args.n >= 200, "need room for n_ref>=100 + trail"
    args.out.mkdir(parents=True, exist_ok=True)

    catalog = {}
    for name, loader in [
        ("halueval", load_halu),
        ("squad", load_squad),
        ("hotpotqa", load_hotpot),
        ("truthfulqa", load_truthfulqa),
    ]:
        rows = loader(args.n)
        if len(rows) < args.n:
            print(f"[warn] {name}: only {len(rows)} < {args.n}", flush=True)
        path = _write(name, rows, args.out)
        catalog[name] = {"n": len(rows), "path": str(path)}

    (args.out / "catalog.json").write_text(
        json.dumps(
            {
                "n_ref_min": 100,
                "recommended": {"n_per": 20, "cut_batch": 5, "n_batches": 10},
                # ref = 20*5 = 100
                "datasets": catalog,
            },
            indent=2,
        )
    )
    (args.out / "README.md").write_text(
        "# Infer datasets (n_ref >= 100)\n\n"
        "| dataset | file | knowledge |\n"
        "|---|---|---|\n"
        "| halueval | `halueval.jsonl` | HaluEval knowledge |\n"
        "| squad | `squad.jsonl` | SQuAD context |\n"
        "| hotpotqa | `hotpotqa.jsonl` | Hotpot paragraphs |\n"
        "| truthfulqa | `truthfulqa.jsonl` | correct answers |\n\n"
        "Recommended clock: `n_per=20`, `cut_batch=5` → **n_ref=100**, "
        "`n_batches=10` → trail=100.\n\n"
        "```bash\n"
        "PYTHONPATH=. python3 scripts/agod/prepare_infer_datasets.py\n"
        "PYTHONPATH=. python3 scripts/agod/online_rfperm_multi_datasets.py "
        "--n-per 20 --cut-batch 5 --n-batches 10\n"
        "```\n"
    )
    print(json.dumps(catalog, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
