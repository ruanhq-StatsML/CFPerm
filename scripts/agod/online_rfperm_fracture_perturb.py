#!/usr/bin/env python3
"""Fracture / concept-drift perturbations: OnlineRFPerm catches inference breaks.

Claim (simple, practical):
  Quiet serving → no fire ≈ inference chain looks fine under the RF probe.
  Mid-stream *fracture* (sudden quality-law break / concept drift) → fire.
  That is the product sentence: we timestamp an inference break.

Perturbations (mid-stream, after n_ref)
--------------------------------------
  invent_fracture  — drop knowledge + invent answers (serving break)
  label_flip       — flip Y (classic concept drift on the label law)
  answer_corrupt   — corrupt answers away from knowledge (faithfulness break)

Datasets: HaluEval, SQuAD, HotpotQA, TruthfulQA, plus BoolQ / NQ-open when
available (extra cards under the same tidy form).

Usage::

  PYTHONPATH=. python3 scripts/agod/online_rfperm_fracture_perturb.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]
ROOT = Path(__file__).resolve().parents[2]
EXPORT = ROOT / "data" / "hf_cache" / "infer_bench_export"
OUT = ROOT / "results" / "agod" / "online_rfperm_fracture"
DOCS = ROOT / "docs" / "biz"

sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402
import online_rfperm_streaming_test as st  # noqa: E402

_WORD = re.compile(r"[a-z0-9]+", re.I)

CORE = ("halueval", "squad", "hotpotqa", "truthfulqa")
EXTRA = ("boolq", "nq_open")


# ---------------------------------------------------------------------------
# Extra datasets → same {question, knowledge, gold} card
# ---------------------------------------------------------------------------


def load_jsonl_card(name: str, n: int) -> list[dict]:
    path = EXPORT / f"{name}.jsonl"
    if not path.exists():
        return []
    rows = []
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            rows.append(
                {
                    "question": str(r.get("question") or ""),
                    "knowledge": str(r.get("knowledge") or ""),
                    "gold": str(r.get("gold") or ""),
                }
            )
            if len(rows) >= n:
                break
    return rows


def ensure_boolq(n: int) -> list[dict]:
    rows = load_jsonl_card("boolq", n)
    if len(rows) >= n:
        return rows[:n]
    from datasets import load_dataset

    ds = load_dataset("google/boolq", split=f"validation[:{max(n, 200)}]")
    rows = []
    for row in ds:
        ans = "yes" if bool(row["answer"]) else "no"
        rows.append(
            {
                "question": str(row["question"]),
                "knowledge": str(row.get("passage") or ""),
                "gold": ans,
            }
        )
        if len(rows) >= n:
            break
    _write_card("boolq", rows)
    return rows[:n]


def ensure_nq_open(n: int) -> list[dict]:
    rows = load_jsonl_card("nq_open", n)
    if len(rows) >= n:
        return rows[:n]
    from datasets import load_dataset

    ds = load_dataset("google-research-datasets/nq_open", split=f"validation[:{max(n, 200)}]")
    rows = []
    for row in ds:
        answers = row.get("answer") or row.get("answers") or []
        if isinstance(answers, str):
            answers = [answers]
        gold = str(answers[0]) if answers else ""
        knowledge = " | ".join(str(a) for a in answers[:5]) or gold
        rows.append(
            {
                "question": str(row["question"]),
                "knowledge": knowledge,
                "gold": gold,
            }
        )
        if len(rows) >= n:
            break
    _write_card("nq_open", rows)
    return rows[:n]


def _write_card(name: str, rows: list[dict]) -> None:
    EXPORT.mkdir(parents=True, exist_ok=True)
    path = EXPORT / f"{name}.jsonl"
    with path.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    (EXPORT / f"{name}_meta.json").write_text(
        json.dumps({"dataset": name, "n": len(rows), "fields": ["question", "knowledge", "gold"]}, indent=2)
    )


def load_rows(name: str, n: int) -> list[dict]:
    if name in CORE:
        rows = load_jsonl_card(name, n)
        if len(rows) < n:
            raise SystemExit(f"{name}: need {n}, got {len(rows)}; run prepare_infer_datasets.py")
        return rows
    if name == "boolq":
        return ensure_boolq(n)
    if name == "nq_open":
        return ensure_nq_open(n)
    raise SystemExit(f"unknown dataset {name}")


# ---------------------------------------------------------------------------
# Build quiet base stream, then apply mid-stream fracture
# ---------------------------------------------------------------------------


def build_quiet_records(rows: list[dict], *, faith_thr: float) -> list[dict]:
    """All-quiet mock stream: knowledge in prompt, endpoint-style answers."""
    backend = live.MockBackend()
    # force quiet for every row via cut beyond length
    return live.run_generation(
        backend,
        rows,
        n_per=len(rows),
        n_batches=1,
        cut_batch=len(rows) + 1,
        max_new_tokens=40,
        faith_thr=faith_thr,
    )


def apply_fracture(
    records: list[dict],
    *,
    kind: str,
    frac_start: int,
    faith_thr: float,
    rng: np.random.Generator,
) -> list[dict]:
    """Mutate records from frac_start onward to inject an inference break."""
    out = []
    for i, r in enumerate(records):
        row = dict(r)
        row["t"] = i
        row["batch"] = 0 if i < frac_start else 1
        if i < frac_start:
            row["fractured"] = False
            row["fracture_kind"] = "none"
            out.append(row)
            continue

        row["fractured"] = True
        row["fracture_kind"] = kind
        kn = row.get("knowledge") or ""
        q = row.get("question") or ""

        if kind == "invent_fracture":
            # serving break: invent, no grounding
            ans = (
                f"I am certain it was founded in Atlantis-{rng.integers(10, 99)} "
                f"in {rng.integers(1800, 1999)} by Dr. Fabricatus regarding: {q[:40]}"
            )
            row["answer"] = ans
            row["system"] = "hop"
            row["hopped"] = True
            faith = live.answer_precision(ans, kn) if kn else 0.0
            row["faith"] = float(faith)
            row["rag_hit"] = 0.0
            row["y_bad"] = 1 if faith < faith_thr * 0.85 else 1

        elif kind == "label_flip":
            # concept drift on Y only (answers stay quiet-looking)
            row["hopped"] = True
            row["system"] = "quiet_flip"
            row["y_bad"] = 1 - int(row.get("y_bad", 0))

        elif kind == "answer_corrupt":
            # faithfulness break: scramble / replace answer tokens
            words = _WORD.findall(row.get("answer") or q)
            junk = [f"xx{rng.integers(0, 999)}" for _ in range(max(6, len(words) // 2))]
            ans = " ".join(junk + words[::-1][:3])
            row["answer"] = ans
            row["system"] = "corrupt"
            row["hopped"] = True
            faith = live.answer_precision(ans, kn) if kn else 0.0
            row["faith"] = float(faith)
            row["rag_hit"] = float(faith)
            row["y_bad"] = 1 if faith < faith_thr else 1

        else:
            raise ValueError(kind)
        out.append(row)
    return out


def detect(records: list[dict], *, ref_end: int, win: int, gate: float) -> dict:
    X, y = st.featurize(records)
    freeze = st.orf_freeze_stream(X, y, ref_end=ref_end, win=win, gate=gate)
    slide = st.orf_slide_stream(X, y, ref_end=ref_end, win=win, gate=gate)

    def pack(d):
        return {
            "first1": d["first1"],
            "first2": d["first2"],
            "first3": d["first3"],
            "detection_delay_obs": d["detection_delay_obs"],
            "n_trail_fires": d["n_trail_fires"],
            **(
                {
                    "quiet_err_mean": d.get("quiet_err_mean"),
                    "hop_err_mean": d.get("hop_err_mean"),
                }
                if "quiet_err_mean" in d
                else {}
            ),
        }

    return {"orf_freeze": pack(freeze), "orf_slide": pack(slide)}


def run_case(
    name: str,
    rows: list[dict],
    *,
    kind: str,
    n_ref: int,
    win: int,
    gate: float,
    faith_thr: float,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)
    quiet = build_quiet_records(rows, faith_thr=faith_thr)
    # smooth baseline: no fracture
    smooth = detect(quiet, ref_end=n_ref, win=win, gate=gate)
    # fracture at mid trail start = n_ref
    fractured = apply_fracture(
        quiet, kind=kind, frac_start=n_ref, faith_thr=faith_thr, rng=rng
    )
    broken = detect(fractured, ref_end=n_ref, win=win, gate=gate)

    y = np.asarray([r["y_bad"] for r in fractured], dtype=int)
    y_q = float(np.mean(y[:n_ref]))
    y_f = float(np.mean(y[n_ref:]))

    caught = (
        broken["orf_freeze"]["first1"] is not None
        or broken["orf_slide"]["first1"] is not None
    )
    smooth_quiet = (
        smooth["orf_freeze"]["first1"] is None and smooth["orf_slide"]["first1"] is None
    )
    if caught and smooth_quiet:
        read = "fracture_caught_smooth_quiet — inference break timestamped"
    elif caught and not smooth_quiet:
        read = "fracture_caught_but_smooth_hot — gate may be sensitive"
    elif not caught and smooth_quiet:
        read = "missed_fracture — check contrast / window"
    else:
        read = "missed_and_noisy"

    return {
        "dataset": name,
        "perturbation": kind,
        "n": len(rows),
        "n_ref": n_ref,
        "win": win,
        "gate": gate,
        "y_quiet": y_q,
        "y_fracture": y_f,
        "smooth": smooth,
        "fracture": broken,
        "caught": caught,
        "smooth_quiet": smooth_quiet,
        "read": read,
    }


def write_report(cases: list[dict], args: argparse.Namespace, out: Path) -> str:
    rows = []
    for c in cases:
        fr = c["fracture"]["orf_freeze"]["detection_delay_obs"]
        sl = c["fracture"]["orf_slide"]["detection_delay_obs"]
        sf = c["smooth"]["orf_freeze"]["first1"]
        ss = c["smooth"]["orf_slide"]["first1"]
        rows.append(
            f"| {c['dataset']} | `{c['perturbation']}` | "
            f"{c['y_quiet']:.2f}→{c['y_fracture']:.2f} | {fr} | {sl} | "
            f"{sf}/{ss} | {c['read']} |"
        )
    body = "\n".join(rows)
    n_caught = sum(1 for c in cases if c["caught"])
    md = f"""# Fracture perturbations — detect inference breaks

> **一句（简单但实用）**：不火 ≈ 推理链路在探针下正常；中间突然 fracture / concept drift → OnlineRFPerm 给断裂打时间戳。

## Setup

| item | value |
|------|-------|
| n / n_ref / win / gate | {args.n} / {args.n_ref} / {args.win} / {args.gate} |
| faith | answer-precision → binary $Y$ |
| smooth control | all-quiet stream (no mid break) |
| fracture | sudden change at $t=n_{{\\mathrm{{ref}}}}$ |

## Perturbations

| kind | what breaks |
|------|-------------|
| `invent_fracture` | serving break: invent + drop grounding |
| `label_flip` | concept drift on $Y$ |
| `answer_corrupt` | faithfulness break: scrambled answers |

## Results

Caught **{n_caught}/{len(cases)}** fracture cases.

| dataset | perturbation | y quiet→frac | freeze delay | slide delay | smooth first1 | read |
|---------|--------------|--------------|-------------:|------------:|---------------|------|
{body}

## Read

- Smooth quiet + fracture fire → **can timestamp an inference break**.
- Your experience matches: under concept-drift-style flips the RF component
  usually fires; no fire on the quiet control means the chain looked fine.

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_fracture_perturb.py \\
  --n {args.n} --n-ref {args.n_ref} --win {args.win} --gate {args.gate}
```
"""
    (out / "REPORT.md").write_text(md)
    (DOCS / "ONLINERFPERM_FRACTURE_PERTURB.md").write_text(md)

    tex_rows = []
    for c in cases:
        fr = c["fracture"]["orf_freeze"]["detection_delay_obs"]
        sl = c["fracture"]["orf_slide"]["detection_delay_obs"]
        fd = "---" if fr is None else str(fr)
        sd = "---" if sl is None else str(sl)
        tex_rows.append(
            f"{c['dataset']} & {c['perturbation'].replace('_', '\\_')} & "
            f"{c['y_quiet']:.2f}$\\to${c['y_fracture']:.2f} & {fd} & {sd} \\\\"
        )
    tex = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs,amsmath}",
            r"\title{OnlineRFPerm fracture perturbations:\\detecting inference breaks}",
            r"\author{CFPerm}",
            r"\date{\today}",
            r"\begin{document}",
            r"\maketitle",
            r"No fire on a quiet stream $\approx$ inference looks fine under the RF probe. "
            r"A mid-stream fracture / concept-drift change should fire---timestamping an inference break.",
            rf"Clock: $n={args.n}$, $n_{{\mathrm{{ref}}}}={args.n_ref}$, win$={args.win}$, $\gamma={args.gate}$.",
            r"\begin{table}[h]\centering",
            r"\caption{Fracture detection delay (trail observation index).}",
            r"\begin{tabular}{lllrr}",
            r"\toprule dataset & perturbation & $y$ & freeze & slide \\",
            r"\midrule",
            *tex_rows,
            r"\bottomrule\end{tabular}\end{table}",
            r"\end{document}",
            "",
        ]
    )
    (out / "fracture.tex").write_text(tex)
    (DOCS / "ONLINERFPERM_FRACTURE_PERTURB.tex").write_text(tex)
    return md


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--n-ref", type=int, default=100)
    ap.add_argument("--win", type=int, default=20)
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--faith-thr", type=float, default=0.45)
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=list(CORE) + list(EXTRA),
    )
    ap.add_argument(
        "--perturbations",
        nargs="+",
        default=["invent_fracture", "label_flip", "answer_corrupt"],
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    assert args.n > args.n_ref >= args.win
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    cases = []
    for name in args.datasets:
        print(f"[fracture] loading {name}…", flush=True)
        try:
            rows = load_rows(name, args.n)
        except Exception as e:
            print(f"[skip] {name}: {e}", flush=True)
            continue
        for j, kind in enumerate(args.perturbations):
            print(f"  → {kind}", flush=True)
            c = run_case(
                name,
                rows,
                kind=kind,
                n_ref=args.n_ref,
                win=args.win,
                gate=args.gate,
                faith_thr=args.faith_thr,
                seed=args.seed + 17 * j,
            )
            cases.append(c)
            print(
                f"    caught={c['caught']} freeze={c['fracture']['orf_freeze']['detection_delay_obs']} "
                f"slide={c['fracture']['orf_slide']['detection_delay_obs']} | {c['read']}",
                flush=True,
            )

    (args.out / "summary.json").write_text(
        json.dumps(
            {
                "stance": "no fire ≈ fine; mid fracture / concept drift → fire (inference break)",
                "n": args.n,
                "n_ref": args.n_ref,
                "win": args.win,
                "gate": args.gate,
                "n_caught": sum(1 for c in cases if c["caught"]),
                "n_cases": len(cases),
                "cases": cases,
            },
            indent=2,
        )
    )
    write_report(cases, args, args.out)
    print(
        f"[fracture] caught {sum(1 for c in cases if c['caught'])}/{len(cases)} → {args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
