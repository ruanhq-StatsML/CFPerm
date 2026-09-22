#!/usr/bin/env python3
"""Multi-step LLM audit: one row per assistant hop.

A conversation is a sequence of product hops, not one Y for the whole thread.

  observe the current pack
    → write this turn's reply
    → auditor writes Y = pass/fail for THIS hop

HH ``chosen`` is only the served trajectory. chosen/rejected is never Y.
Raw prompts are not stored. X is serving geometry of this hop
(the usual 13 reply features plus ``x_step``).

Usage::

    PYTHONPATH=. python3 scripts/build_llm_audit_multistep_xy.py
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.build_llm_audit_xy import CACHE, OUT, X_DIMS, download, style_vector
from scripts.llm_audit_online_bootstrap_prototype import (
    CUT_BATCH,
    FLIP_RATE,
    SEED,
    apply_preference_hop,
)

HH_TEST = "https://huggingface.co/datasets/Anthropic/hh-rlhf/resolve/refs%2Fconvert%2Fparquet/default/test/0000.parquet"

STEP_DIMS = [*X_DIMS, "step"]
STEP_COLS = [f"x_{d}" for d in STEP_DIMS]
_TURN = re.compile(r"\n\n(?=Human:|Assistant:)")
# Frozen policy pack: ship if this hop has enough content and is not shouty.
AUDITOR_CUT = 0.3
N_HOPS = 1200
N_PER = 80
MIN_TURNS = 2


def assistant_turns(conversation: str) -> list[str]:
    """Split a served HH conversation into assistant hops. Drop empty turns."""
    text = (conversation or "").strip()
    if not text:
        return []
    out = []
    for chunk in _TURN.split(text):
        chunk = chunk.strip()
        if not chunk.startswith("Assistant:"):
            continue
        body = chunk[len("Assistant:") :].strip()
        if body:
            out.append(body)
    return out


def prototype_auditor(text: str) -> int:
    """Frozen policy pack on this hop. Not HH chosen."""
    v = style_vector(text)
    score = float(v[0] + v[1] - 1.5 * v[9] - 0.3 * v[4])
    return int(score > AUDITOR_CUT)


def hop_vector(text: str, step: int) -> np.ndarray:
    x_step = min(float(step) / 8.0, 1.0)
    return np.append(style_vector(text), x_step)


def explode_conversations(conversations: list[str]) -> list[dict]:
    """One row per assistant hop. Keep only multi-step episodes."""
    rows = []
    episode = 0
    for conv in conversations:
        turns = assistant_turns(conv)
        if len(turns) < MIN_TURNS:
            continue
        for step, body in enumerate(turns):
            x = hop_vector(body, step)
            rows.append(
                {
                    "y": prototype_auditor(body),
                    "episode": episode,
                    "step": step,
                    **{c: float(x[i]) for i, c in enumerate(STEP_COLS)},
                }
            )
        episode += 1
    return rows


def pack_hops(rows: list[dict], n: int = N_HOPS, n_per: int = N_PER) -> list[dict]:
    n_use = (min(len(rows), n) // n_per) * n_per
    packed = []
    for i, row in enumerate(rows[:n_use]):
        rec = dict(row)
        rec["batch"] = i // n_per
        packed.append(rec)
    return packed


def write_multistep(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["y", "batch", "episode", "step", *STEP_COLS]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            out = {
                "y": int(row["y"]),
                "batch": int(row["batch"]),
                "episode": int(row["episode"]),
                "step": int(row["step"]),
            }
            for c in STEP_COLS:
                out[c] = f"{float(row[c]):.6g}"
            w.writerow(out)


def load_hh_test(n_hops: int = N_HOPS) -> list[dict]:
    path = download(HH_TEST, CACHE / "hh_rlhf_test.parquet")
    tbl = pq.read_table(path, columns=["chosen"])
    rows = explode_conversations([str(c or "") for c in tbl.column("chosen").to_pylist()])
    return pack_hops(rows, n=n_hops)


def overlay_hop(rows: list[dict], *, cut_batch: int = CUT_BATCH, seed: int = SEED) -> list[dict]:
    y = np.asarray([int(r["y"]) for r in rows], dtype=int)
    batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
    y_h, _ = apply_preference_hop(y, batch, cut_batch=cut_batch, flip_rate=FLIP_RATE, seed=seed)
    out = []
    for rec, yi in zip(rows, y_h):
        hopped = dict(rec)
        hopped["y"] = int(yi)
        out.append(hopped)
    return out


def summarize(rows: list[dict]) -> dict:
    y = np.asarray([int(r["y"]) for r in rows], dtype=int)
    steps = np.asarray([int(r["step"]) for r in rows], dtype=int)
    episodes = {int(r["episode"]) for r in rows}
    return {
        "n": int(len(rows)),
        "n_episodes": int(len(episodes)),
        "n_batches": int(max(int(r["batch"]) for r in rows) + 1) if rows else 0,
        "y_pass_rate": float(np.mean(y)) if len(y) else None,
        "mean_step": float(np.mean(steps)) if len(steps) else None,
        "frac_step_ge_1": float(np.mean(steps >= 1)) if len(steps) else None,
        "schema": ["y", "batch", "episode", "step", *STEP_COLS],
        "y_meaning": "1=pass on this hop, 0=fail. Not HH chosen.",
        "y_is_chosen": False,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    quiet = load_hh_test()
    hopped = overlay_hop(quiet)
    q_path = OUT / "xy_hh_multistep_consistent.csv"
    h_path = OUT / "xy_hh_multistep_hop.csv"
    write_multistep(q_path, quiet)
    write_multistep(h_path, hopped)
    info = {
        "note": "Multi-step audit: one row per assistant hop. HH chosen is not Y. Raw text is not stored.",
        "auditor": f"frozen pack score = n_toks + n_chars - 1.5*upper - 0.3*bang; pass iff score > {AUDITOR_CUT}",
        "cite": "bai2022hh-rlhf",
        "files": {
            q_path.name: {
                **summarize(quiet),
                "path": str(q_path.relative_to(ROOT)),
                "regime": "quiet",
            },
            h_path.name: {
                **summarize(hopped),
                "path": str(h_path.relative_to(ROOT)),
                "regime": "hop",
            },
        },
    }
    man = OUT / "MULTISTEP.json"
    man.write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))
    print("wrote", q_path)
    print("wrote", h_path)


if __name__ == "__main__":
    main()
