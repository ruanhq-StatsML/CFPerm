#!/usr/bin/env python3
"""Prototype: formulate Hotpot's (1 query × 10 wiki paras) tensor shape.

Raw distractor example
    question:  str
    titles:    (10,)
    docs:      (10,)          # paragraph text, not stored in xy
    gold:      2 titles among the 10

Scores (still length 10)
    sparse BM25:  (10,)
    dense cosine: (10,)
    RRF fusion:   (10,)

Pair table (what to monitor)
    Y: (10,)  1 iff this title is a supporting fact
    X: (10, d)  channel scores / ranks / overlap / graph seed, no gold flag
    batch: query index, repeated 10 times  → one arrival = one candidate pool

Query-level fused Y is a collapse of this tensor, not the native shape.

Usage::

    PYTHONPATH=. python3 scripts/prototype_hotpot_10para_shape.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_hybrid_retrieval_xy import (  # noqa: E402
    N_CAND,
    N_ROWS,
    OUT,
    load_hotpot,
    query_pool_tensors,
)

OUT.mkdir(parents=True, exist_ok=True)


def show_one(i: int, question, supporting, context) -> dict:
    ten = query_pool_tensors(question, context, supporting, fracture_dense=False)
    y = ten["y_pair"].tolist()
    rec = {
        "query_index": i,
        "shape": {
            "titles": [N_CAND],
            "sparse": list(ten["sparse"].shape),
            "dense": list(ten["dense"].shape),
            "rrf": list(ten["rrf"].shape),
            "y_pair": [N_CAND],
            "X_pair": [N_CAND, 9],
            "batch": "query_id repeated 10 times",
        },
        "n_cand": ten["n_cand"],
        "n_gold_in_pool": ten["n_gold"],
        "y_pair": y,
        "sparse": [round(float(x), 4) for x in ten["sparse"]],
        "dense": [round(float(x), 4) for x in ten["dense"]],
        "rrf": [round(float(x), 4) for x in ten["rrf"]],
        "note": "wiki paragraph text and question are not written",
    }
    print(f"query {i}: tensor (1 × {N_CAND})")
    print("  y_pair ", y, "  # 1 = supporting title in this slot")
    print("  sparse ", rec["sparse"])
    print("  dense  ", rec["dense"])
    print("  rrf    ", rec["rrf"])
    print("  X      ", rec["shape"]["X_pair"], "  batch =", i, "×", N_CAND)
    return rec


def main() -> int:
    questions, supporting, contexts = load_hotpot(max(3, N_ROWS))
    examples = [
        show_one(i, questions[i], supporting[i], contexts[i]) for i in range(3)
    ]
    dest = OUT / "SHAPE_EXAMPLE.json"
    dest.write_text(json.dumps(examples, indent=2) + "\n")
    print(f"wrote {dest}")
    print(
        "pair table: n = n_query × 10 ="
        f" {N_ROWS} × {N_CAND} = {N_ROWS * N_CAND}  → xy_hotpot_pairs.csv"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
