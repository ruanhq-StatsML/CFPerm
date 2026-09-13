#!/usr/bin/env python3
"""Dump tensor shapes and one inspectable prototype hop.

  python3 scripts/show_prototype_board.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_multimodal_attribution import write_json
from prototype_board import (
    inspect_shapes,
    plot_prototype_board,
    prototype_hop,
    write_prototype_npz,
)
from typed_shift_stepsize import make_typed_stream

OUT = ROOT / "results" / "prototype_board"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--head", default="video")
    ap.add_argument("--t", type=int, default=3)
    args = ap.parse_args()
    stream = make_typed_stream(
        n_batches=8,
        n_per=48,
        n_classes=6,
        seed=2026,
        cov={"video": 0.18, "audio": 0.05, "text": 0.0},
    )
    hop = prototype_hop(stream, t=args.t)
    shapes = inspect_shapes(stream, t=args.t)
    OUT.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "shapes.json", shapes)
    write_prototype_npz(hop, OUT / "prototypes_video.npz", head=args.head)
    plot_prototype_board(hop, shapes, OUT / "prototype_board.png", head=args.head)
    b = shapes["blocks"][args.head]
    print("layout", shapes["layout"]["s"], flush=True)
    print("stream X", shapes["stream"]["X"], "classes", shapes["stream"]["n_classes"], flush=True)
    print(
        "%s  X=%s  μ=%s  V=%s  logits=%s  M=%s  cos_stale=%.3f  cos_sdc=%.3f"
        % (
            args.head,
            b["X_curr"],
            b["mu_old"],
            b["instdisc"]["V"],
            b["instdisc"]["logits"],
            b["gpm_M"],
            b["cos_stale"],
            b["cos_sdc"],
        ),
        flush=True,
    )
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
