#!/usr/bin/env python3
"""Batch-export review-agent cards from a results tree (审出加速·波次批处理).

Walks --root for **/summary.json, writes one card dir per summary under --out-dir.
Does not change Drill gates. Gray flags apply per source path.

  PYTHONPATH=. python3 scripts/tencent_gr/export_review_agent_cards_batch.py \\
    --root results --out-dir results/tencent_gr_review_agent_cards_batch
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from export_review_agent_card import (  # noqa: E402
    DEFAULT_FLAGS,
    build_card,
    card_to_md,
    load_flags,
)

ROOT = Path(__file__).resolve().parents[2]


def _safe_slug(path: Path, root: Path) -> str:
    try:
        rel = path.parent.relative_to(root)
    except ValueError:
        rel = path.parent
    raw = str(rel)
    if raw in (".", ""):
        raw = path.parent.name or "card"
    slug = re.sub(r"[^a-zA-Z0-9._-]+", "_", raw).strip("._")
    return slug or "card"


def find_summaries(root: Path) -> List[Path]:
    return sorted(p for p in root.rglob("summary.json") if p.is_file())


def export_one(
    summary: Path,
    *,
    out_parent: Path,
    root: Path,
    flags: Dict[str, Any],
) -> Path:
    blob = json.loads(summary.read_text())
    card = build_card(blob, source=str(summary), flags=flags)
    dest = out_parent / _safe_slug(summary, root)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "review_agent_card.json").write_text(
        json.dumps(card, indent=2, ensure_ascii=False) + "\n"
    )
    (dest / "review_agent_card.md").write_text(card_to_md(card))
    (dest / "ticket_custom_fields.json").write_text(
        json.dumps(card["ticket_custom_fields"], indent=2, ensure_ascii=False) + "\n"
    )
    (dest / "source.txt").write_text(str(summary) + "\n")
    return dest


def main() -> None:
    ap = argparse.ArgumentParser(description="Batch export review-agent cards")
    ap.add_argument("--root", type=Path, default=ROOT / "results")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_agent_cards_batch",
    )
    ap.add_argument("--flags", type=Path, default=DEFAULT_FLAGS)
    ap.add_argument("--limit", type=int, default=0, help="0 = no limit")
    args = ap.parse_args()

    flags = load_flags(args.flags)
    summaries = find_summaries(args.root)
    if args.limit > 0:
        summaries = summaries[: args.limit]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    index: List[Dict[str, Any]] = []
    for s in summaries:
        dest = export_one(s, out_parent=args.out_dir, root=args.root, flags=flags)
        allow = json.loads((dest / "review_agent_card.json").read_text())["gray_flags"][
            "allow"
        ]
        index.append(
            {
                "source": str(s),
                "out_dir": str(dest),
                "gray_allow": allow,
            }
        )
        print(f"{'ok' if allow else 'gray_off'}\t{s}\t→\t{dest}")
    (args.out_dir / "index.json").write_text(
        json.dumps({"n": len(index), "cards": index}, indent=2, ensure_ascii=False)
        + "\n"
    )
    print(f"wrote {len(index)} cards → {args.out_dir}")


if __name__ == "__main__":
    main()
