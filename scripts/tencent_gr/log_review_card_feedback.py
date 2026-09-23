#!/usr/bin/env python3
"""Append useful / not_useful feedback for a review-agent card.

Closes the 审出加速 loop without a full UI: auditors (or review agents)
stamp whether the clue card helped. Feeds later ranking weights.

  PYTHONPATH=. python3 scripts/tencent_gr/log_review_card_feedback.py \\
    --card results/tencent_gr_review_agent_card/review_agent_card.json \\
    --label useful --note "末跳队列对上了"
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]


def build_feedback(
    card: Dict[str, Any],
    *,
    label: str,
    note: str,
    reviewer: str,
) -> Dict[str, Any]:
    d = card.get("direction") or {}
    return {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "label": label,
        "note": note,
        "reviewer": reviewer,
        "source_summary": card.get("source_summary"),
        "sign_Dy": d.get("sign_Dy"),
        "queue_bucket": (card.get("review_hint") or {}).get("queue_bucket"),
        "tips_top3": (card.get("tips") or [])[:3],
        "ticket_fields_echo": {
            k: (card.get("ticket_custom_fields") or {}).get(k)
            for k in (
                "graph_shift_sign_dy",
                "graph_shift_queue_bucket",
                "graph_shift_tip_top3",
            )
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Log useful/not_useful on review card")
    ap.add_argument(
        "--card",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_agent_card" / "review_agent_card.json",
    )
    ap.add_argument(
        "--label",
        choices=("useful", "not_useful"),
        required=True,
    )
    ap.add_argument("--note", type=str, default="")
    ap.add_argument("--reviewer", type=str, default="human_or_agent")
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_agent_card" / "review_feedback.jsonl",
    )
    args = ap.parse_args()
    card = json.loads(args.card.read_text())
    row = build_feedback(
        card, label=args.label, note=args.note, reviewer=args.reviewer
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row, ensure_ascii=False) + "\n")
    # also mirror under docs artifacts for the repo trail
    art = ROOT / "docs" / "summaries" / "artifacts" / "review_feedback.jsonl"
    with art.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row, ensure_ascii=False) + "\n")
    print(json.dumps(row, ensure_ascii=False, indent=2))
    print("appended", args.out)


if __name__ == "__main__":
    main()
