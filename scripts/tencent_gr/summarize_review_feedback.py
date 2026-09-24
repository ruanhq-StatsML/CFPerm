#!/usr/bin/env python3
"""Roll up review-card feedback into ops-facing useful rates.

  PYTHONPATH=. python3 scripts/tencent_gr/summarize_review_feedback.py \\
    --feedback docs/summaries/artifacts/review_feedback.jsonl
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]


def load_rows(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append(json.loads(line))
    return rows


def summarize(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    by_bucket: Dict[str, Dict[str, int]] = defaultdict(lambda: {"useful": 0, "not_useful": 0})
    total = {"useful": 0, "not_useful": 0}
    for r in rows:
        lab = r.get("label")
        if lab not in ("useful", "not_useful"):
            continue
        b = str(r.get("queue_bucket") or "unknown")
        by_bucket[b][lab] += 1
        total[lab] += 1

    def rate(u: int, n: int) -> float | None:
        s = u + n
        return round(u / s, 4) if s else None

    buckets = []
    for b, c in sorted(by_bucket.items(), key=lambda x: -(x[1]["useful"] + x[1]["not_useful"])):
        buckets.append(
            {
                "queue_bucket": b,
                "useful": c["useful"],
                "not_useful": c["not_useful"],
                "useful_rate": rate(c["useful"], c["not_useful"]),
            }
        )
    return {
        "n_feedback": len(rows),
        "overall": {
            "useful": total["useful"],
            "not_useful": total["not_useful"],
            "useful_rate": rate(total["useful"], total["not_useful"]),
        },
        "by_queue_bucket": buckets,
        "kpi_note": "人分钟代理：useful_rate↑ 且件均耗时可另接工单系统",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize review-card feedback")
    ap.add_argument(
        "--feedback",
        type=Path,
        default=ROOT / "docs" / "summaries" / "artifacts" / "review_feedback.jsonl",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "summaries" / "artifacts" / "review_feedback_summary.json",
    )
    args = ap.parse_args()
    summary = summarize(load_rows(args.feedback))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
