#!/usr/bin/env python3
"""Suggest gray sample_rate / kill_switch from review-card useful_rate.

System-interface only (constraint 1): feedback summary → flags suggestion.
Does NOT auto-apply; operator copies suggested JSON. No Drill changes.

  PYTHONPATH=. python3 scripts/tencent_gr/suggest_review_card_sample_rate.py \\
    --summary docs/summaries/artifacts/review_feedback_summary.json \\
    --flags configs/review_agent_card_flags.json
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]


def suggest_flags(
    feedback_summary: Dict[str, Any],
    current_flags: Dict[str, Any],
    *,
    min_n: int = 5,
    low_rate: float = 0.4,
    high_rate: float = 0.7,
    kill_rate: float = 0.2,
) -> Dict[str, Any]:
    overall = dict(feedback_summary.get("overall") or {})
    useful = int(overall.get("useful") or 0)
    not_useful = int(overall.get("not_useful") or 0)
    n = useful + not_useful
    rate = overall.get("useful_rate")
    cur_rate = float(current_flags.get("sample_rate", 1.0))
    cur_rate = max(0.0, min(1.0, cur_rate))
    suggested = {
        "enabled": bool(current_flags.get("enabled", True)),
        "sample_rate": cur_rate,
        "kill_switch": bool(current_flags.get("kill_switch", False)),
        "note": current_flags.get("note")
        or "审出加速灰度：由 useful_rate 建议，需人工确认后覆盖",
    }
    action = "hold"
    reason = "insufficient_feedback"
    if n < min_n or rate is None:
        reason = f"need_n>={min_n} (have {n})"
    elif float(rate) < kill_rate:
        suggested["kill_switch"] = True
        suggested["sample_rate"] = min(cur_rate, 0.1)
        action = "kill_switch_on"
        reason = f"useful_rate={rate}<{kill_rate}"
    elif float(rate) < low_rate:
        suggested["sample_rate"] = round(max(0.1, cur_rate * 0.5), 4)
        action = "decrease_sample_rate"
        reason = f"useful_rate={rate}<{low_rate}"
    elif float(rate) >= high_rate and cur_rate < 1.0:
        suggested["sample_rate"] = round(min(1.0, max(cur_rate * 1.25, cur_rate + 0.1)), 4)
        suggested["kill_switch"] = False
        action = "increase_sample_rate"
        reason = f"useful_rate={rate}>={high_rate}"
    else:
        action = "hold"
        reason = f"useful_rate={rate} in band"
    return {
        "suggest_type": "review_card_gray_from_useful_rate",
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "disclaimer": "建议文件，不自动覆盖现网 flags；人工确认后 cp。",
        "inputs": {
            "n_labeled": n,
            "useful_rate": rate,
            "current_sample_rate": cur_rate,
            "current_kill_switch": bool(current_flags.get("kill_switch", False)),
            "min_n": min_n,
        },
        "action": action,
        "reason": reason,
        "suggested_flags": suggested,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Suggest gray flags from useful_rate")
    ap.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "docs" / "summaries" / "artifacts" / "review_feedback_summary.json",
    )
    ap.add_argument(
        "--flags",
        type=Path,
        default=ROOT / "configs" / "review_agent_card_flags.json",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "configs" / "review_agent_card_flags.suggested.json",
    )
    ap.add_argument("--min-n", type=int, default=5)
    args = ap.parse_args()
    fb = json.loads(args.summary.read_text()) if args.summary.exists() else {}
    flags = (
        json.loads(args.flags.read_text())
        if args.flags.exists()
        else {"enabled": True, "sample_rate": 1.0, "kill_switch": False}
    )
    report = suggest_flags(fb, flags, min_n=args.min_n)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(report["suggested_flags"], indent=2, ensure_ascii=False) + "\n"
    )
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print("wrote", args.out, report_path)


if __name__ == "__main__":
    main()
