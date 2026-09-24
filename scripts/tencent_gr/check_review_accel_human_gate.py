#!/usr/bin/env python3
"""Human-use readiness gate for「审出加速包」.

Commercial stop-line: don't keep inventing B29+/ETA until a real auditor
stamps useful|not_useful (reviewer != smoke). Display/ops only; no Drill change.

  PYTHONPATH=. python3 scripts/tencent_gr/check_review_accel_human_gate.py \\
    --feedback docs/summaries/artifacts/review_feedback.jsonl
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
SMOKE_REVIEWERS = {"smoke", "bot", "ci", "test", "ood_timer", "agent", "cursor"}


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


def evaluate(rows: List[Dict[str, Any]], *, min_human: int = 1) -> Dict[str, Any]:
    human = [
        r
        for r in rows
        if r.get("label") in ("useful", "not_useful")
        and str(r.get("reviewer") or "").strip().lower() not in SMOKE_REVIEWERS
    ]
    smoke = [
        r
        for r in rows
        if str(r.get("reviewer") or "").strip().lower() in SMOKE_REVIEWERS
    ]
    ready = len(human) >= min_human
    next_action = (
        "解冻 ETA/规则 backlog 可议；继续只做审出卡现网接通"
        if ready
        else "把 review_agent_card.md 贴给 1 位真人审核，打 useful/not_useful（--reviewer 真名）"
    )
    return {
        "gate_type": "review_accel_human_use",
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "ready": ready,
        "status": "READY" if ready else "NOT_READY",
        "min_human": min_human,
        "n_human_feedback": len(human),
        "n_smoke_feedback": len(smoke),
        "n_total": len(rows),
        "human_reviewers": sorted(
            {str(r.get("reviewer")) for r in human if r.get("reviewer")}
        ),
        "feature_freeze": not ready,
        "next_action": next_action,
        "kpi_note": "件均审核耗时↓ 需真人有用率；未 READY 前禁止新开 ETA/规则包",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="审出加速包 human-use gate")
    ap.add_argument(
        "--feedback",
        type=Path,
        default=ROOT / "docs" / "summaries" / "artifacts" / "review_feedback.jsonl",
    )
    ap.add_argument("--min-human", type=int, default=1)
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "summaries" / "artifacts" / "review_accel_human_gate.json",
    )
    args = ap.parse_args()
    report = evaluate(load_rows(args.feedback), min_human=args.min_human)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print("wrote", args.out)
    raise SystemExit(0 if report["ready"] else 2)


if __name__ == "__main__":
    main()
