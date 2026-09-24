#!/usr/bin/env python3
"""End-to-end smoke for「审出加速包」(export→feedback→rollup→suggest→dry-run POST).

Constraint 5: reuses summary.json / direction only. No Drill changes.

  PYTHONPATH=. python3 scripts/tencent_gr/smoke_review_accel_pack.py \\
    --summary results/tencent_gr_w1w2_mmd_po_fsds/summary.json
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
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
    load_tip_overlay,
)
from log_review_card_feedback import build_feedback  # noqa: E402
from post_review_card_ticket_fields import build_payload  # noqa: E402
from suggest_review_card_sample_rate import suggest_flags  # noqa: E402
from summarize_review_feedback import summarize  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def run_smoke(
    summary_path: Path,
    *,
    out_dir: Path,
    flags_path: Path,
    tip_overlay_path: Path | None,
    label: str,
) -> Dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    blob = json.loads(summary_path.read_text())
    flags = load_flags(flags_path)
    overlay = load_tip_overlay(tip_overlay_path)
    card = build_card(
        blob, source=str(summary_path), flags=flags, tip_overlay=overlay
    )
    (out_dir / "review_agent_card.json").write_text(
        json.dumps(card, indent=2, ensure_ascii=False) + "\n"
    )
    (out_dir / "review_agent_card.md").write_text(card_to_md(card))
    (out_dir / "ticket_custom_fields.json").write_text(
        json.dumps(card["ticket_custom_fields"], indent=2, ensure_ascii=False) + "\n"
    )

    fb = build_feedback(card, label=label, note="smoke_review_accel_pack", reviewer="smoke")
    fb_path = out_dir / "review_feedback.jsonl"
    with fb_path.open("a") as f:
        f.write(json.dumps(fb, ensure_ascii=False) + "\n")
    rows: List[Dict[str, Any]] = []
    for line in fb_path.read_text().splitlines():
        if line.strip():
            rows.append(json.loads(line))
    fb_summary = summarize(rows)
    (out_dir / "review_feedback_summary.json").write_text(
        json.dumps(fb_summary, indent=2, ensure_ascii=False) + "\n"
    )

    suggest = suggest_flags(fb_summary, flags, min_n=5)
    (out_dir / "flags.suggested.json").write_text(
        json.dumps(suggest["suggested_flags"], indent=2, ensure_ascii=False) + "\n"
    )
    (out_dir / "flags.suggested.report.json").write_text(
        json.dumps(suggest, indent=2, ensure_ascii=False) + "\n"
    )

    payload = build_payload(card["ticket_custom_fields"], source=str(summary_path))
    receipt = {
        "dry_run": True,
        "url": None,
        "payload": payload,
        "result": {"ok": True, "status": "dry_run", "body": "smoke only"},
    }
    (out_dir / "ticket_post_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n"
    )

    report = {
        "smoke_type": "review_accel_pack",
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "summary": str(summary_path),
        "out_dir": str(out_dir),
        "sign_Dy": (card.get("direction") or {}).get("sign_Dy"),
        "queue_bucket": (card.get("review_hint") or {}).get("queue_bucket"),
        "sla_level": ((card.get("review_hint") or {}).get("sla_urgency") or {}).get(
            "level"
        ),
        "eta_soft_hint": (card.get("shared_context") or {}).get("eta_soft_hint"),
        "gray_allow": (card.get("gray_flags") or {}).get("allow"),
        "feedback_label": label,
        "suggest_action": suggest.get("action"),
        "ticket_field_count": len(card.get("ticket_custom_fields") or {}),
        "steps_ok": [
            "export_card",
            "log_feedback",
            "summarize_feedback",
            "suggest_flags",
            "dry_run_ticket_post",
        ],
    }
    (out_dir / "smoke_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    )
    md = [
        "# 审出加速包 Smoke",
        "",
        f"- ts: `{report['ts']}`",
        f"- summary: `{summary_path}`",
        f"- sign_Dy: **{report['sign_Dy']}**",
        f"- queue: **{report['queue_bucket']}**",
        f"- SLA: `{report['sla_level']}`",
        f"- eta_soft_hint: `{report['eta_soft_hint']}`",
        f"- gray_allow: `{report['gray_allow']}`",
        f"- feedback: `{label}`",
        f"- suggest_action: `{report['suggest_action']}`",
        f"- ticket_fields: {report['ticket_field_count']}",
        "",
        "Steps: " + " → ".join(report["steps_ok"]),
        "",
    ]
    (out_dir / "smoke_report.md").write_text("\n".join(md))
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description="Smoke 审出加速包 end-to-end")
    ap.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "results" / "tencent_gr_w1w2_mmd_po_fsds" / "summary.json",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_accel_smoke",
    )
    ap.add_argument("--flags", type=Path, default=DEFAULT_FLAGS)
    ap.add_argument("--tip-overlay", type=Path, default=None)
    ap.add_argument(
        "--label",
        choices=("useful", "not_useful"),
        default="useful",
        help="smoke feedback stamp",
    )
    args = ap.parse_args()
    report = run_smoke(
        args.summary,
        out_dir=args.out_dir,
        flags_path=args.flags,
        tip_overlay_path=args.tip_overlay,
        label=args.label,
    )
    print((args.out_dir / "smoke_report.md").read_text())
    print("wrote", args.out_dir)
    print(json.dumps({k: report[k] for k in ("sign_Dy", "sla_level", "suggest_action", "steps_ok")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
