#!/usr/bin/env python3
"""POST (or dry-run) review-card ticket_custom_fields to a SOAR/工单 sink.

Constraint 2: cross the ticket system without touching Drill.
Default --dry-run: print payload only. Real POST needs --no-dry-run + --url.

  PYTHONPATH=. python3 scripts/tencent_gr/post_review_card_ticket_fields.py \\
    --card results/tencent_gr_review_agent_card/review_agent_card.json
"""
from __future__ import annotations

import argparse
import json
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]


def load_ticket_fields(path: Path) -> Dict[str, Any]:
    blob = json.loads(path.read_text())
    if "ticket_custom_fields" in blob:
        return dict(blob["ticket_custom_fields"] or {})
    # already a flat ticket fields JSON
    return dict(blob)


def build_payload(fields: Dict[str, Any], *, source: str) -> Dict[str, Any]:
    return {
        "event": "graph_shift_review_card",
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": source,
        "disclaimer": "clue_not_conviction",
        "custom_fields": fields,
    }


def post_json(url: str, payload: Dict[str, Any], *, timeout: float = 10.0) -> Dict[str, Any]:
    data = json.dumps(payload, ensure_ascii=False).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json; charset=utf-8"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8", errors="replace")
            return {
                "ok": True,
                "status": getattr(resp, "status", 200),
                "body": body[:2000],
            }
    except urllib.error.HTTPError as e:
        return {"ok": False, "status": e.code, "body": e.read().decode("utf-8", errors="replace")[:2000]}
    except Exception as e:  # noqa: BLE001 — surface sink errors to ops
        return {"ok": False, "status": None, "body": str(e)}


def main() -> None:
    ap = argparse.ArgumentParser(description="Dry-run / POST ticket_custom_fields")
    ap.add_argument(
        "--card",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_agent_card" / "review_agent_card.json",
        help="review_agent_card.json or ticket_custom_fields.json",
    )
    ap.add_argument("--url", type=str, default="", help="SOAR/工单 webhook URL")
    ap.add_argument(
        "--dry-run",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="default: print payload only (safe)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "summaries" / "artifacts" / "ticket_post_receipt.json",
    )
    args = ap.parse_args()
    fields = load_ticket_fields(args.card)
    payload = build_payload(fields, source=str(args.card))
    receipt: Dict[str, Any] = {
        "dry_run": bool(args.dry_run),
        "url": args.url or None,
        "payload": payload,
        "result": None,
    }
    if args.dry_run or not args.url:
        receipt["result"] = {
            "ok": True,
            "status": "dry_run",
            "body": "no POST; pass --no-dry-run --url https://...",
        }
    else:
        receipt["result"] = post_json(args.url, payload)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(receipt["result"], indent=2, ensure_ascii=False))
    print("fields", len(fields), "→", args.out)


if __name__ == "__main__":
    main()
