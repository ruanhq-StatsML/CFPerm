#!/usr/bin/env python3
"""Diff two review-agent cards (波次→波次展示差分，不改 Drill).

Helps auditors see what flipped between waves without re-reading full cards.
Display / paste only — 审出加速 scope.

  PYTHONPATH=. python3 scripts/tencent_gr/diff_review_agent_cards.py \\
    --prev results/tencent_gr_review_agent_card/review_agent_card.json \\
    --cur  results/tencent_gr_review_agent_cards_batch/agod_grad_rfperm/review_agent_card.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from export_review_agent_card import build_card, load_flags, DEFAULT_FLAGS  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def _load_card(path: Path, *, flags: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    blob = json.loads(path.read_text())
    if blob.get("card_type") == "review_agent_context":
        return blob
    # treat as localize summary.json
    return build_card(blob, source=str(path), flags=flags or load_flags(DEFAULT_FLAGS))


def diff_cards(prev: Dict[str, Any], cur: Dict[str, Any]) -> Dict[str, Any]:
    p_tips = list(prev.get("tips") or [])
    c_tips = list(cur.get("tips") or [])
    p_set, c_set = set(p_tips), set(c_tips)
    p_signs = dict((prev.get("direction") or {}).get("tip_signs") or {})
    c_signs = dict((cur.get("direction") or {}).get("tip_signs") or {})
    sign_flips = []
    for t in sorted(p_set & c_set):
        a, b = p_signs.get(t, "0"), c_signs.get(t, "0")
        if a != b:
            sign_flips.append({"feature": t, "prev": a, "cur": b})
    p_dy = (prev.get("direction") or {}).get("sign_Dy")
    c_dy = (cur.get("direction") or {}).get("sign_Dy")
    p_q = (prev.get("review_hint") or {}).get("queue_bucket")
    c_q = (cur.get("review_hint") or {}).get("queue_bucket")
    out: Dict[str, Any] = {
        "diff_type": "review_agent_card_wave_diff",
        "disclaimer": "展示差分，非定罪；不改 Drill 门。",
        "prev_source": prev.get("source_summary"),
        "cur_source": cur.get("source_summary"),
        "sign_Dy": {"prev": p_dy, "cur": c_dy, "flipped": p_dy != c_dy},
        "queue_bucket": {"prev": p_q, "cur": c_q, "flipped": p_q != c_q},
        "tips_added": sorted(c_set - p_set),
        "tips_removed": sorted(p_set - c_set),
        "tip_sign_flips": sign_flips,
        "paste_for_agent": "",
        "ticket_custom_fields": {
            "graph_shift_diff_sign_dy_prev": p_dy,
            "graph_shift_diff_sign_dy_cur": c_dy,
            "graph_shift_diff_queue_prev": p_q,
            "graph_shift_diff_queue_cur": c_q,
            "graph_shift_diff_tips_added": ",".join(sorted(c_set - p_set)[:5]),
            "graph_shift_diff_tips_removed": ",".join(sorted(p_set - c_set)[:5]),
            "graph_shift_disclaimer": "clue_not_conviction_diff",
        },
    }
    lines: List[str] = [
        "【审核上下文·波次差分】",
        f"sign_Dy: {p_dy} → {c_dy}" + (" (翻转)" if p_dy != c_dy else ""),
        f"队列: {p_q} → {c_q}" + (" (换队)" if p_q != c_q else ""),
    ]
    if out["tips_added"]:
        lines.append("新增 tip: " + ", ".join(out["tips_added"][:8]))
    if out["tips_removed"]:
        lines.append("消失 tip: " + ", ".join(out["tips_removed"][:8]))
    if sign_flips:
        lines.append(
            "符号翻转: "
            + ", ".join(f"{x['feature']}:{x['prev']}→{x['cur']}" for x in sign_flips[:8])
        )
    if not (out["tips_added"] or out["tips_removed"] or sign_flips or p_dy != c_dy or p_q != c_q):
        lines.append("无显著差分（展示层）")
    lines.append(f"声明: {out['disclaimer']}")
    out["paste_for_agent"] = "\n".join(lines)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Diff two review-agent cards (display only)")
    ap.add_argument("--prev", type=Path, required=True)
    ap.add_argument("--cur", type=Path, required=True)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_agent_card_diff",
    )
    ap.add_argument("--flags", type=Path, default=DEFAULT_FLAGS)
    args = ap.parse_args()
    flags = load_flags(args.flags)
    prev = _load_card(args.prev, flags=flags)
    cur = _load_card(args.cur, flags=flags)
    diff = diff_cards(prev, cur)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "review_card_diff.json").write_text(
        json.dumps(diff, indent=2, ensure_ascii=False) + "\n"
    )
    (args.out_dir / "review_card_diff.md").write_text(
        "# 审出卡波次差分\n\n```\n" + diff["paste_for_agent"] + "\n```\n"
    )
    (args.out_dir / "ticket_custom_fields.json").write_text(
        json.dumps(diff["ticket_custom_fields"], indent=2, ensure_ascii=False) + "\n"
    )
    print(diff["paste_for_agent"])
    print("wrote", args.out_dir)


if __name__ == "__main__":
    main()
