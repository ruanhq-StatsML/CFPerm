#!/usr/bin/env python3
"""Export a review-agent context card from localize summary.json.

Commercial MVP for「审出加速包」: pasteable JSON/MD for audit agents.
Does not change Drill gates; does not claim fraud conviction.

  PYTHONPATH=. python3 scripts/tencent_gr/export_review_agent_card.py \\
    --summary results/tencent_gr_w1w2_mmd_po_fsds/summary.json \\
    --out-dir results/tencent_gr_review_agent_card
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FLAGS = ROOT / "configs" / "review_agent_card_flags.json"

TIP_BUCKETS = {
    "u_span_sec": "活跃跨度异常（短刷/长挂）",
    "i_credit_last": "末次归因偏高（末跳嫌疑）",
    "i_share_last": "末次份额偏高（末跳嫌疑）",
    "i_credit_linear": "路径线性归因变动",
    "i_share_linear": "路径线性份额变动",
    "i_credit_first": "首次归因偏高",
    "i_share_first": "首次份额偏高",
    "i_n_covisit_neighbors": "共现邻域变动（团伙共点）",
    "ui_pop_mismatch": "热度-活跃错配",
    "i_n_users": "触达用户数变动",
    "i_n_exp": "曝光规模变动",
    "i_log1p_n_users": "触达用户规模(log)",
    "i_log1p_n_exp": "曝光规模(log)",
    "i_log1p_n_covisit": "共现规模(log)",
}


def load_flags(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {"enabled": True, "sample_rate": 1.0, "kill_switch": False}
    return json.loads(path.read_text())


def flags_allow(flags: Dict[str, Any], *, source: str) -> Dict[str, Any]:
    """Gray / rollback gate for card injection (no Drill change)."""
    kill = bool(flags.get("kill_switch"))
    enabled = bool(flags.get("enabled", True)) and not kill
    rate = float(flags.get("sample_rate", 1.0))
    rate = max(0.0, min(1.0, rate))
    # deterministic sample by source path (stable across processes)
    h = int(hashlib.md5(source.encode("utf-8")).hexdigest()[:8], 16) % 10_000
    sampled = (h / 10_000.0) < rate
    allow = enabled and sampled
    return {
        "allow": allow,
        "enabled": enabled,
        "kill_switch": kill,
        "sample_rate": rate,
        "sampled": sampled,
        "reason": (
            "kill_switch"
            if kill
            else ("disabled" if not flags.get("enabled", True) else ("not_sampled" if not sampled else "ok"))
        ),
    }


def _tips_from_summary(blob: Dict[str, Any]) -> List[str]:
    for key in ("fsds_W1_holdout", "fsds_W1", "fsds_W2_temporal", "fsds_W2"):
        block = blob.get(key) or {}
        tops = block.get("top_features")
        if isinstance(tops, list) and tops:
            return [str(x) for x in tops]
    direction = blob.get("direction") or {}
    tip_signs = direction.get("tip_signs") or {}
    if tip_signs:
        return list(tip_signs.keys())
    return []


def _support_summary(blob: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "localize_k": blob.get("localize_k") or (blob.get("k") or {}).get("order"),
        "localized_items_head": (blob.get("localized_items_head") or [])[:15],
        "n_localized_edges": blob.get("n_localized_edges") or blob.get("n_edges"),
        "gap_days": (blob.get("timeline") or {}).get("gap_days", blob.get("gap_days")),
        "k": blob.get("k"),
    }
    return out


def _bucketize(
    tips: List[str],
    tip_signs: Dict[str, str],
    *,
    buckets_map: Optional[Dict[str, str]] = None,
) -> List[Dict[str, str]]:
    mapping = dict(TIP_BUCKETS)
    if buckets_map:
        mapping.update({str(k): str(v) for k, v in buckets_map.items()})
    rows = []
    for t in tips:
        rows.append(
            {
                "feature": t,
                "sign": tip_signs.get(t, "0"),
                "bucket": mapping.get(t, "其它图特征 tip"),
            }
        )
    return rows


def load_tip_overlay(path: Optional[Path]) -> Dict[str, Any]:
    if path is None or not path.exists():
        return {"industry": "default", "buckets": {}}
    blob = json.loads(path.read_text())
    return {
        "industry": str(blob.get("industry") or path.stem),
        "buckets": dict(blob.get("buckets") or {}),
        "note": blob.get("note"),
    }


def build_card(
    blob: Dict[str, Any],
    *,
    source: str,
    flags: Optional[Dict[str, Any]] = None,
    tip_overlay: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    gate = flags_allow(flags or {"enabled": True, "sample_rate": 1.0}, source=source)
    overlay = tip_overlay or {"industry": "default", "buckets": {}}
    direction = dict(blob.get("direction") or {})
    tips = _tips_from_summary(blob)
    tip_signs = dict(direction.get("tip_signs") or {})
    buckets = _bucketize(
        tips[:12], tip_signs, buckets_map=overlay.get("buckets") or None
    )
    sign_dy = direction.get("sign_Dy", "flat")
    primary_bucket = buckets[0]["bucket"] if buckets else "未选 tip"
    card = {
        "card_type": "review_agent_context",
        "disclaimer": "图谱分布变动线索，非定罪结论；供审核分流与上下文，不自动封禁。",
        "source_summary": source,
        "gray_flags": gate,
        "tip_overlay": {
            "industry": overlay.get("industry", "default"),
            "note": overlay.get("note"),
        },
        "support": _support_summary(blob),
        "direction": {
            "Dy": direction.get("Dy"),
            "sign_Dy": sign_dy,
            "y_ref": direction.get("y_ref"),
            "y_cur": direction.get("y_cur"),
            "tip_signs": tip_signs,
            "report": direction.get("report")
            or f"[Direction] sign_Dy={sign_dy}; tips={tips[:5]}",
        },
        "tips": tips[:12],
        "tip_buckets": buckets,
        "review_hint": {
            "queue_bucket": primary_bucket,
            "outcome_read": {
                "pos": "块上成功率上升：优先刷量/互点/末跳队列",
                "neg": "块上成功率下降：优先劣质灌入/劫持残留队列",
                "flat": "X 漂了但 click 率平：先当供给/分布漂移，慎升强动作",
            }.get(sign_dy, "flat"),
            "suggested_action_level": "L1_watch" if gate["allow"] else "L0_observe",
        },
        "paste_for_agent": "",
    }
    if not gate["allow"]:
        card["disclaimer"] = (
            "【灰度关闭/未命中采样】本卡不建议写入工单。" + card["disclaimer"]
        )
        card["review_hint"]["queue_bucket"] = "GRAY_DISABLED"
    lines = [
        "【审核上下文·图谱变动线索】",
        f"灰度: allow={gate['allow']} reason={gate['reason']}",
        f"词典: industry={overlay.get('industry', 'default')}",
        f"方向: sign_Dy={sign_dy} Dy={direction.get('Dy')}",
        f"支撑: localize_k={card['support'].get('localize_k')} "
        f"edges={card['support'].get('n_localized_edges')}",
        f"建议队列: {card['review_hint']['queue_bucket']}",
        f"读法: {card['review_hint']['outcome_read']}",
        "Tips:",
    ]
    for b in buckets[:8]:
        lines.append(f"  - {b['feature']} ({b['sign']}): {b['bucket']}")
    lines.append(f"声明: {card['disclaimer']}")
    card["paste_for_agent"] = "\n".join(lines)
    # SOAR / 工单自定义字段：审出加速接通现网的最小映射
    card["ticket_custom_fields"] = {
        "graph_shift_sign_dy": sign_dy,
        "graph_shift_dy": direction.get("Dy"),
        "graph_shift_queue_bucket": card["review_hint"]["queue_bucket"],
        "graph_shift_action_level": card["review_hint"]["suggested_action_level"],
        "graph_shift_tip_top3": ",".join(tips[:3]),
        "graph_shift_tip_signs_top3": ",".join(
            f"{t}:{tip_signs.get(t, '0')}" for t in tips[:3]
        ),
        "graph_shift_localize_k": card["support"].get("localize_k"),
        "graph_shift_disclaimer": "clue_not_conviction",
        "graph_shift_gray_allow": gate["allow"],
        "graph_shift_gray_reason": gate["reason"],
        "graph_shift_tip_industry": overlay.get("industry", "default"),
    }
    return card


def card_to_md(card: Dict[str, Any]) -> str:
    d = card["direction"]
    tf = card.get("ticket_custom_fields") or {}
    lines = [
        "# 审核 Agent 上下文卡",
        "",
        f"> {card['disclaimer']}",
        "",
        f"- source: `{card['source_summary']}`",
        f"- sign_Dy: **{d.get('sign_Dy')}** (Dy={d.get('Dy')})",
        f"- 建议队列: **{card['review_hint']['queue_bucket']}**",
        f"- 读法: {card['review_hint']['outcome_read']}",
        f"- 建议动作级: `{card['review_hint']['suggested_action_level']}`",
        "",
        "## Tips",
        "",
        "| feature | sign | bucket |",
        "|---|:---:|---|",
    ]
    for b in card["tip_buckets"]:
        lines.append(f"| `{b['feature']}` | {b['sign']} | {b['bucket']} |")
    lines += [
        "",
        "## 工单自定义字段（可直接 POST）",
        "",
        "```json",
        json.dumps(tf, ensure_ascii=False, indent=2),
        "```",
        "",
        "## 粘贴给审核 Agent",
        "",
        "```",
        card["paste_for_agent"],
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description="Export review-agent card from summary.json")
    ap.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "results" / "tencent_gr_w1w2_mmd_po_fsds" / "summary.json",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_review_agent_card",
    )
    ap.add_argument(
        "--flags",
        type=Path,
        default=DEFAULT_FLAGS,
        help="gray/rollback flags JSON",
    )
    ap.add_argument(
        "--tip-overlay",
        type=Path,
        default=None,
        help="optional tip bucket overlay JSON (e.g. content_farm dictionary)",
    )
    args = ap.parse_args()
    blob = json.loads(args.summary.read_text())
    flags = load_flags(args.flags)
    overlay = load_tip_overlay(args.tip_overlay)
    card = build_card(
        blob, source=str(args.summary), flags=flags, tip_overlay=overlay
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "review_agent_card.json").write_text(
        json.dumps(card, indent=2, ensure_ascii=False) + "\n"
    )
    (args.out_dir / "review_agent_card.md").write_text(card_to_md(card))
    (args.out_dir / "ticket_custom_fields.json").write_text(
        json.dumps(card["ticket_custom_fields"], indent=2, ensure_ascii=False) + "\n"
    )
    print(card["paste_for_agent"])
    print("wrote", args.out_dir)


if __name__ == "__main__":
    main()
