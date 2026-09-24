#!/usr/bin/env python3
"""Auto-generate reason codes from sample-chunk adjacent board summary.

Reason codes are *routing / claim-control labels*, not convictions and not a
ship decision for HGB. They turn board readings (AUC, Jaccard, ops−content gap,
tops) into pasteable codes for review cards / tickets.

  PYTHONPATH=. python3 scripts/generate_board_reason_codes.py \\
    --summary results/sample_chunk_adjacent_board/summary.json \\
    --out results/sample_chunk_adjacent_board/reason_codes

Or called from ``run_sample_chunk_adjacent_board`` after smoke.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]

# Align with scripts/tencent_gr/export_review_agent_card.TIP_BUCKETS (subset + volume)
FEATURE_COPY = {
    "u_span_sec": "活跃跨度异常（短刷/长挂）",
    "i_credit_last": "末次归因偏高（末跳嫌疑）",
    "i_share_last": "末次份额偏高（末跳嫌疑）",
    "i_credit_linear": "路径线性归因变动",
    "i_share_linear": "路径线性份额变动",
    "i_credit_first": "首次归因偏高",
    "i_share_first": "首次份额偏高",
    "i_item_credit_rank": "商户归因秩变动",
    "i_n_covisit_neighbors": "共现邻域变动（团伙共点）",
    "ui_pop_mismatch": "热度-活跃错配",
    "e_log1p_exp": "边曝光强度(log)",
    "e_n_exp": "边曝光次数",
    "e_log1p_clk": "边点击强度(log)",
    "i_log1p_n_exp": "商户曝光规模(log)",
    "i_n_exp": "商户曝光规模",
    "u_n_exp": "用户曝光规模",
    "u_log1p_n_exp": "用户曝光规模(log)",
}

VOLUME_FEATS = {
    "e_log1p_exp",
    "e_n_exp",
    "e_log1p_clk",
    "e_n_clk",
    "i_log1p_n_exp",
    "i_n_exp",
    "u_n_exp",
    "u_log1p_n_exp",
    "u_n_events",
    "i_item_pop_rank",
}
CREDIT_FEATS = {
    "i_credit_last",
    "i_share_last",
    "i_credit_linear",
    "i_share_linear",
    "i_credit_first",
    "i_share_first",
    "i_item_credit_rank",
    "i_log1p_credit_linear",
}

# Catalog: stable codes reviewers / agents can filter on
CATALOG: Dict[str, Dict[str, str]] = {
    "RC_BOARD_NOT_SHIP": {
        "severity": "block_claim",
        "title": "看板禁止直接上线 HGB",
        "family": "ship_gate",
    },
    "RC_SHORTLIST_NOT_DRIVER": {
        "severity": "block_claim",
        "title": "shortlist ≠ 真驱动",
        "family": "driver_gate",
    },
    "RC_WEAK_TRANSFER": {
        "severity": "info",
        "title": "关联下一段传不过去",
        "family": "transfer",
    },
    "RC_DIFFDB_TOKEN_NO_TRAVEL": {
        "severity": "watch",
        "title": "DiffDB prompt token 对 NSFW 不稳传",
        "family": "diffusiondb",
    },
    "RC_INTENSITY_BASELINE": {
        "severity": "info",
        "title": "强度基线（曝光量在传）",
        "family": "tencent_ops",
    },
    "RC_CONTENT_CREDIT_CANDIDATE": {
        "severity": "investigate",
        "title": "content 面板 credit/share 候选",
        "family": "tencent_content",
    },
    "RC_OPS_CONTENT_GAP": {
        "severity": "watch",
        "title": "ops−content 缺口：强度解释大半 transfer",
        "family": "tencent_gap",
    },
    "RC_SHIFTING_DRIVERS": {
        "severity": "watch",
        "title": "高 transfer 但驱动名单在换",
        "family": "stability",
    },
    "RC_PERSISTENT_ASSOCIATION": {
        "severity": "info",
        "title": "关联稳传且 top 稳定",
        "family": "stability",
    },
    "RC_PROBE_DISAGREE": {
        "severity": "watch",
        "title": "HGB 与 LogReg transfer 分歧",
        "family": "probe",
    },
    "RC_RANK_MEAN_SPLIT": {
        "severity": "info",
        "title": "FSDS 排名与 cmean 漂移不一致",
        "family": "selection",
    },
    "RC_PLANTED_SANITY_OK": {
        "severity": "info",
        "title": "种植漂移包探针正常",
        "family": "sanity",
    },
}


def _finite(x: Any) -> Optional[float]:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if v != v:  # NaN
        return None
    return v


def _feat_copy(name: str) -> str:
    if name in FEATURE_COPY:
        return FEATURE_COPY[name]
    # prefix soft match
    for k, v in FEATURE_COPY.items():
        if name.startswith(k) or k in name:
            return v
    if re.match(r"^[a-z0-9_]+$", name) and len(name) < 40:
        return f"特征 `{name}`"
    return f"token `{name}`"


def _tops_from_pack(ds_blk: Dict[str, Any], chunk: str) -> List[str]:
    rows = (ds_blk.get("by_chunk") or {}).get(str(chunk), {}).get("rows") or []
    if not rows:
        return []
    return list(rows[0].get("top_fsds") or [])[:5]


def _code(
    code: str,
    *,
    evidence: Dict[str, Any],
    human_copy: str,
    next_actions: Sequence[str],
    pack: Optional[str] = None,
    chunk: Optional[int] = None,
    tip_feats: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    meta = CATALOG[code]
    tips = [
        {"feature": f, "copy": _feat_copy(f)}
        for f in (tip_feats or [])
        if f
    ]
    return {
        "code": code,
        "family": meta["family"],
        "severity": meta["severity"],
        "title": meta["title"],
        "pack": pack,
        "chunk": chunk,
        "evidence": evidence,
        "human_copy": human_copy,
        "next_actions": list(next_actions),
        "tip_feats": tips,
        "allows_ship_model": False,
        "allows_driver_claim": False,
    }


def generate_reason_codes(summary: Dict[str, Any]) -> Dict[str, Any]:
    """Rule engine: board summary.json → ordered reason-code list."""
    codes: List[Dict[str, Any]] = []
    cross = list(summary.get("cross_pack") or [])
    datasets = dict(summary.get("datasets") or {})
    gaps = list(summary.get("ops_content_gap") or [])
    ship = dict(summary.get("ship_gate") or {})

    # --- always-on claim controls ---
    codes.append(
        _code(
            "RC_BOARD_NOT_SHIP",
            evidence={
                "promote_HGB_to_production": ship.get(
                    "promote_HGB_to_production", False
                ),
                "reason": ship.get("reason"),
            },
            human_copy=(
                "本卡来自相邻切窗 transfer 探针，不是线上打分器。"
                "禁止因 AUC 高直接上线 HGB；缺 label delay / 校准 / Acc·时延 / 灰发回滚。"
            ),
            next_actions=[
                "保持 promote_HGB_to_production=false",
                "若要上线另走 holdout Acc + calibration + gray rollback",
            ],
        )
    )
    codes.append(
        _code(
            "RC_SHORTLIST_NOT_DRIVER",
            evidence={"note": "transfer shortlist only"},
            human_copy=(
                "看板只做特征族 shortlist，不做真驱动认定。"
                "真驱动需 ablation / PO-risk / 干预实验。"
            ),
            next_actions=[
                "shortlist 进 localize/PO 队列",
                "禁止把 top_fsds 写成定罪话术",
            ],
        )
    )

    # --- per pack@N ---
    for row in cross:
        ds = str(row.get("dataset") or "")
        chunk = int(row.get("chunk") or 0)
        auc = _finite(row.get("mean_auc"))
        lr = _finite(row.get("mean_logreg_auc"))
        jac = _finite(row.get("mean_jaccard"))
        fc = _finite(row.get("mean_fsds_cmean_jaccard"))
        dy = _finite(row.get("mean_delta_Y"))
        reading = str(row.get("reading") or "")
        tops = _tops_from_pack(datasets.get(ds) or {}, str(chunk))
        ev_base = {
            "mean_auc": auc,
            "mean_logreg_auc": lr,
            "mean_jaccard": jac,
            "mean_fsds_cmean_jaccard": fc,
            "mean_delta_Y": dy,
            "reading": reading,
            "top_fsds": tops,
        }

        if auc is not None and auc < 0.65:
            codes.append(
                _code(
                    "RC_WEAK_TRANSFER",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        f"`{ds}` @N={chunk}: 下一段 AUC≈{auc:.2f}，关联基本不传。"
                        "不要把本窗 top 特征当成稳定过滤器。"
                    ),
                    next_actions=[
                        "降权该特征族的 tip 叙事",
                        "若业务仍关心 Y，换 X（如 DiffDB 加超参/CLIP）再探针",
                    ],
                    tip_feats=tops[:3],
                )
            )

        if (
            auc is not None
            and auc >= 0.85
            and jac is not None
            and jac >= 0.5
        ):
            codes.append(
                _code(
                    "RC_PERSISTENT_ASSOCIATION",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        f"`{ds}` @N={chunk}: 高 transfer 且 top Jaccard≈{jac:.2f}，"
                        "关联稳传——先当 persistent correlate，勿直接当因果驱动。"
                    ),
                    next_actions=["记入 shortlist", "上 PO/ablation 前保持 L1"],
                    tip_feats=tops[:3],
                )
            )

        if (
            auc is not None
            and auc >= 0.85
            and jac is not None
            and jac < 0.35
        ):
            codes.append(
                _code(
                    "RC_SHIFTING_DRIVERS",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        f"`{ds}` @N={chunk}: AUC 高但驱动名单 Jaccard≈{jac:.2f}，"
                        "预测还在、名单在换——按 regime/composition 读，勿锁死单一 tip。"
                    ),
                    next_actions=[
                        "对比 N=1000 vs 2000 粒度",
                        "人审看当窗 top，不复用上窗话术",
                    ],
                    tip_feats=tops[:3],
                )
            )

        if auc is not None and lr is not None and abs(auc - lr) >= 0.15:
            codes.append(
                _code(
                    "RC_PROBE_DISAGREE",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        f"`{ds}` @N={chunk}: HGB={auc:.2f} vs LogReg={lr:.2f} 分歧≥0.15。"
                        "别把非线性探针优势误写成业务必然；对照 content/线性面板。"
                    ),
                    next_actions=[
                        "并列报告双探针",
                        "优先信任两探针都同意的特征族",
                    ],
                    tip_feats=tops[:3],
                )
            )

        if fc is not None and fc < 0.15:
            codes.append(
                _code(
                    "RC_RANK_MEAN_SPLIT",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        f"`{ds}` @N={chunk}: FSDS∩cmean≈{fc:.2f}，"
                        "排名选中的特征和均值漂移特征不是一回事——查法分开写。"
                    ),
                    next_actions=[
                        "卡上同时贴 top_fsds 与 top_cmean",
                        "禁止只用 cmean 讲「预测力」",
                    ],
                    tip_feats=tops[:3],
                )
            )

        # pack-specific
        if ds == "diffusiondb" and auc is not None and auc < 0.70:
            codes.append(
                _code(
                    "RC_DIFFDB_TOKEN_NO_TRAVEL",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        "DiffusionDB：prompt TF-IDF 对 image_nsfw 弱传。"
                        "风格 token 常进 SelectKBest，但不构成稳定 NSFW 过滤。"
                    ),
                    next_actions=[
                        "NSFW 故事不要押在单一风格 token",
                        "可选：concat CFG/step 或换 embedding 后再跑探针",
                    ],
                    tip_feats=tops[:3],
                )
            )

        if ds == "tencent_gr" and tops:
            vol_hits = [t for t in tops if t in VOLUME_FEATS or "exp" in t]
            if vol_hits and auc is not None and auc >= 0.9:
                codes.append(
                    _code(
                        "RC_INTENSITY_BASELINE",
                        pack=ds,
                        chunk=chunk,
                        evidence={**ev_base, "volume_hits": vol_hits},
                        human_copy=(
                            "Tencent ops：transfer 几乎由曝光/点击强度贡献。"
                            "这是强度基线，不是 content tip；勿写成商户内容问题。"
                        ),
                        next_actions=[
                            "对照 tencent_gr_content 面板看缺口",
                            "强度特征只做 ops 分流，不进定罪桶",
                        ],
                        tip_feats=vol_hits[:3],
                    )
                )

        if ds == "tencent_gr_content" and tops:
            cred = [t for t in tops if t in CREDIT_FEATS or "credit" in t or "share" in t]
            if cred:
                codes.append(
                    _code(
                        "RC_CONTENT_CREDIT_CANDIDATE",
                        pack=ds,
                        chunk=chunk,
                        evidence={**ev_base, "credit_hits": cred},
                        human_copy=(
                            "Tencent content（已去量）：credit/share 进入 top——"
                            "作为末跳/路径归因 **候选** 进审出卡，不做真驱动/L2。"
                        ),
                        next_actions=[
                            "映射 TIP_BUCKETS 末跳/份额话术",
                            "进 localize / PO shortlist",
                            "人审 useful 票后再谈升桶",
                        ],
                        tip_feats=cred[:3],
                    )
                )

        if ds == "waymo_proxy" and auc is not None and auc >= 0.9:
            codes.append(
                _code(
                    "RC_PLANTED_SANITY_OK",
                    pack=ds,
                    chunk=chunk,
                    evidence=ev_base,
                    human_copy=(
                        "Waymo proxy 种植漂移上探针高 transfer——"
                        "看板 sanity 通过；若此处失败，先修探针再读业务包。"
                    ),
                    next_actions=["保留为回归对照包"],
                    tip_feats=tops[:3],
                )
            )

    # --- ops − content gap ---
    for g in gaps:
        gap = _finite(g.get("auc_gap_ops_minus_content"))
        if gap is None:
            continue
        chunk = int(g.get("chunk") or 0)
        if gap >= 0.1:
            codes.append(
                _code(
                    "RC_OPS_CONTENT_GAP",
                    pack="tencent_gr",
                    chunk=chunk,
                    evidence=dict(g),
                    human_copy=(
                        f"ops−content gap≈{gap:.2f} @N={chunk}："
                        "大半 transfer 是强度；content 残差才是归因候选空间。"
                    ),
                    next_actions=[
                        "汇报时先报 gap 再报 content tops",
                        "强度基线与 content 候选分两行写进工单",
                    ],
                    tip_feats=list(g.get("content_top0") or [])[:3],
                )
            )

    # de-dupe by (code, pack, chunk) keeping first
    seen = set()
    uniq: List[Dict[str, Any]] = []
    for c in codes:
        key = (c["code"], c.get("pack"), c.get("chunk"))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(c)

    # severity order for ticket: block_claim → investigate → watch → info
    rank = {"block_claim": 0, "investigate": 1, "watch": 2, "info": 3}
    uniq.sort(key=lambda c: (rank.get(c["severity"], 9), c["code"], str(c.get("pack"))))

    paste_lines = [
        "【board reason-codes】",
        "来源: sample-chunk adjacent transfer board（非定罪 / 非上线门）",
        "",
    ]
    for c in uniq:
        loc = ""
        if c.get("pack"):
            loc = f" `{c['pack']}`"
            if c.get("chunk"):
                loc += f"@N={c['chunk']}"
        paste_lines.append(f"- `{c['code']}` [{c['severity']}]{loc}: {c['human_copy']}")
        for t in c.get("tip_feats") or []:
            paste_lines.append(f"    · {t['feature']}: {t['copy']}")

    return {
        "note": (
            "Reason codes auto-generated from adjacent-board summary. "
            "Routing / claim-control only; allows_ship_model=false always."
        ),
        "n_codes": len(uniq),
        "codes": uniq,
        "paste_for_agent": "\n".join(paste_lines) + "\n",
        "catalog_version": "board_rc_v1",
    }


def render_md(payload: Dict[str, Any]) -> str:
    lines = [
        "# Board reason codes",
        "",
        payload["note"],
        "",
        f"catalog={payload['catalog_version']} · n={payload['n_codes']}",
        "",
        "| code | sev | pack@N | title |",
        "|---|---|---|---|",
    ]
    for c in payload["codes"]:
        loc = c.get("pack") or "—"
        if c.get("chunk"):
            loc = f"{loc}@{c['chunk']}"
        lines.append(
            f"| `{c['code']}` | {c['severity']} | {loc} | {c['title']} |"
        )
    lines += ["", "## Paste for agent", "", "```", payload["paste_for_agent"].rstrip(), "```", ""]
    lines += ["", "## Detail", ""]
    for c in payload["codes"]:
        lines += [
            f"### `{c['code']}` — {c['title']}",
            f"- severity: **{c['severity']}** · family: {c['family']}",
            f"- pack: `{c.get('pack')}` · chunk: `{c.get('chunk')}`",
            f"- copy: {c['human_copy']}",
            f"- next: {', '.join(c['next_actions'])}",
            f"- evidence: `{json.dumps(c['evidence'], ensure_ascii=False)[:240]}`",
            "",
        ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "results/sample_chunk_adjacent_board/summary.json",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results/sample_chunk_adjacent_board/reason_codes",
    )
    args = ap.parse_args()
    summary = json.loads(args.summary.read_text())
    payload = generate_reason_codes(summary)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "reason_codes.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    )
    (args.out / "REASON_CODES.md").write_text(render_md(payload))
    (args.out / "paste_for_agent.txt").write_text(payload["paste_for_agent"])
    print(
        json.dumps(
            {
                "n_codes": payload["n_codes"],
                "codes": [c["code"] for c in payload["codes"]],
                "out": str(args.out),
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
