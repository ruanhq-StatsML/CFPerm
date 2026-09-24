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
    # --- advertising funnel (Tencent-GR as concrete pack) ---
    "RC_AD_FUNNEL_CONTEXT": {
        "severity": "info",
        "title": "广告漏斗语境：曝→点→转化",
        "family": "ad_funnel",
    },
    "RC_AD_BUY_INTENSITY": {
        "severity": "info",
        "title": "投放/曝光强度基线（买量侧）",
        "family": "ad_funnel",
    },
    "RC_AD_LAST_TOUCH_CANDIDATE": {
        "severity": "investigate",
        "title": "末跳/路径归因候选（落地侧）",
        "family": "ad_funnel",
    },
    "RC_AD_CONVERT_DIP": {
        "severity": "watch",
        "title": "相邻窗转化率下行（灌入/差流倾向）",
        "family": "ad_funnel",
    },
    "RC_AD_CONVERT_SPIKE": {
        "severity": "watch",
        "title": "相邻窗转化率上行（刷量/爆款倾向）",
        "family": "ad_funnel",
    },
    "RC_AD_STRUCTURE_DRIFT": {
        "severity": "info",
        "title": "转化近平 + 结构 tip：投放结构漂移",
        "family": "ad_funnel",
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


def _sign_from_delta(mean_dy: Optional[float], *, thr: float = 0.005) -> str:
    if mean_dy is None:
        return "flat"
    if mean_dy > thr:
        return "pos"
    if mean_dy < -thr:
        return "neg"
    return "flat"


def build_ad_scenario(
    summary: Dict[str, Any], codes: Sequence[Dict[str, Any]]
) -> Dict[str, Any]:
    """Map Tencent-GR board readings → advertising funnel scenario card.

    Concrete pack = TencentGR edges (exp/clk/convert); merchant ≈ advertiser via
    item_feat.122 in upstream docs. Board grain = every N samples along e_last_ts.
    """
    ds = (summary.get("datasets") or {}).get("tencent_gr") or {}
    dsc = (summary.get("datasets") or {}).get("tencent_gr_content") or {}
    if not ds.get("ok"):
        return {
            "ok": False,
            "reason": "tencent_gr panel missing — ad scenario needs convert pack",
        }

    # Prefer N=1000 for finer business grain; fall back
    by = ds.get("by_chunk") or {}
    chunk_key = "1000" if "1000" in by else (next(iter(by), None))
    blk = by.get(chunk_key) or {}
    mean_dy = _finite(blk.get("mean_delta_Y"))
    if mean_dy is None and blk.get("rows"):
        dys = [_finite(r.get("delta_Y")) for r in blk["rows"]]
        dys = [d for d in dys if d is not None]
        if dys:
            mean_dy = float(sum(dys) / len(dys))
    # fall back to cross_pack row
    if mean_dy is None:
        for row in summary.get("cross_pack") or []:
            if row.get("dataset") == "tencent_gr" and str(row.get("chunk")) == str(
                chunk_key
            ):
                mean_dy = _finite(row.get("mean_delta_Y"))
                break
    sign = _sign_from_delta(mean_dy)
    tops_ops = _tops_from_pack(ds, str(chunk_key))
    tops_content = _tops_from_pack(dsc, str(chunk_key)) if dsc.get("ok") else []
    gaps = list(summary.get("ops_content_gap") or [])
    gap1000 = next((g for g in gaps if int(g.get("chunk") or 0) == 1000), None)
    gap = _finite((gap1000 or (gaps[0] if gaps else {})).get("auc_gap_ops_minus_content"))

    code_set = {c["code"] for c in codes}
    has_intensity = "RC_INTENSITY_BASELINE" in code_set or "RC_AD_BUY_INTENSITY" in code_set
    has_credit = "RC_CONTENT_CREDIT_CANDIDATE" in code_set

    if sign == "pos":
        family, family_code = "刷量族", "S1_brush"
        if has_credit:
            sub, sub_code = "末跳/归因嫌疑", "ad_last_touch"
            read = (
                "转化升 + content 末跳 tip：投放侧末跳/刷量队（L1；辨真爆款 vs 刷热）"
            )
        else:
            sub, sub_code = "弱刷量或真爆款", "ad_weak_or_viral"
            read = "转化升但无 credit tip：弱刷量或活动爆款 → 必须人工辨"
    elif sign == "neg":
        family, family_code = "灌入族", "S2_inject"
        sub, sub_code = "转化下行/差流", "ad_convert_dip"
        read = (
            "相邻窗转化掉：先查差流灌入、落地劣化、账户被打压残留——"
            "勿直接判商户/创意「变差」"
        )
    else:
        family, family_code = "漂移族", "S3_drift"
        sub, sub_code = "投放结构漂移", "ad_structure_shift"
        read = (
            "转化近平：结构 tip（强度/归因）在漂，成功率未证实联动 → "
            "查投放策略/定向/创意轮换，慎升强动作"
        )

    funnel = {
        "layers": [
            {"name": "曝光 exp", "feats": ["e_n_exp", "e_log1p_exp", "i_n_exp"]},
            {"name": "点击 clk", "feats": ["e_n_clk", "e_ctr", "u_ctr"]},
            {"name": "转化 convert", "y": "y_convert", "note": "终端成功标签"},
        ],
        "grain": f"every {chunk_key} edges by e_last_ts",
        "entity": "user × item(ad/creative) × merchant≈advertiser",
    }

    paste = [
        "【广告漏斗场景卡 · TencentGR】",
        f"族: {family} (`{family_code}`) / 子类: {sub} (`{sub_code}`)",
        f"读法: {read}",
        f"sign_Dy(board)={sign} · mean_Δȳ={mean_dy} · N={chunk_key}",
        f"ops−content gap={gap}",
        f"ops tops: {tops_ops[:3]}",
        f"content tops: {tops_content[:3]}",
        "动作默认: L1_watch · auto_ban=false · 禁止 HGB 直接上线",
        "强度 = 买量基线；credit/share = 末跳候选 shortlist，非真驱动证书",
        "",
    ]

    return {
        "ok": True,
        "scenario_id": "tencent_gr_ad_funnel",
        "pack": "tencent_gr",
        "chunk": int(chunk_key) if chunk_key else None,
        "funnel": funnel,
        "sign_Dy_board": sign,
        "mean_delta_Y": mean_dy,
        "family": family,
        "family_code": family_code,
        "sub": sub,
        "sub_code": sub_code,
        "read": read,
        "ops_content_gap": gap,
        "has_intensity_baseline": has_intensity,
        "has_last_touch_candidate": has_credit,
        "tops_ops": tops_ops[:5],
        "tops_content": tops_content[:5],
        "action_level_default": "L1_watch",
        "auto_ban": False,
        "allows_ship_model": False,
        "allows_driver_claim": False,
        "paste_for_agent": "\n".join(paste),
        "aligned_with": "direction_report.scenario_from_direction (S1/S2/S3)",
    }


def _append_ad_codes(
    codes: List[Dict[str, Any]], summary: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """Extra AD_* codes when Tencent panels present."""
    ds = (summary.get("datasets") or {}).get("tencent_gr") or {}
    if not ds.get("ok"):
        return codes
    by = ds.get("by_chunk") or {}
    chunk_key = "1000" if "1000" in by else (next(iter(by), None))
    if not chunk_key:
        return codes
    blk = by.get(chunk_key) or {}
    mean_dy = _finite(blk.get("mean_delta_Y"))
    if mean_dy is None:
        for row in summary.get("cross_pack") or []:
            if row.get("dataset") == "tencent_gr" and str(row.get("chunk")) == str(
                chunk_key
            ):
                mean_dy = _finite(row.get("mean_delta_Y"))
                break
    sign = _sign_from_delta(mean_dy)
    tops = _tops_from_pack(ds, str(chunk_key))
    dsc = (summary.get("datasets") or {}).get("tencent_gr_content") or {}
    tops_c = _tops_from_pack(dsc, str(chunk_key)) if dsc.get("ok") else []

    out = list(codes)
    out.append(
        _code(
            "RC_AD_FUNNEL_CONTEXT",
            pack="tencent_gr",
            chunk=int(chunk_key),
            evidence={"funnel": "exp→clk→y_convert", "mean_delta_Y": mean_dy},
            human_copy=(
                "广告场景：边=曝光/点击流，Y=y_convert，商户键≈广告主。"
                "看板按每 N 条边切窗做 transfer 探针，服务审出分流而非出价模型。"
            ),
            next_actions=[
                "用工单话术读 S1/S2/S3，不把 AUC 当 CTR 提升证明",
                "强度与末跳分两行写",
            ],
            tip_feats=tops[:2],
        )
    )
    if any(t in VOLUME_FEATS or "exp" in t for t in tops):
        out.append(
            _code(
                "RC_AD_BUY_INTENSITY",
                pack="tencent_gr",
                chunk=int(chunk_key),
                evidence={"tops": tops[:5]},
                human_copy=(
                    "买量/曝光强度在相邻窗稳传：这是投放量基线。"
                    "可解释「量在不在」，不能单独解释「创意/落地好不好」。"
                ),
                next_actions=["ops 行只报曝光强度", "创意问题看 content 面板"],
                tip_feats=[t for t in tops if t in VOLUME_FEATS or "exp" in t][:3],
            )
        )
    cred = [t for t in tops_c if t in CREDIT_FEATS or "credit" in t or "share" in t]
    if cred:
        out.append(
            _code(
                "RC_AD_LAST_TOUCH_CANDIDATE",
                pack="tencent_gr_content",
                chunk=int(chunk_key),
                evidence={"credit_hits": cred},
                human_copy=(
                    "去量后 credit/share 进 top：广告路径末跳/份额 **候选**。"
                    "对齐审出卡末跳桶；仍 L1，不做自动限投。"
                ),
                next_actions=[
                    "映射 tip 桶 i_share_last / i_credit_last",
                    "人审辨末跳操控 vs 正常回收",
                ],
                tip_feats=cred[:3],
            )
        )
    if sign == "neg":
        out.append(
            _code(
                "RC_AD_CONVERT_DIP",
                pack="tencent_gr",
                chunk=int(chunk_key),
                evidence={"sign_Dy_board": sign, "mean_delta_Y": mean_dy},
                human_copy=(
                    f"相邻窗 mean_Δȳ≈{mean_dy:.4f}（neg）：转化下行倾向。"
                    "广告读法=灌入/差流/落地劣化候选，先 L1 盯梢。"
                ),
                next_actions=["对 S2_inject 队列", "核对定向变更与打压残留"],
                tip_feats=tops_c[:3] or tops[:3],
            )
        )
    elif sign == "pos":
        out.append(
            _code(
                "RC_AD_CONVERT_SPIKE",
                pack="tencent_gr",
                chunk=int(chunk_key),
                evidence={"sign_Dy_board": sign, "mean_delta_Y": mean_dy},
                human_copy=(
                    f"相邻窗 mean_Δȳ≈{mean_dy:.4f}（pos）：转化上行倾向。"
                    "广告读法=刷量/末跳或真爆款，必须人工辨。"
                ),
                next_actions=["对 S1_brush 队列", "有 credit tip 则优先末跳队"],
                tip_feats=tops_c[:3] or tops[:3],
            )
        )
    else:
        out.append(
            _code(
                "RC_AD_STRUCTURE_DRIFT",
                pack="tencent_gr",
                chunk=int(chunk_key),
                evidence={"sign_Dy_board": sign, "mean_delta_Y": mean_dy},
                human_copy=(
                    "转化近平：投放结构（量/归因）可能在漂，成功率未联动证实 → S3 漂移族。"
                ),
                next_actions=["查策略/定向/创意轮换", "禁在 flat 上自称刷量"],
                tip_feats=tops[:3],
            )
        )
    return out


def generate_reason_codes(summary: Dict[str, Any]) -> Dict[str, Any]:
    """Rule engine: board summary.json → ordered reason-code list + ad scenario."""
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

    codes = _append_ad_codes(codes, summary)

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

    ad = build_ad_scenario(summary, uniq)

    paste_lines = [
        "【board reason-codes】",
        "来源: sample-chunk adjacent transfer board（非定罪 / 非上线门）",
        "",
    ]
    if ad.get("ok"):
        paste_lines += [
            "—— 广告漏斗场景 ——",
            ad["paste_for_agent"].rstrip(),
            "",
            "—— 全量 reason-codes ——",
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
            "Routing / claim-control only; allows_ship_model=false always. "
            "Includes TencentGR advertising funnel scenario when pack present."
        ),
        "n_codes": len(uniq),
        "codes": uniq,
        "ad_scenario": ad,
        "paste_for_agent": "\n".join(paste_lines) + "\n",
        "catalog_version": "board_rc_v2_ad",
    }


def render_md(payload: Dict[str, Any]) -> str:
    lines = [
        "# Board reason codes",
        "",
        payload["note"],
        "",
        f"catalog={payload['catalog_version']} · n={payload['n_codes']}",
        "",
    ]
    ad = payload.get("ad_scenario") or {}
    if ad.get("ok"):
        lines += [
            "## Advertising funnel scenario (TencentGR)",
            "",
            f"- family: **{ad.get('family')}** (`{ad.get('family_code')}`)",
            f"- sub: {ad.get('sub')} (`{ad.get('sub_code')}`)",
            f"- sign_Dy(board)={ad.get('sign_Dy_board')} · mean_Δȳ={ad.get('mean_delta_Y')}",
            f"- read: {ad.get('read')}",
            f"- ops−content gap: {ad.get('ops_content_gap')}",
            "",
            "```",
            ad.get("paste_for_agent", "").rstrip(),
            "```",
            "",
        ]
    lines += [
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
    ad = payload.get("ad_scenario") or {}
    if ad.get("ok"):
        (args.out / "ad_scenario.json").write_text(
            json.dumps(ad, indent=2, ensure_ascii=False) + "\n"
        )
        (args.out / "ad_scenario_paste.txt").write_text(ad["paste_for_agent"])
    print(
        json.dumps(
            {
                "n_codes": payload["n_codes"],
                "codes": [c["code"] for c in payload["codes"]],
                "ad_scenario": {
                    "ok": ad.get("ok"),
                    "family_code": ad.get("family_code"),
                    "sign_Dy_board": ad.get("sign_Dy_board"),
                    "read": ad.get("read"),
                },
                "out": str(args.out),
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
