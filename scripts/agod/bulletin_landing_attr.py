#!/usr/bin/env python3
"""Bulletin landing — agent path / data audit / style drift → text+feature attr.

Maps the polished bulletin to one runnable board:

  1) AI agent continuous inference path  (HaluEval regime + path metrics)
  2) Human/synthetic data quality audit  (HH RFPerm-as-Judge + Top-k)
  3) Style / register drift              (style domain AUC, separate ledger)

Degradation is attributed to **text-metric families** and **feature families**
via RF-domain mass + LOGO Δ (FSDS-lite), then mapped to production actions.

Usage::

    PYTHONPATH=. python3 scripts/agod/bulletin_landing_attr.py
    PYTHONPATH=. python3 scripts/agod/bulletin_landing_attr.py --synth
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "agod" / "bulletin_landing"
DOCS = ROOT / "docs" / "biz"

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

from scripts.agod.hf_landing_protos import (  # noqa: E402
    build_halu_stream,
    build_hh_stream,
    ensure_halu,
    ensure_hh,
    jsonable,
    load_jsonl,
    run_halluc_demo,
    run_judge_demo,
    style_vector,
)

STYLE_NAMES = (
    "tok_len",
    "char_len",
    "avg_word",
    "qmark",
    "bang",
    "hedge",
    "formal",
    "first_person",
    "newlines",
    "upper_ratio",
)

# Feature-family layout used for attribution (aligned with packed X).
# HH pack:  [hash(n_text) | style(10)]
# Halu pack:[hash(n_text) | rag(1) | style(10)]  — we append path(4) here.


def _rf_domain_vimp(X: np.ndarray, w: np.ndarray, *, seed: int = 0) -> tuple[float, np.ndarray]:
    clf = RandomForestClassifier(
        n_estimators=40,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )
    n = len(w)
    idx = np.arange(n)
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    cut = max(n * 3 // 4, 1)
    tr, te = idx[:cut], idx[cut:]
    if len(np.unique(w[tr])) < 2 or len(te) == 0 or len(np.unique(w[te])) < 2:
        clf.fit(X, w)
        return float("nan"), clf.feature_importances_.astype(float)
    clf.fit(X[tr], w[tr])
    auc = float(roc_auc_score(w[te], clf.predict_proba(X[te])[:, 1]))
    return auc, clf.feature_importances_.astype(float)


def _po_risk(X: np.ndarray, y: np.ndarray, w: np.ndarray, *, seed: int = 0) -> dict:
    """Lightweight PO-risk + τ̂ VIMP (concept-relevant axis)."""
    n = len(y)
    yf = y.astype(float)
    # Single-fit surrogate (board speed); production uses CV like po_fs_logo.
    m = RandomForestRegressor(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    e = RandomForestClassifier(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 1, n_jobs=1
    )
    m.fit(X, yf)
    e.fit(X, w)
    m_hat = m.predict(X)
    e_hat = np.clip(e.predict_proba(X)[:, 1], 0.05, 0.95)
    po = (yf - m_hat) * (w.astype(float) - e_hat)
    tau = RandomForestRegressor(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 7, n_jobs=1
    )
    tau.fit(X, po)
    tau_hat = tau.predict(X)
    return {
        "risk": float(np.mean(tau_hat**2)),
        "vimp": tau.feature_importances_.astype(float),
    }


def _family_mass(vimp: np.ndarray, families: list[str], groups: dict[str, list[int]]) -> dict:
    out = {g: 0.0 for g in families}
    for g, ix in groups.items():
        out[g] = float(np.sum(vimp[ix])) if ix else 0.0
    s = sum(out.values()) + 1e-12
    return {k: out[k] / s for k in families}


def _logo_share(
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    families: list[str],
    groups: dict[str, list[int]],
    *,
    seed: int,
) -> dict:
    full = _po_risk(X, y, w, seed=seed)
    logo = {}
    for i, g in enumerate(families):
        drop = set(groups[g])
        keep = [j for j in range(X.shape[1]) if j not in drop]
        if not keep:
            logo[g] = {"delta": 0.0, "R_minus": full["risk"]}
            continue
        rm = _po_risk(X[:, keep], y, w, seed=seed + 11 + i)["risk"]
        logo[g] = {"delta": float(full["risk"] - rm), "R_minus": float(rm)}
    pos = {g: max(logo[g]["delta"], 0.0) for g in families}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in families}
    return {"po_risk": full["risk"], "logo": logo, "logo_share": share, "po_vimp": full["vimp"]}


def _top_named(vimp: np.ndarray, names: list[str], k: int = 8) -> list[dict]:
    order = np.argsort(-vimp)
    out = []
    for r, j in enumerate(order[:k]):
        out.append({"rank": r + 1, "name": names[j], "vimp": float(vimp[j])})
    return out


def _metric_shift(pre: np.ndarray, post: np.ndarray, names: tuple[str, ...]) -> list[dict]:
    rows = []
    for i, name in enumerate(names):
        a, b = float(np.mean(pre[:, i])), float(np.mean(post[:, i]))
        rows.append(
            {
                "metric": name,
                "mean_pre": a,
                "mean_post": b,
                "delta": b - a,
                "abs_delta": abs(b - a),
            }
        )
    rows.sort(key=lambda r: -r["abs_delta"])
    return rows


def _pick_top_family(logo_share: dict, rf_mass: dict) -> str:
    """Prefer LOGO share; fall back to RF-domain mass when LOGO is flat."""
    if max(logo_share.values(), default=0.0) > 1e-6:
        return max(logo_share, key=logo_share.get)
    return max(rf_mass, key=rf_mass.get)


def _action_for(family: str, scenario: str) -> str:
    table = {
        ("agent", "rag"): "刷检索 / 补知识 / 换索引（别先回滚整模）",
        ("agent", "path"): "审工具调用链 / 限步深 / 修路由专家",
        ("agent", "style"): "改 decoding / 素材口径（另账，不进客服主账）",
        ("agent", "text_hash"): "抽审 Top-k 会话；必要时 canary / 回滚生成侧",
        ("audit", "style"): "文风工单：素材/register，不挡偏好合并",
        ("audit", "text_hash"): "偏好映射门禁：拒合并 / 限量 / Top-k 人工复核",
        ("style", "style"): "画风/文风刷新：改模板、decoding、创意批次",
        ("style", "text_hash"): "语义内容轴另开；勿与画风混账",
    }
    return table.get((scenario, family), "子集定位后只动该族")


# ---------------------------------------------------------------------------
# Scenario builders
# ---------------------------------------------------------------------------


def _path_metrics(n: int, batch: np.ndarray, cut: int, rng: np.random.Generator) -> np.ndarray:
    """Synthetic agent-path metrics that hop with the inference regime.

    Columns: step_depth, tool_call_rate, retry_rate, route_entropy.
    After cut: deeper chains + more tools + more retries (path decay signature).
    """
    base = np.column_stack(
        [
            rng.uniform(0.1, 0.35, n),
            rng.uniform(0.05, 0.2, n),
            rng.uniform(0.0, 0.1, n),
            rng.uniform(0.2, 0.45, n),
        ]
    )
    after = batch >= cut
    base[after, 0] += rng.uniform(0.25, 0.45, after.sum())
    base[after, 1] += rng.uniform(0.2, 0.4, after.sum())
    base[after, 2] += rng.uniform(0.15, 0.35, after.sum())
    base[after, 3] += rng.uniform(0.15, 0.3, after.sum())
    return np.clip(base, 0.0, 1.5)


def run_agent_path_landing(halu_rows: list[dict], *, gate: float, seed: int) -> dict:
    pack = build_halu_stream(halu_rows, n_per=80, cut_batch=4, seed=seed, n_features=96)
    rng = np.random.default_rng(seed + 3)
    path = _path_metrics(len(pack["y"]), pack["batch"], pack["cut_batch"], rng)
    pack["X"] = np.hstack([pack["X"], path])
    pack["path"] = path

    demo = run_halluc_demo(pack, gate=gate, seed=seed)
    X, y, batch = pack["X"], pack["y"], pack["batch"]
    cut = int(pack["cut_batch"])
    w = (batch >= cut).astype(int)
    n_text = 96
    # Layout: hash | rag(1) | style(10) | path(4)
    names = (
        [f"hash_{i}" for i in range(n_text)]
        + ["rag_hit"]
        + list(STYLE_NAMES)
        + ["step_depth", "tool_call_rate", "retry_rate", "route_entropy"]
    )
    families = ["text_hash", "rag", "style", "path"]
    groups = {
        "text_hash": list(range(0, n_text)),
        "rag": [n_text],
        "style": list(range(n_text + 1, n_text + 11)),
        "path": list(range(n_text + 11, n_text + 15)),
    }
    auc, v_rf = _rf_domain_vimp(X, w, seed=seed)
    logo = _logo_share(X, y, w, families, groups, seed=seed)
    rf_mass = _family_mass(v_rf, families, groups)
    style_cols = X[:, groups["style"]]
    path_cols = X[:, groups["path"]]
    pre, post = batch < cut, batch >= cut
    text_shifts = _metric_shift(style_cols[pre], style_cols[post], STYLE_NAMES)
    path_shifts = _metric_shift(
        path_cols[pre],
        path_cols[post],
        ("step_depth", "tool_call_rate", "retry_rate", "route_entropy"),
    )
    top_family = _pick_top_family(logo["logo_share"], rf_mass)
    return {
        "bulletin_domain": "AI智能体连续推理路径",
        "method": "OnlineRFPerm + FSDS-lite (RF-domain mass / LOGO)",
        "detection": {
            "first_fire_t": demo["first_fire_t"],
            "fire_rate": demo["fire_rate"],
            "rate_before": demo["rate_before"],
            "rate_after": demo["rate_after"],
            "hop_at_cut": demo.get("hop_at_cut"),
        },
        "attribution": {
            "rf_domain_auc": auc,
            "rf_mass": rf_mass,
            "po_risk": logo["po_risk"],
            "logo_share": logo["logo_share"],
            "logo": logo["logo"],
            "top_family": top_family,
            "top_features_rf": _top_named(v_rf, names, 8),
            "text_metric_shifts": text_shifts[:5],
            "path_metric_shifts": path_shifts,
        },
        "action": {
            "primary": _action_for(top_family, "agent"),
            "ledger": "客服/助手主账（工单·退款·承接）",
            "not": "不是新的 LLM 训练法；不是全空间唯一因果",
        },
    }


def run_audit_landing(hh_rows: list[dict], *, gate: float, seed: int) -> dict:
    pack = build_hh_stream(hh_rows, n_per=60, cut_batch=4, seed=seed, n_features=96)
    demo = run_judge_demo(pack, gate=gate, seed=seed)
    X, y, batch = pack["X"], pack["y"], pack["batch"]
    cut = int(pack["cut_batch"])
    w = (batch >= cut).astype(int)
    n_text = 96
    names = [f"hash_{i}" for i in range(n_text)] + list(STYLE_NAMES)
    families = ["text_hash", "style"]
    groups = {
        "text_hash": list(range(0, n_text)),
        "style": list(range(n_text, n_text + 10)),
    }
    auc, v_rf = _rf_domain_vimp(X, w, seed=seed)
    logo = _logo_share(X, y, w, families, groups, seed=seed)
    rf_mass = _family_mass(v_rf, families, groups)
    # Preference map hop lives in Y|X — style AUC is the separate portrait axis.
    ranking = (demo.get("hop_at_cut") or {}).get("ranking") or {}
    top_family = _pick_top_family(logo["logo_share"], rf_mass)
    # Prefer preference-axis action when judge ratio clearly fired.
    if demo["judge_err_ratio"] > 1.5:
        primary_action = _action_for("text_hash", "audit")
        primary_family = "text_hash (preference map)"
    else:
        primary_action = _action_for(top_family, "audit")
        primary_family = top_family
    return {
        "bulletin_domain": "人工/合成数据质量检测与审核",
        "method": "RFPerm-as-Judge + PO Top-k + style/preference split",
        "detection": {
            "first_fire_t": demo["first_fire_t"],
            "fire_rate": demo["fire_rate"],
            "judge_err_pre": demo["judge_err_pre"],
            "judge_err_post": demo["judge_err_post"],
            "judge_err_ratio": demo["judge_err_ratio"],
            "style_domain_auc": demo["style_domain_auc"],
            "hop_ranking": ranking,
        },
        "attribution": {
            "rf_domain_auc": auc,
            "rf_mass": rf_mass,
            "po_risk": logo["po_risk"],
            "logo_share": logo["logo_share"],
            "logo": logo["logo"],
            "top_family": top_family,
            "preference_vs_style": {
                "preference_signal": "judge_err_ratio / OnlineRFPerm fire",
                "style_signal": "style_domain_auc (P(X) portrait)",
                "split": "两本账：偏好门禁 ≠ 文风工单",
            },
            "top_features_rf": _top_named(v_rf, names, 8),
        },
        "action": {
            "primary": primary_action,
            "primary_family": primary_family,
            "topk_note": "fire 后按 po_risk0 Top-k 人工复核；安静回均匀抽检",
            "ledger": "对齐/合并门禁防损账",
        },
    }


def run_style_landing(hh_rows: list[dict], *, seed: int) -> dict:
    """Style/register drift board — covariant portrait only, separate ledger."""
    styles = []
    for r in hh_rows:
        for key in ("chosen", "rejected"):
            from scripts.agod.hf_landing_protos import assistant_reply

            styles.append(style_vector(assistant_reply(r[key])))
    styles = np.vstack(styles)
    # Explicit register contrast: low vs high formality terciles.
    formality = styles[:, 6] + 0.5 * styles[:, 1]
    q_lo, q_hi = np.quantile(formality, [0.33, 0.67])
    lo, hi = formality <= q_lo, formality >= q_hi
    X = styles
    w = np.zeros(len(styles), dtype=int)
    w[hi] = 1
    # Keep only lo/hi for clean domain contrast.
    mask = lo | hi
    X_m, w_m = X[mask], w[mask]
    auc, v_rf = _rf_domain_vimp(X_m, w_m, seed=seed)
    # LOGO on style metric families: length / punct / register
    names = list(STYLE_NAMES)
    families = ["length", "punct", "register"]
    groups = {
        "length": [0, 1, 2],
        "punct": [3, 4, 8],
        "register": [5, 6, 7, 9],
    }
    # Fake Y = formality band for concept-lite LOGO (style-only board).
    y_m = w_m.copy()
    logo = _logo_share(X_m, y_m, w_m, families, groups, seed=seed)
    rf_mass = _family_mass(v_rf, families, groups)
    shifts = _metric_shift(X[lo], X[hi], STYLE_NAMES)
    top_family = _pick_top_family(logo["logo_share"], rf_mass)
    top_metric = shifts[0]["metric"] if shifts else "formal"
    return {
        "bulletin_domain": "画风/文风漂移监测",
        "method": "style domain AUC + metric LOGO (covariate board)",
        "detection": {
            "style_domain_auc": auc,
            "n_lo": int(lo.sum()),
            "n_hi": int(hi.sum()),
        },
        "attribution": {
            "rf_mass": rf_mass,
            "logo_share": logo["logo_share"],
            "logo": logo["logo"],
            "top_family": top_family,
            "top_metric": top_metric,
            "text_metric_shifts": shifts[:6],
            "top_features_rf": _top_named(v_rf, names, 8),
        },
        "action": {
            "primary": f"改 {top_metric}/{top_family} 对应的模板或 decoding",
            "ledger": "CTR/品牌另账 — 永不并进客服主账",
            "not": "不要把文风信号拿去误伤偏好头或客服回滚",
        },
    }


def _synth_hh(n_pairs: int = 280) -> list[dict]:
    rows = []
    for i in range(n_pairs):
        chosen = (
            f"\n\nHuman: Q{i} about topic {i % 7}?\n\n"
            f"Assistant: Therefore I recommend option {i} with details regarding the plan."
        )
        rejected = (
            f"\n\nHuman: Q{i} about topic {i % 7}?\n\n"
            f"Assistant: maybe idk lol {i}!!! not sure"
        )
        rows.append({"chosen": chosen, "rejected": rejected})
    return rows


def _synth_halu(n: int = 480) -> list[dict]:
    rows = []
    for i in range(n):
        know = f"entity {i} was founded in {1900 + (i % 50)} in city {i % 11}"
        q = f"When was entity {i} founded?"
        if i % 2 == 0:
            ans = f"entity {i} was founded in {1900 + (i % 50)}"
            hall = "no"
        else:
            ans = f"entity {i} opened a mall in 2099"
            hall = "yes"
        rows.append(
            {"knowledge": know, "question": q, "answer": ans, "hallucination": hall}
        )
    return rows


def write_report(summary: dict, path: Path) -> None:
    a = summary["agent_path"]
    u = summary["data_audit"]
    s = summary["style_drift"]
    lines = [
        "# Bulletin 落地跑板 — 检测 → 文本/特征归因 → 动作",
        "",
        "## 对外主句",
        "",
        "> " + summary["bulletin"],
        "",
        "## 1) AI 智能体连续推理路径",
        "",
        f"- 火情 first_fire_t=`{a['detection']['first_fire_t']}` "
        f"rate `{a['detection']['rate_before']:.3f}`→`{a['detection']['rate_after']:.3f}`",
        f"- RF-domain AUC=`{a['attribution']['rf_domain_auc']}`",
        f"- 特征族份额 (LOGO): `{a['attribution']['logo_share']}`",
        f"- 主归因族: **{a['attribution']['top_family']}** → {a['action']['primary']}",
        f"- 路径指标漂移 Top: `{a['attribution']['path_metric_shifts'][:3]}`",
        f"- 文本指标漂移 Top: `{a['attribution']['text_metric_shifts'][:3]}`",
        "",
        "## 2) 人工/合成数据质检与审核",
        "",
        f"- judge err `{u['detection']['judge_err_pre']:.3f}`→`{u['detection']['judge_err_post']:.3f}` "
        f"(ratio `{u['detection']['judge_err_ratio']:.2f}`)",
        f"- style_domain_auc=`{u['detection']['style_domain_auc']:.3f}`（画像轴另账）",
        f"- LOGO 份额: `{u['attribution']['logo_share']}`",
        f"- 动作: {u['action']['primary']}",
        "",
        "## 3) 画风/文风漂移",
        "",
        f"- style_domain_auc=`{s['detection']['style_domain_auc']}`",
        f"- 指标族份额: `{s['attribution']['logo_share']}`",
        f"- Top 文本指标: `{s['attribution']['top_metric']}`",
        f"- 动作: {s['action']['primary']}（{s['action']['ledger']}）",
        "",
        "## 口径",
        "",
        "检测打开对照窗 → 子集/特征族定位 → 只动该子集的刷新/审核/门禁。",
        "不是新 LLM 训练法，不是全空间唯一因果根因。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synth", action="store_true", help="offline synthetic rows")
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    bulletin = (
        "将分布漂移检测与子集定位应用于 AI 智能体连续推理路径、"
        "人工/合成数据质量检测与审核，以及画风/文风漂移监测；"
        "在推理质量衰减时，能把问题归因到文本指标与特征层面，"
        "并直接接到刷新、审核与经营动作，具备可落地的生产价值。"
    )

    if args.synth:
        hh_rows, halu_rows = _synth_hh(320), _synth_halu(560)
        source = "synth"
    else:
        hh_path = ensure_hh()
        hu_path = ensure_halu()
        hh_rows = load_jsonl(hh_path, n=1600)
        halu_rows = load_jsonl(hu_path, n=1600)
        source = "hf_cache"

    agent = run_agent_path_landing(halu_rows, gate=args.gate, seed=args.seed)
    audit = run_audit_landing(hh_rows, gate=args.gate, seed=args.seed)
    style = run_style_landing(hh_rows, seed=args.seed)

    summary = {
        "bulletin": bulletin,
        "source": source,
        "gate": args.gate,
        "seed": args.seed,
        "agent_path": agent,
        "data_audit": audit,
        "style_drift": style,
        "landing_map": {
            "检测": "OnlineRFPerm / RFPerm / style domain AUC",
            "归因": "文本指标漂移表 + 特征族 RF-mass / LOGO share",
            "动作": "刷新检索·审路径·Top-k 审核·文风另账",
        },
    }
    (OUT / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2), encoding="utf-8"
    )
    write_report(summary, OUT / "REPORT.md")
    # Mirror a short pointer into the bulletin polish doc section is done by docs edit.

    print(
        f"[agent] fire={agent['detection']['first_fire_t']} "
        f"top_family={agent['attribution']['top_family']} "
        f"logo={agent['attribution']['logo_share']}"
    )
    print(
        f"[audit] ratio={audit['detection']['judge_err_ratio']:.2f} "
        f"style_auc={audit['detection']['style_domain_auc']:.3f} "
        f"action={audit['action']['primary_family']}"
    )
    print(
        f"[style] auc={style['detection']['style_domain_auc']} "
        f"top_metric={style['attribution']['top_metric']} "
        f"logo={style['attribution']['logo_share']}"
    )
    print(f"wrote {OUT / 'REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
