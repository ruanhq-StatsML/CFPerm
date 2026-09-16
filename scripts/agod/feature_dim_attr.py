#!/usr/bin/env python3
"""Unified feature-dimension attribution — style + graph fold into one layer.

Stance
------
Style attribution and graph attribution are **not** two products.
Both first collapse into a shared **feature-dimension** catalog; then the same
L1 localization → L2 FSDS/LOGO drill runs on that catalog.

  style metrics  ──┐
  graph nodes    ──┼──► feature dimensions D = {d1..dk} ──► mass / LOGO ──► action
  path / rag     ──┘

Usage::

    PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py
    PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py --synth
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "agod" / "feature_dim_attr"

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

from scripts.agod.hf_landing_protos import (  # noqa: E402
    assistant_reply,
    build_halu_stream,
    build_hh_stream,
    ensure_halu,
    ensure_hh,
    jsonable,
    load_jsonl,
    style_vector,
)


# ---------------------------------------------------------------------------
# Unified feature-dimension catalog
# ---------------------------------------------------------------------------

STYLE_METRIC_NAMES = (
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

# Style metrics themselves are feature dims under the Style graph node.
STYLE_SUBFAMILIES = {
    "length": ("tok_len", "char_len", "avg_word"),
    "punct": ("qmark", "bang", "newlines"),
    "register": ("hedge", "formal", "first_person", "upper_ratio"),
}

# Graph nodes (business or LLM-isomorphic) = top-level feature dimensions.
# Merchant/Author/Product/Order map to the same slots as LLM towers/paths.
GRAPH_DIMS = (
    "merchant",  # supply / knowledge / retrieval side
    "author",  # creator / style / register side
    "product",  # SKU / content unit / answer payload
    "order",  # conversion / path / tool-chain outcome
)

# LLM landing aliases — same dimension indices, different names on the board.
LLM_DIM_ALIAS = {
    "merchant": "rag_retrieval",
    "author": "style_register",
    "product": "text_payload",
    "order": "agent_path",
}

ACTION_BY_DIM = {
    "merchant": "刷检索 / 补知识 / 修供给侧特征",
    "author": "改素材·decoding / 作者侧限流（文风另账）",
    "product": "审内容单元 / Top-k 抽检 payload",
    "order": "审工具链·路径时延 / 转化粒策略",
    "rag_retrieval": "刷检索 / 换索引（别先回滚整模）",
    "style_register": "改模板·decoding（另账，不进客服主账）",
    "text_payload": "抽审会话文本 / canary 生成侧",
    "agent_path": "限步深 / 修路由专家 / 审 tool 调用",
}


@dataclass
class FeatureDimPack:
    """One matrix X with a dimension → column-index map (the unified layer)."""

    X: np.ndarray
    y: np.ndarray
    w: np.ndarray  # pre/post or domain contrast
    names: list[str]
    dims: dict[str, list[int]]  # dimension → column indices
    meta: dict = field(default_factory=dict)

    @property
    def dim_names(self) -> list[str]:
        return list(self.dims.keys())


def rf_domain_vimp(X: np.ndarray, w: np.ndarray, *, seed: int = 0) -> tuple[float, np.ndarray]:
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=6, min_samples_leaf=2, random_state=seed, n_jobs=1
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


def po_risk(X: np.ndarray, y: np.ndarray, w: np.ndarray, *, seed: int = 0) -> dict:
    yf = y.astype(float)
    m = RandomForestRegressor(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    e = RandomForestClassifier(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 1, n_jobs=1
    )
    m.fit(X, yf)
    e.fit(X, w)
    po = (yf - m.predict(X)) * (w.astype(float) - np.clip(e.predict_proba(X)[:, 1], 0.05, 0.95))
    tau = RandomForestRegressor(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 7, n_jobs=1
    )
    tau.fit(X, po)
    return {"risk": float(np.mean(tau.predict(X) ** 2)), "vimp": tau.feature_importances_.astype(float)}


def dim_mass(vimp: np.ndarray, dims: dict[str, list[int]]) -> dict[str, float]:
    out = {d: float(np.sum(vimp[ix])) if ix else 0.0 for d, ix in dims.items()}
    s = sum(out.values()) + 1e-12
    return {k: v / s for k, v in out.items()}


def logo_by_dim(
    X: np.ndarray, y: np.ndarray, w: np.ndarray, dims: dict[str, list[int]], *, seed: int
) -> dict:
    full = po_risk(X, y, w, seed=seed)
    logo = {}
    for i, (d, ix) in enumerate(dims.items()):
        drop = set(ix)
        keep = [j for j in range(X.shape[1]) if j not in drop]
        if not keep:
            logo[d] = {"delta": 0.0, "R_minus": full["risk"]}
            continue
        rm = po_risk(X[:, keep], y, w, seed=seed + 11 + i)["risk"]
        logo[d] = {"delta": float(full["risk"] - rm), "R_minus": float(rm)}
    pos = {d: max(v["delta"], 0.0) for d, v in logo.items()}
    s = sum(pos.values()) + 1e-12
    share = {d: pos[d] / s for d in dims}
    return {"po_risk": full["risk"], "logo": logo, "logo_share": share, "po_vimp": full["vimp"]}


def pick_top_dims(logo_share: dict, rf_mass: dict, k: int = 2) -> list[str]:
    """L1: pick top-k feature dimensions (LOGO first, RF-mass fallback)."""
    use = logo_share if max(logo_share.values(), default=0.0) > 1e-6 else rf_mass
    return [d for d, _ in sorted(use.items(), key=lambda kv: -kv[1])[:k]]


def metric_shift(pre: np.ndarray, post: np.ndarray, names: tuple[str, ...] | list[str]) -> list[dict]:
    rows = []
    for i, name in enumerate(names):
        a, b = float(np.mean(pre[:, i])), float(np.mean(post[:, i]))
        rows.append({"metric": name, "mean_pre": a, "mean_post": b, "delta": b - a, "abs_delta": abs(b - a)})
    rows.sort(key=lambda r: -r["abs_delta"])
    return rows


def attribute_pack(pack: FeatureDimPack, *, seed: int = 0, top_k: int = 2) -> dict:
    """L1 dim localization + L2 within-dim mass — one pipeline for style & graph."""
    auc, v_rf = rf_domain_vimp(pack.X, pack.w, seed=seed)
    logo = logo_by_dim(pack.X, pack.y, pack.w, pack.dims, seed=seed)
    rf_mass = dim_mass(v_rf, pack.dims)
    top = pick_top_dims(logo["logo_share"], rf_mass, k=top_k)
    # L2: within each top dim, rank named features
    within = {}
    for d in top:
        ix = pack.dims[d]
        local = v_rf[ix]
        local_names = [pack.names[j] for j in ix]
        order = np.argsort(-local)
        within[d] = [
            {"rank": r + 1, "name": local_names[j], "vimp": float(local[j])}
            for r, j in enumerate(order[:5])
        ]
    actions = [ACTION_BY_DIM.get(d, f"只刷新维度 {d}") for d in top]
    return {
        "rf_domain_auc": auc,
        "rf_mass": rf_mass,
        "logo_share": logo["logo_share"],
        "logo": logo["logo"],
        "po_risk": logo["po_risk"],
        "top_dims": top,
        "within_dim": within,
        "actions": actions,
        "meta": pack.meta,
    }


# ---------------------------------------------------------------------------
# Build packs: style folded into graph dims; graph dims = feature dims
# ---------------------------------------------------------------------------


def build_style_as_feature_dims(hh_rows: list[dict]) -> FeatureDimPack:
    """Style attribution → feature dimensions (length / punct / register)."""
    styles = []
    for r in hh_rows:
        for key in ("chosen", "rejected"):
            styles.append(style_vector(assistant_reply(r[key])))
    styles = np.vstack(styles)
    formality = styles[:, 6] + 0.5 * styles[:, 1]
    q_lo, q_hi = np.quantile(formality, [0.33, 0.67])
    lo, hi = formality <= q_lo, formality >= q_hi
    mask = lo | hi
    X = styles[mask]
    w = np.zeros(mask.sum(), dtype=int)
    # hi among masked rows
    w[formality[mask] >= q_hi] = 1
    name_to_i = {n: i for i, n in enumerate(STYLE_METRIC_NAMES)}
    dims = {fam: [name_to_i[n] for n in names] for fam, names in STYLE_SUBFAMILIES.items()}
    return FeatureDimPack(
        X=X,
        y=w.copy(),
        w=w,
        names=list(STYLE_METRIC_NAMES),
        dims=dims,
        meta={"board": "style_as_feature_dims", "ledger": "CTR/品牌另账"},
    )


def build_graph_feature_dims(
    halu_rows: list[dict], *, seed: int = 0, alias: str = "llm"
) -> FeatureDimPack:
    """Graph nodes → feature dimensions; style is the Author dim, not a side channel.

    Layout per row:
      merchant/rag : rag_hit (1)
      author/style : style_vector (10)
      product/text : hash text proxy — here we use answer length stats from hash
                     board: take first 8 hashing dims as payload proxy
      order/path   : step_depth, tool_call_rate, retry_rate, route_entropy (4)
    """
    pack = build_halu_stream(halu_rows, n_per=80, cut_batch=4, seed=seed, n_features=32)
    rng = np.random.default_rng(seed + 3)
    n = len(pack["y"])
    batch = pack["batch"]
    cut = int(pack["cut_batch"])
    # From packed X: hash(32) | rag(1) | style(10)
    X_base = pack["X"]
    rag = X_base[:, 32:33]
    style = X_base[:, 33:43]
    text = X_base[:, :8]  # payload proxy
    path = np.column_stack(
        [
            rng.uniform(0.1, 0.35, n),
            rng.uniform(0.05, 0.2, n),
            rng.uniform(0.0, 0.1, n),
            rng.uniform(0.2, 0.45, n),
        ]
    )
    after = batch >= cut
    path[after] += rng.uniform(0.2, 0.4, size=(after.sum(), 4))
    path = np.clip(path, 0.0, 1.5)

    X = np.hstack([rag, style, text, path])
    names = (
        ["rag_hit"]
        + list(STYLE_METRIC_NAMES)
        + [f"payload_h{i}" for i in range(8)]
        + ["step_depth", "tool_call_rate", "retry_rate", "route_entropy"]
    )
    raw_dims = {
        "merchant": [0],
        "author": list(range(1, 11)),
        "product": list(range(11, 19)),
        "order": list(range(19, 23)),
    }
    if alias == "llm":
        dims = {LLM_DIM_ALIAS[k]: v for k, v in raw_dims.items()}
    else:
        dims = raw_dims
    w = (batch >= cut).astype(int)
    return FeatureDimPack(
        X=X,
        y=pack["y"],
        w=w,
        names=names,
        dims=dims,
        meta={
            "board": "graph_as_feature_dims",
            "alias": alias,
            "graph_dims": list(GRAPH_DIMS),
            "llm_alias": LLM_DIM_ALIAS,
            "note": "style folded into author/style_register dimension — not a parallel product",
        },
    )


def build_audit_feature_dims(hh_rows: list[dict], *, seed: int = 0) -> FeatureDimPack:
    """Preference audit: text_hash vs style as two feature dimensions."""
    pack = build_hh_stream(hh_rows, n_per=60, cut_batch=4, seed=seed, n_features=64)
    X, y, batch = pack["X"], pack["y"], pack["batch"]
    cut = int(pack["cut_batch"])
    names = [f"hash_{i}" for i in range(64)] + list(STYLE_METRIC_NAMES)
    dims = {
        "preference_text": list(range(0, 64)),
        "style_register": list(range(64, 74)),
    }
    return FeatureDimPack(
        X=X,
        y=y,
        w=(batch >= cut).astype(int),
        names=names,
        dims=dims,
        meta={
            "board": "audit_feature_dims",
            "split": "偏好维 vs 文风维 — 同一特征维度层，两本账",
        },
    )


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
    g = summary["graph_dims"]
    s = summary["style_dims"]
    a = summary["audit_dims"]
    lines = [
        "# 特征维度统一归因 — 风格 + 图谱先并维",
        "",
        "> 风格归因与图谱归因 **先整合为特征维度**，再跑同一套 L1→L2。",
        "",
        "## 并维口径",
        "",
        "| 来源 | 收成的特征维度 |",
        "|------|----------------|",
        "| 画风/文风指标 | length / punct / register（Style 节点下的子维） |",
        "| 商户·作者·商品·订单图 | merchant / author / product / order |",
        "| LLM 同构别名 | rag_retrieval / style_register / text_payload / agent_path |",
        "",
        "## Graph → 特征维度",
        "",
        f"- Top dims: `{g['top_dims']}`",
        f"- LOGO share: `{g['logo_share']}`",
        f"- RF mass: `{g['rf_mass']}`",
        f"- Actions: {g['actions']}",
        f"- Within: `{g['within_dim']}`",
        "",
        "## Style → 特征维度（不是旁路产品）",
        "",
        f"- Top dims: `{s['top_dims']}`",
        f"- LOGO/RF: logo=`{s['logo_share']}` rf=`{s['rf_mass']}`",
        f"- Actions: {s['actions']}",
        "",
        "## Audit → 特征维度（偏好维 vs 文风维）",
        "",
        f"- Top dims: `{a['top_dims']}`",
        f"- share: `{a['logo_share']}` / `{a['rf_mass']}`",
        f"- Actions: {a['actions']}",
        "",
        "## 一句",
        "",
        "先并到特征维度，再 localization；风格与图谱是同一层目录上的不同节点，不是两套归因系统。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synth", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alias", choices=["llm", "graph"], default="llm")
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    if args.synth:
        hh_rows, hu_rows = _synth_hh(300), _synth_halu(560)
        source = "synth"
    else:
        hh_rows = load_jsonl(ensure_hh(), n=1600)
        hu_rows = load_jsonl(ensure_halu(), n=1600)
        source = "hf_cache"

    style_pack = build_style_as_feature_dims(hh_rows)
    graph_pack = build_graph_feature_dims(hu_rows, seed=args.seed, alias=args.alias)
    audit_pack = build_audit_feature_dims(hh_rows, seed=args.seed)

    style_attr = attribute_pack(style_pack, seed=args.seed)
    graph_attr = attribute_pack(graph_pack, seed=args.seed, top_k=2)
    audit_attr = attribute_pack(audit_pack, seed=args.seed)

    # Attach text-metric shifts for style board readability
    formality = style_pack.X[:, 6] + 0.5 * style_pack.X[:, 1]
    med = np.median(formality)
    style_attr["text_metric_shifts"] = metric_shift(
        style_pack.X[formality <= med],
        style_pack.X[formality > med],
        STYLE_METRIC_NAMES,
    )[:6]

    summary = {
        "stance": "style + graph first unify into feature dimensions, then L1→L2",
        "source": source,
        "seed": args.seed,
        "catalog": {
            "graph_dims": list(GRAPH_DIMS),
            "llm_alias": LLM_DIM_ALIAS,
            "style_subfamilies": {k: list(v) for k, v in STYLE_SUBFAMILIES.items()},
        },
        "style_dims": style_attr,
        "graph_dims": graph_attr,
        "audit_dims": audit_attr,
    }
    (OUT / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2), encoding="utf-8"
    )
    write_report(summary, OUT / "REPORT.md")
    print(f"[style] top={style_attr['top_dims']} mass={style_attr['rf_mass']}")
    print(f"[graph] top={graph_attr['top_dims']} logo={graph_attr['logo_share']}")
    print(f"[audit] top={audit_attr['top_dims']} logo={audit_attr['logo_share']}")
    print(f"wrote {OUT / 'REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
