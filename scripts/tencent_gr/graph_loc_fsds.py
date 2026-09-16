#!/usr/bin/env python3
"""TencentGR: graph localization → FSDS, with inject Hit@K evaluation.

Uses ``results/tencent_gr_fs150/user_feats_selected.parquet`` (150 feats).

Graph dimensions (business-shaped nodes on this board)::

  user     ← funnel + session + decay   (who / history)
  merchant ← money / price / pay         (supply / ARPU side)
  product  ← diversity / item entropy    (SKU mix)
  order    ← attr + markov               (path / transition)

L1: dimension mass (RF-domain + LOGO) + leave-one-dim-out ΔAUC
L2: FSDS/LOGO only inside S* = Top-2 dims
Eval: inject a mean-shift hop on one dim in the late clock half; report Hit@K
      vs flat all-feature FSDS mapped back to dims.

Usage::

  PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path("results/tencent_gr_fs150")
OUT = ROOT / "graph_loc_fsds"
SEED = 0
TREES = 25
DEPTH = 5
FOLDS = 2
TOP_K = 2


# ---------------------------------------------------------------------------
# Feature → graph dimension
# ---------------------------------------------------------------------------


def family_of(name: str) -> str:
    if name.startswith("attr_"):
        return "attr"
    if name.startswith("trans"):
        return "markov"
    if name.startswith("sess_"):
        return "session"
    if name.startswith("dec_"):
        return "decay"
    if name.startswith("x_") or name.startswith("sq_") or name.startswith("log1p_abs_"):
        return "cross"
    if name.startswith("item_entropy") or name.startswith("n_uniq_"):
        return "diversity"
    money = (
        "pay_cnt",
        "log1p_pay_cnt",
        "pay_user",
        "arpu_sum_proxy",
        "arpu_mean_proxy",
        "arpu_p50_proxy",
        "log1p_arpu_sum",
        "arpu_per_active_day_proxy",
        "life_cnv_price_mean",
        "life_cnv_price_sum",
        "30d_cnv_price_mean",
    )
    if name in money or "price" in name or name.startswith("arpu_"):
        return "money"
    return "funnel"


def dim_of(name: str) -> str | None:
    """Map feature family → graph node. Cross is edge-like; excluded from L1 nodes."""
    fam = family_of(name)
    if fam in ("funnel", "session", "decay"):
        return "user"
    if fam == "money":
        return "merchant"
    if fam == "diversity":
        return "product"
    if fam in ("attr", "markov"):
        return "order"
    return None  # cross etc.


DIMS = ("user", "merchant", "product", "order")


# ---------------------------------------------------------------------------
# FSDS-lite
# ---------------------------------------------------------------------------


def rf_domain(X: np.ndarray, W: np.ndarray, *, seed: int) -> tuple[float, np.ndarray]:
    clf = RandomForestClassifier(
        n_estimators=TREES,
        max_depth=DEPTH,
        min_samples_leaf=4,
        random_state=seed,
        n_jobs=-1,
    )
    n = len(W)
    idx = np.arange(n)
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    cut = n * 3 // 4
    tr, te = idx[:cut], idx[cut:]
    clf.fit(X[tr], W[tr])
    auc = float(roc_auc_score(W[te], clf.predict_proba(X[te])[:, 1]))
    return auc, clf.feature_importances_.astype(float)


def po_risk_fit(X: np.ndarray, Y: np.ndarray, W: np.ndarray, *, seed: int) -> dict:
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    cv = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=TREES,
            max_depth=DEPTH,
            min_samples_leaf=4,
            random_state=seed + fold,
            n_jobs=-1,
        )
        e = RandomForestClassifier(
            n_estimators=TREES,
            max_depth=DEPTH,
            min_samples_leaf=4,
            random_state=seed + 40 + fold,
            n_jobs=-1,
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, 0.05, 0.95)
    po = (Yf - m_hat) * (W.astype(float) - e_hat)
    tau = RandomForestRegressor(
        n_estimators=TREES,
        max_depth=DEPTH,
        min_samples_leaf=4,
        random_state=seed + 7,
        n_jobs=-1,
    )
    tau.fit(X, po)
    tau_hat = tau.predict(X)
    return {
        "risk": float(np.mean(tau_hat**2)),
        "vimp": tau.feature_importances_.astype(float),
        "po": po,
    }


def mass_by_groups(vimp: np.ndarray, groups: dict[str, list[int]]) -> dict[str, float]:
    out = {g: float(np.sum(vimp[ix])) if ix else 0.0 for g, ix in groups.items()}
    s = sum(out.values()) + 1e-12
    return {k: out[k] / s for k in out}


def logo_by_groups(
    X: np.ndarray, Y: np.ndarray, W: np.ndarray, groups: dict[str, list[int]], *, seed: int
) -> dict:
    full = po_risk_fit(X, Y, W, seed=seed)
    logo = {}
    for i, (g, ix) in enumerate(groups.items()):
        drop = set(ix)
        keep = [j for j in range(X.shape[1]) if j not in drop]
        if len(keep) < 2:
            logo[g] = {"delta": 0.0, "R_minus": full["risk"]}
            continue
        rm = po_risk_fit(X[:, keep], Y, W, seed=seed + 11 + i)["risk"]
        logo[g] = {"delta": float(full["risk"] - rm), "R_minus": float(rm)}
    pos = {g: max(v["delta"], 0.0) for g, v in logo.items()}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in groups}
    return {"po_risk": full["risk"], "logo": logo, "logo_share": share, "po_vimp": full["vimp"]}


def leave_one_dim_delta_auc(
    X: np.ndarray, W: np.ndarray, groups: dict[str, list[int]], *, seed: int
) -> dict[str, float]:
    """How much RF-domain AUC drops when dim g is removed (larger = more responsible)."""
    full_auc, _ = rf_domain(X, W, seed=seed)
    out = {}
    for i, (g, ix) in enumerate(groups.items()):
        drop = set(ix)
        keep = [j for j in range(X.shape[1]) if j not in drop]
        if len(keep) < 2:
            out[g] = 0.0
            continue
        auc_m, _ = rf_domain(X[:, keep], W, seed=seed + 3 + i)
        out[g] = float(full_auc - auc_m)
    # shift to non-neg for ranking share
    m = min(out.values()) if out else 0.0
    shifted = {k: max(v - m, 0.0) for k, v in out.items()}
    s = sum(shifted.values()) + 1e-12
    return {k: shifted[k] / s for k in shifted}


def combine_l1_scores(
    rf_mass: dict[str, float],
    logo_share: dict[str, float],
    lodo: dict[str, float],
) -> dict[str, float]:
    """Equal blend of three L1 signals → graph localization score."""
    keys = list(rf_mass)
    raw = {
        k: (rf_mass.get(k, 0.0) + logo_share.get(k, 0.0) + lodo.get(k, 0.0)) / 3.0
        for k in keys
    }
    s = sum(raw.values()) + 1e-12
    return {k: raw[k] / s for k in keys}


def top_k(scores: dict[str, float], k: int) -> list[str]:
    return [d for d, _ in sorted(scores.items(), key=lambda kv: -kv[1])[:k]]


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def load_board() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], dict]:
    df = pd.read_parquet(ROOT / "user_feats_selected.parquet")
    fs = json.loads((ROOT / "selected_features.json").read_text())
    names = [x["name"] for x in fs if x["name"] in df.columns]
    # Prefer graph dims; keep cross for flat baseline only
    graph_names = [n for n in names if dim_of(n) is not None]
    X = df[graph_names].fillna(0).to_numpy(np.float64)
    Y = df["_label_future_cnv"].to_numpy(np.int64)
    clock = df["_seq_t_end"].to_numpy(np.float64)
    W = (clock > np.median(clock)).astype(int)
    groups: dict[str, list[int]] = {d: [] for d in DIMS}
    for j, n in enumerate(graph_names):
        groups[dim_of(n)].append(j)
    meta = {
        "n": int(len(df)),
        "p_graph": len(graph_names),
        "p_all_selected": len(names),
        "dim_sizes": {d: len(groups[d]) for d in DIMS},
        "pos_rate": float(Y.mean()),
        "w1_rate": float(W.mean()),
    }
    return X, Y, W, graph_names, {"groups": groups, "meta": meta, "all_names": names, "df": df}


def run_l1(X, Y, W, groups, *, seed: int, heavy: bool = False) -> dict:
    auc, v_rf = rf_domain(X, W, seed=seed)
    rf_mass = mass_by_groups(v_rf, groups)
    lodo = leave_one_dim_delta_auc(X, W, groups, seed=seed)
    if heavy:
        logo = logo_by_groups(X, Y, W, groups, seed=seed)
        logo_share = logo["logo_share"]
        po_risk = logo["po_risk"]
        logo_raw = logo["logo"]
    else:
        # Fast path: use RF mass as PO proxy for ranking (full LOGO on natural board only)
        logo_share = rf_mass
        po_risk = float("nan")
        logo_raw = {}
    blended = combine_l1_scores(rf_mass, logo_share, lodo)
    s_star = top_k(blended, TOP_K)
    return {
        "rf_domain_auc": auc,
        "rf_mass": rf_mass,
        "logo_share": logo_share,
        "logo": logo_raw,
        "lodo_share": lodo,
        "blended": blended,
        "S_star": s_star,
        "po_risk": po_risk,
    }


def run_l2(X, Y, W, names, groups, s_star: list[str], *, seed: int) -> dict:
    """FSDS inside S* only."""
    ix = []
    local_groups: dict[str, list[int]] = {}
    local_names = []
    for d in s_star:
        start = len(ix)
        for j in groups[d]:
            ix.append(j)
            local_names.append(names[j])
        local_groups[d] = list(range(start, len(ix)))
    if len(ix) < 2:
        return {"empty": True}
    Xs = X[:, ix]
    auc, v_rf = rf_domain(Xs, W, seed=seed)
    logo = logo_by_groups(Xs, Y, W, local_groups, seed=seed + 5)
    # within-dim feature ranks by RF vimp
    within = {}
    for d in s_star:
        sub = local_groups[d]
        local_v = v_rf[sub]
        local_n = [local_names[j] for j in sub]
        order = np.argsort(-local_v)
        within[d] = [
            {"rank": r + 1, "name": local_n[j], "vimp": float(local_v[j])}
            for r, j in enumerate(order[:8])
        ]
    return {
        "rf_domain_auc": auc,
        "rf_mass": mass_by_groups(v_rf, local_groups),
        "logo_share": logo["logo_share"],
        "within_dim": within,
        "po_risk": logo["po_risk"],
    }


def flat_fsds_dim_scores(X, Y, W, names, *, seed: int) -> dict[str, float]:
    """Baseline: all-feature RF-domain mass aggregated to dims (no LODO/LOGO)."""
    groups = {d: [] for d in DIMS}
    for j, n in enumerate(names):
        d = dim_of(n)
        if d:
            groups[d].append(j)
    auc, v_rf = rf_domain(X, W, seed=seed)
    _ = auc
    rf_mass = mass_by_groups(v_rf, groups)
    return rf_mass


def inject_dim_hop(
    X: np.ndarray,
    W: np.ndarray,
    groups: dict[str, list[int]],
    dim: str,
    *,
    strength: float,
    seed: int,
    confound: float = 0.0,
) -> np.ndarray:
    """Shift late-window columns of ``dim``; optional small noise on other dims."""
    rng = np.random.default_rng(seed)
    X2 = X.copy()
    late = W == 1
    for g, ix in groups.items():
        if not ix:
            continue
        for j in ix:
            col = X2[:, j]
            sd = float(np.std(col)) + 1e-6
            if g == dim:
                X2[late, j] = (
                    col[late]
                    + strength * sd
                    + rng.normal(0, 0.05 * sd, int(late.sum()))
                )
            elif confound > 0:
                X2[late, j] = col[late] + rng.normal(
                    0, confound * sd, int(late.sum())
                )
    return X2


def eval_inject(
    X,
    Y,
    W,
    names,
    groups,
    *,
    strengths=(0.8, 1.5, 2.5),
    confounds=(0.0, 0.35),
    seed: int = 0,
) -> dict:
    rows = []
    for dim in DIMS:
        if len(groups[dim]) < 1:
            continue
        for strength in strengths:
            for confound in confounds:
                Xh = inject_dim_hop(
                    X,
                    W,
                    groups,
                    dim,
                    strength=strength,
                    seed=seed + hash((dim, strength, confound)) % 10000,
                    confound=confound,
                )
                l1 = run_l1(Xh, Y, W, groups, seed=seed, heavy=False)
                flat = flat_fsds_dim_scores(Xh, Y, W, names, seed=seed)
                s_star = l1["S_star"]
                flat_top = top_k(flat, TOP_K)
                flat_top1 = top_k(flat, 1)
                graph_top1 = top_k(l1["blended"], 1)
                hit_graph = dim in s_star
                hit_flat = dim in flat_top
                hit_graph1 = dim in graph_top1
                hit_flat1 = dim in flat_top1
                rank_graph = (
                    sorted(l1["blended"], key=lambda k: -l1["blended"][k]).index(dim) + 1
                )
                rank_flat = sorted(flat, key=lambda k: -flat[k]).index(dim) + 1
                # Skip L2 inside inject loop for speed; run once on natural board
                l2 = {}
                rows.append(
                    {
                        "injected_dim": dim,
                        "strength": strength,
                        "confound": confound,
                        "S_star": s_star,
                        "flat_top": flat_top,
                        "hit_graph@2": hit_graph,
                        "hit_flat@2": hit_flat,
                        "hit_graph@1": hit_graph1,
                        "hit_flat@1": hit_flat1,
                        "rank_graph": rank_graph,
                        "rank_flat": rank_flat,
                        "l1_blended": l1["blended"],
                        "flat_blended": flat,
                        "l2_within_top": {
                            d: (l2.get("within_dim") or {}).get(d, [])[:3] for d in s_star
                        }
                        if l2
                        else {},
                    }
                )
    g_hits = [r["hit_graph@2"] for r in rows]
    f_hits = [r["hit_flat@2"] for r in rows]
    g1 = [r["hit_graph@1"] for r in rows]
    f1 = [r["hit_flat@1"] for r in rows]
    # hard slice: weak + confound
    hard = [r for r in rows if r["strength"] <= 1.0 and r["confound"] > 0]
    summary = {
        "n_trials": len(rows),
        "hit_at_2_graph": float(np.mean(g_hits)) if g_hits else 0.0,
        "hit_at_2_flat": float(np.mean(f_hits)) if f_hits else 0.0,
        "hit_at_1_graph": float(np.mean(g1)) if g1 else 0.0,
        "hit_at_1_flat": float(np.mean(f1)) if f1 else 0.0,
        "mean_rank_graph": float(np.mean([r["rank_graph"] for r in rows])) if rows else None,
        "mean_rank_flat": float(np.mean([r["rank_flat"] for r in rows])) if rows else None,
        "graph_beats_or_ties_flat": float(
            np.mean([r["rank_graph"] <= r["rank_flat"] for r in rows])
        )
        if rows
        else None,
        "hard_n": len(hard),
        "hard_hit_at_2_graph": float(np.mean([r["hit_graph@2"] for r in hard])) if hard else None,
        "hard_hit_at_2_flat": float(np.mean([r["hit_flat@2"] for r in hard])) if hard else None,
        "hard_hit_at_1_graph": float(np.mean([r["hit_graph@1"] for r in hard])) if hard else None,
        "hard_hit_at_1_flat": float(np.mean([r["hit_flat@1"] for r in hard])) if hard else None,
    }
    return {"trials": rows, "summary": summary}


def write_report(natural: dict, inject: dict, meta: dict) -> str:
    s = inject["summary"]
    trial_lines = []
    for r in inject["trials"]:
        trial_lines.append(
            f"| {r['injected_dim']} | {r['strength']} | {r.get('confound', 0)} | "
            f"`{r['S_star']}` | `{r['flat_top']}` | "
            f"{'Y' if r['hit_graph@2'] else 'N'} | {'Y' if r['hit_flat@2'] else 'N'} | "
            f"{'Y' if r.get('hit_graph@1') else 'N'} | {'Y' if r.get('hit_flat@1') else 'N'} | "
            f"{r['rank_graph']} | {r['rank_flat']} |"
        )
    nat = natural
    md = f"""# TencentGR：Graph localization → FSDS

数据：`user_feats_selected.parquet` n={meta['n']} p_graph={meta['p_graph']}  
维大小：`{meta['dim_sizes']}`  Y=future_cnv  W=clock median

## 句式

```text
L1 图维 localization（user / merchant / product / order）
  → S* = Top-{TOP_K}
L2 仅在 S* 内 FSDS / LOGO / 特征排序
```

## Natural board（无注入）

- RF-domain AUC: **{nat['rf_domain_auc']:.3f}**
- S*: **{nat['S_star']}**
- blended mass: `{ {k: round(v,3) for k,v in nat['blended'].items()} }`
- LOGO share: `{ {k: round(v,3) for k,v in nat['logo_share'].items()} }`

## Inject Hit@K（精确归因主考）

late 半窗平移目标维；`confound>0` 时其它维加干扰噪声。

| 注入维 | strength | confound | Graph S* | Flat Top-2 | H@2_g | H@2_f | H@1_g | H@1_f | rank_g | rank_f |
|--------|----------|----------|----------|------------|-------|-------|-------|-------|--------|--------|
{chr(10).join(trial_lines)}

### Summary

| 指标 | 值 |
|------|-----|
| trials | {s['n_trials']} |
| **Hit@2 graph** | **{s['hit_at_2_graph']:.3f}** |
| Hit@2 flat FSDS | {s['hit_at_2_flat']:.3f} |
| **Hit@1 graph** | **{s['hit_at_1_graph']:.3f}** |
| Hit@1 flat | {s['hit_at_1_flat']:.3f} |
| mean rank graph / flat | {s['mean_rank_graph']:.2f} / {s['mean_rank_flat']:.2f} |
| hard (weak+confound) Hit@2 g/f | {s.get('hard_hit_at_2_graph')} / {s.get('hard_hit_at_2_flat')} |
| hard Hit@1 g/f | {s.get('hard_hit_at_1_graph')} / {s.get('hard_hit_at_1_flat')} |

对外一句：注入维 Hit@2（graph）= **{s['hit_at_2_graph']:.0%}**，Hit@1= **{s['hit_at_1_graph']:.0%}**
（flat @2/@1 = {s['hit_at_2_flat']:.0%}/{s['hit_at_1_flat']:.0%}）；
自然窗 S*={nat['S_star']}。

## 读法

- Hit@K 高 ⇒ 图维 L1 **能点到被污染的维**（维级精确归因）。
- hard 切片（弱注入+干扰）更接近线上噪窗。
- L2 只在 S* 内排特征；不宣称全空间唯一因果。

```bash
PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds.py
```
"""
    return md


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    X, Y, W, names, pack = load_board()
    groups, meta = pack["groups"], pack["meta"]
    print(f"n={meta['n']} p={meta['p_graph']} dims={meta['dim_sizes']} pos={meta['pos_rate']:.3f}")

    natural = run_l1(X, Y, W, groups, seed=SEED, heavy=True)
    natural_l2 = run_l2(X, Y, W, names, groups, natural["S_star"], seed=SEED)
    print(f"natural S*={natural['S_star']} auc={natural['rf_domain_auc']:.3f}")

    inject = eval_inject(
        X, Y, W, names, groups, strengths=(0.8, 1.5, 2.5), confounds=(0.0, 0.35), seed=SEED
    )
    s = inject["summary"]
    print(
        f"inject Hit@2 graph={s['hit_at_2_graph']:.3f} flat={s['hit_at_2_flat']:.3f} "
        f"Hit@1 g/f={s['hit_at_1_graph']:.3f}/{s['hit_at_1_flat']:.3f} "
        f"hard@2 g/f={s.get('hard_hit_at_2_graph')}/{s.get('hard_hit_at_2_flat')}"
    )

    summary = {
        "stance": "graph localization → FSDS on TencentGR 150; inject Hit@K for precise dim attribution",
        "dims": list(DIMS),
        "meta": meta,
        "natural_l1": natural,
        "natural_l2": natural_l2,
        "inject": inject,
        "external_one_liner_cn": (
            f"腾讯GR图定位→FSDS：注入维 Hit@2={s['hit_at_2_graph']:.0%}"
            f"/Hit@1={s['hit_at_1_graph']:.0%}"
            f"（flat @2/@1={s['hit_at_2_flat']:.0%}/{s['hit_at_1_flat']:.0%}）；"
            f"hard弱+扰 Hit@2={s.get('hard_hit_at_2_graph')}；"
            f"自然窗 S*={natural['S_star']}。"
        ),
    }
    (OUT / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2), encoding="utf-8"
    )
    (OUT / "natural_l2.json").write_text(
        json.dumps(jsonable(natural_l2), indent=2), encoding="utf-8"
    )
    md = write_report(natural, inject, meta)
    (OUT / "REPORT.md").write_text(md, encoding="utf-8")
    # also drop a short biz pointer
    docs = Path("docs/biz/TENCENT_GRAPH_LOC_FSDS.md")
    docs.write_text(
        md
        + "\n\n评估协议对照：`GRAPH_USECASE_EVAL.md`。\n",
        encoding="utf-8",
    )
    print(summary["external_one_liner_cn"])
    print(f"wrote {OUT / 'REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
