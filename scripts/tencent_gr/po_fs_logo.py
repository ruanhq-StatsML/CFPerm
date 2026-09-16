#!/usr/bin/env python3
"""PO-risk / FSDS feature-selection board on the TencentGR 150.

Clock = prefix last-ts. Batch W = 1 if later than median clock.
Y = future_cnv. F-score ranking is in-sample association with Y.
This board is OOD VIMP on (X, Y, W):

  RF-domain  : P(W|X) VIMP          — covariate / P(X)
  PO-risk    : R = E[τ̂(X)²], τ̂≈E[φ|X], φ=(Y-μ)(W-e)  — concept-relevant
  LOGO family: Δ_g = R − R^{(-g)}

  PYTHONPATH=. python3 scripts/tencent_gr/po_fs_logo.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path("results/tencent_gr_fs150")
SEED = 0
TREES = 40
DEPTH = 6
FOLDS = 3


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


def po_risk_fit(X, Y, W, *, seed: int):
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


def rf_domain(X, W, *, seed: int):
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


def mass(vimp, names, families):
    out = {g: 0.0 for g in families}
    for v, n in zip(vimp, names):
        out[family_of(n)] += float(v)
    s = sum(out.values()) + 1e-12
    return {k: out[k] / s for k in families}


def main() -> None:
    df = pd.read_parquet(ROOT / "user_feats_selected.parquet")
    fs = json.loads((ROOT / "selected_features.json").read_text())
    fscore = {x["name"]: float(x["score"]) for x in fs}
    names = [x["name"] for x in fs]
    drop = {"user_id", "_seq_t_end", "_label_future_cnv", "_label_pay_user"}
    # keep F-score order, one copy
    cols = [c for c in names if c in df.columns]
    X = df[cols].fillna(0).to_numpy(np.float64)
    Y = df["_label_future_cnv"].to_numpy(np.int64)
    clock = df["_seq_t_end"].to_numpy(np.float64)
    W = (clock > np.median(clock)).astype(int)
    families = sorted({family_of(n) for n in cols})
    groups = {g: [i for i, n in enumerate(cols) if family_of(n) == g] for g in families}

    print(f"n={len(df)} p={len(cols)} pos={Y.mean():.3f} W1={W.mean():.3f}")
    print("family sizes", {g: len(ix) for g, ix in groups.items()})

    auc, v_rf = rf_domain(X, W, seed=SEED)
    print(f"RF-domain AUC={auc:.3f}")
    full = po_risk_fit(X, Y, W, seed=SEED)
    print(f"PO-risk={full['risk']:.6f}")

    logo = {}
    for g, ix in groups.items():
        keep = [j for j in range(X.shape[1]) if j not in set(ix)]
        rm = po_risk_fit(X[:, keep], Y, W, seed=SEED + 11 + families.index(g))["risk"]
        logo[g] = {"R_minus": rm, "delta": float(full["risk"] - rm)}
        print(f"  LOGO {g:10s} Δ={logo[g]['delta']:+.6f}  R(-g)={rm:.6f}")

    pos = {g: max(logo[g]["delta"], 0.0) for g in families}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in families}

    po_v = full["vimp"]
    rf_v = v_rf
    f_v = np.array([fscore.get(n, 0.0) for n in cols])
    tau_po_f = float(kendalltau(po_v, f_v).correlation)
    tau_rf_f = float(kendalltau(rf_v, f_v).correlation)
    tau_po_rf = float(kendalltau(po_v, rf_v).correlation)

    order_po = np.argsort(-po_v)
    top_po = [{"rank": r + 1, "name": cols[j], "family": family_of(cols[j]), "po_vimp": float(po_v[j]), "fscore": float(f_v[j])} for r, j in enumerate(order_po[:20])]
    order_f = np.argsort(-f_v)
    top_f = [{"rank": r + 1, "name": cols[j], "family": family_of(cols[j])} for r, j in enumerate(order_f[:20])]

    summary = {
        "n": int(len(df)),
        "p": int(len(cols)),
        "pos_rate": float(Y.mean()),
        "w1_rate": float(W.mean()),
        "clock": "seq_t_end median split",
        "rf_domain_auc": auc,
        "po_risk": full["risk"],
        "family_n": {g: len(ix) for g, ix in groups.items()},
        "rf_mass": mass(rf_v, cols, families),
        "po_mass": mass(po_v, cols, families),
        "logo": logo,
        "logo_share": share,
        "kendall": {
            "po_vs_fscore": tau_po_f,
            "rfdomain_vs_fscore": tau_rf_f,
            "po_vs_rfdomain": tau_po_rf,
        },
        "top20_po": top_po,
        "top20_fscore": top_f,
    }
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT / "PO_FS_LOGO.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    def fmt(d):
        return " ".join(f"{k}={v:.2f}" for k, v in sorted(d.items(), key=lambda kv: -kv[1]))

    md = f"""# PO-risk feature selection vs F-score (TencentGR 150)

Clock = `_seq_t_end` median. W=later half. Y=`future_cnv`.
RF-domain AUC **{auc:.3f}**. PO-risk R=E[τ̂²] **{full['risk']:.6f}**.

Kendall τ: PO-VIMP vs F-score **{tau_po_f:.3f}**; RF-domain vs F-score **{tau_rf_f:.3f}**; PO vs RF-domain **{tau_po_rf:.3f}**.

## Family mass

| family | n | RF-domain (P(X)) | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
"""
    for g in sorted(families, key=lambda x: -share[x]):
        md += (
            f"| {g} | {len(groups[g])} | {summary['rf_mass'][g]:.3f} | "
            f"{summary['po_mass'][g]:.3f} | {logo[g]['delta']:+.5f} | {share[g]:.3f} |\n"
        )
    md += (
        "\nF-score ranks in-sample Y association. PO-LOGO ranks which "
        "family carries the *batch map* φ=(Y-μ)(W-e). Low Kendall ⇒ the 150 "
        "F-board is not an OOD board.\n"
    )
    (ROOT / "PO_FS_LOGO.md").write_text(md, encoding="utf-8")
    print(json.dumps({"kendall": summary["kendall"], "logo_share": share, "rf_auc": auc}, indent=2))
    print("wrote", ROOT / "PO_FS_LOGO.md")


if __name__ == "__main__":
    main()
