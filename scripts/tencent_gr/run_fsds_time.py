#!/usr/bin/env python3
"""FSDS on TencentGR 150, time clocks only.

Clocks (no shuffle):
  median   W = 1{seq_t_end > median}
  q1q4     W on Q1 vs Q4 tails
Holdout: first 75% by clock → last 25% (HGB).

Boards: RF-domain (P(W|X)), PO-risk LOGO, time-OOS HGB
(F-150 vs PO-32 vs RF-domain-32).

  PYTHONPATH=. python3 scripts/tencent_gr/run_fsds_time.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score

sys.path.insert(0, str(Path("scripts/tencent_gr").resolve()))
from po_fs_logo import family_of, mass, po_risk_fit, rf_domain

ROOT = Path("results/tencent_gr_fs150")
SEED = 0
K = 32


def load():
    df = pd.read_parquet(ROOT / "user_feats_selected.parquet")
    fs = json.loads((ROOT / "selected_features.json").read_text())
    po32 = json.loads((ROOT / "FEATURES_PO32.json").read_text())
    names = [x["name"] for x in fs if x["name"] in df.columns]
    X = df[names].fillna(0).to_numpy(np.float64)
    y = df["_label_future_cnv"].to_numpy(np.int64)
    clock = df["_seq_t_end"].to_numpy(np.float64)
    fscore = np.array([next(x["score"] for x in fs if x["name"] == n) for n in names])
    return df, names, X, y, clock, fscore, [x["name"] for x in po32]


def clock_w(clock, kind: str):
    if kind == "median":
        mask = np.ones(len(clock), dtype=bool)
        W = (clock > np.median(clock)).astype(int)
        return mask, W
    q1, q3 = np.quantile(clock, [0.25, 0.75])
    mask = (clock <= q1) | (clock >= q3)
    W = np.zeros(len(clock), dtype=int)
    W[clock >= q3] = 1
    return mask, W


def hgb_split(Xtr, ytr, Xte, yte):
    if len(np.unique(yte)) < 2:
        return {
            "auc": float("nan"),
            "ap": float(average_precision_score(yte, np.zeros(len(yte)))) if len(yte) else float("nan"),
            "acc": float(accuracy_score(yte, np.zeros(len(yte), dtype=int))) if len(yte) else float("nan"),
            "pos_tr": float(ytr.mean()),
            "pos_te": float(yte.mean()),
            "n_tr": int(len(ytr)),
            "n_te": int(len(yte)),
            "note": "single-class test",
        }
    m = HistGradientBoostingClassifier(
        max_depth=4, max_iter=120, learning_rate=0.08, random_state=42
    )
    m.fit(Xtr, ytr)
    p = m.predict_proba(Xte)[:, 1]
    return {
        "auc": float(roc_auc_score(yte, p)),
        "ap": float(average_precision_score(yte, p)),
        "acc": float(accuracy_score(yte, (p >= 0.5).astype(int))),
        "pos_tr": float(ytr.mean()),
        "pos_te": float(yte.mean()),
        "n_tr": int(len(ytr)),
        "n_te": int(len(yte)),
    }


def col_take(X, cols_idx):
    return X if cols_idx is None else X[:, cols_idx]


def logo_families(X, y, W, names, seed):
    families = sorted({family_of(n) for n in names})
    groups = {g: [i for i, n in enumerate(names) if family_of(n) == g] for g in families}
    full = po_risk_fit(X, y, W, seed=seed)
    logo = {}
    for g, ix in groups.items():
        keep = [j for j in range(X.shape[1]) if j not in set(ix)]
        rm = po_risk_fit(X[:, keep], y, W, seed=seed + 11 + families.index(g))["risk"]
        logo[g] = {"delta": float(full["risk"] - rm), "R_minus": rm, "n": len(ix)}
    pos = {g: max(v["delta"], 0.0) for g, v in logo.items()}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in families}
    return full, logo, share, families


def board_clock(kind, names, X, y, clock, fscore, seed):
    mask, Wfull = clock_w(clock, kind)
    Xm, ym, Wm = X[mask], y[mask], Wfull[mask]
    auc, v_rf = rf_domain(Xm, Wm, seed=seed)
    full, logo, share, families = logo_families(Xm, ym, Wm, names, seed)
    tau_po_f = float(kendalltau(full["vimp"], fscore).correlation)
    tau_rf_f = float(kendalltau(v_rf, fscore).correlation)
    tau_po_rf = float(kendalltau(full["vimp"], v_rf).correlation)
    print(
        f"[{kind}] n={len(ym)} pos={ym.mean():.3f} W1={Wm.mean():.3f} "
        f"RF-AUC={auc:.3f} PO-R={full['risk']:.5f} τ(PO,F)={tau_po_f:.3f}"
    )
    return {
        "kind": kind,
        "n": int(len(ym)),
        "pos": float(ym.mean()),
        "w1": float(Wm.mean()),
        "rf_domain_auc": auc,
        "po_risk": full["risk"],
        "kendall": {
            "po_vs_fscore": tau_po_f,
            "rfdomain_vs_fscore": tau_rf_f,
            "po_vs_rfdomain": tau_po_rf,
        },
        "rf_mass": mass(v_rf, names, families),
        "po_mass": mass(full["vimp"], names, families),
        "logo": logo,
        "logo_share": share,
        "po_top": [
            {"rank": r + 1, "name": names[j], "family": family_of(names[j]), "vimp": float(full["vimp"][j])}
            for r, j in enumerate(np.argsort(-full["vimp"])[:10])
        ],
        "rf_top": [
            {"rank": r + 1, "name": names[j], "family": family_of(names[j]), "vimp": float(v_rf[j])}
            for r, j in enumerate(np.argsort(-v_rf)[:10])
        ],
        "rf_vimp": v_rf,
    }


def main() -> None:
    df, names, X, y, clock, fscore, po32_names = load()
    qs = np.quantile(clock, [0, 0.25, 0.5, 0.75, 1.0])
    quart = np.digitize(clock, qs[1:-1], right=True)
    pos_q = {f"Q{i+1}": float(y[quart == i].mean()) for i in range(4)}
    n_q = {f"Q{i+1}": int((quart == i).sum()) for i in range(4)}
    print("pos by time quartile", pos_q)

    clocks = {}
    for kind in ("median", "q1q4"):
        clocks[kind] = board_clock(kind, names, X, y, clock, fscore, SEED)

    idx = {n: i for i, n in enumerate(names)}
    po_idx = [idx[n] for n in po32_names if n in idx]
    rf_idx = list(np.argsort(-clocks["median"]["rf_vimp"])[:K])
    sets = {
        "Q1→Q2": (quart == 0, quart == 1),
        "Q12→Q3": (quart <= 1, quart == 2),
        "Q12→Q4": (quart <= 1, quart == 3),
    }
    hold = {}
    for split, (tr, te) in sets.items():
        hold[split] = {
            "f150": hgb_split(X[tr], y[tr], X[te], y[te]),
            "po32": hgb_split(col_take(X[tr], po_idx), y[tr], col_take(X[te], po_idx), y[te]),
            "rfdomain32": hgb_split(col_take(X[tr], rf_idx), y[tr], col_take(X[te], rf_idx), y[te]),
        }
        print(
            split,
            {k: (round(v["auc"], 4) if v["auc"] == v["auc"] else None) for k, v in hold[split].items()},
            "pos",
            hold[split]["f150"]["pos_tr"],
            "→",
            hold[split]["f150"]["pos_te"],
        )

    # drop vimp arrays from json
    for c in clocks.values():
        c.pop("rf_vimp", None)

    summary = {
        "clock": "_seq_t_end",
        "n": int(len(y)),
        "quartile_pos": pos_q,
        "quartile_n": n_q,
        "clocks": {k: {kk: vv for kk, vv in rec.items()} for k, rec in clocks.items()},
        "time_holdout_hgb": hold,
    }
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT / "FSDS_TIME.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    def row(kind):
        r = clocks[kind]
        return (
            f"| {kind} | {r['n']} | {r['pos']:.3f} | {r['rf_domain_auc']:.3f} | "
            f"{r['po_risk']:.5f} | {r['kendall']['po_vs_fscore']:.3f} |"
        )

    md = [
        "# FSDS time split (TencentGR 150)",
        "",
        "Clock = `_seq_t_end`. No shuffle.",
        f"Pos rate by quartile: {pos_q}",
        "",
        "| clock | n | pos | RF-domain AUC | PO-risk | τ(PO, F-score) |",
        "|---|---:|---:|---:|---:|---:|",
        row("median"),
        row("q1q4"),
        "",
        "Time holdout HGB (by quartile, not a dumb last-25% — Q4 has no positives):",
        "",
        "| split | set | AUC | AP | Acc | pos_tr | pos_te |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for split, recs in hold.items():
        for k, v in recs.items():
            auc = "nan" if v["auc"] != v["auc"] else f"{v['auc']:.4f}"
            md.append(
                f"| {split} | {k} | {auc} | {v['ap']:.4f} | {v['acc']:.4f} | "
                f"{v['pos_tr']:.3f} | {v['pos_te']:.3f} |"
            )
    md += [
        "",
        "RF-domain = covariate / P(X). PO-risk = φ=(Y-μ)(W-e).",
        "Q4 pos=0 is right-censoring on this clock, not a concept hop.",
        "Q1→Q2 is the comparable time OOS. Q12→Q3 is the label-rate drop.",
        "τ(PO, F) stays ~0: F board is not FSDS.",
        "",
    ]
    (ROOT / "FSDS_TIME.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote", ROOT / "FSDS_TIME.md")


if __name__ == "__main__":
    main()
