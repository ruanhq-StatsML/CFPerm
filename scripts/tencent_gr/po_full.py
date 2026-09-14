#!/usr/bin/env python3
"""Direct PO-risk VIMP on the full one-pass user matrix (no X_KEEP, no F-cut).

X = tables/user.parquet 全部数值列（one-pass 漏斗/衰减/场/归因/买后/转移/cross）。
Y = future_cnv（prefix 外）。W = 1{_seq_t_end > median}。
不做 LOGO、不做 Y-swap、不做 bootstrap。φ=(Y-μ)(W-e)，VIMP 全列排序。

顺带同一套 (Y,W) 打 F-score 150，方便对照。

  PYTHONPATH=. python3 scripts/tencent_gr/po_full.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from po_fs_logo import family_of, po_risk_fit, rf_domain  # noqa: E402

ROOT = Path(__file__).resolve().parents[2] / "results/tencent_gr_fs150"
DROP = {"user_id", "_seq_t_end", "_label_future_cnv", "_label_pay_user"}


def _xy(df: pd.DataFrame, cols: list[str]):
    X = df[cols].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(np.float64)
    Y = df["_label_future_cnv"].to_numpy(np.int64)
    clock = df["_seq_t_end"].to_numpy(np.float64)
    W = (clock > np.median(clock)).astype(int)
    return X, Y, W


def _board(name: str, df: pd.DataFrame, cols: list[str], *, seed: int) -> dict:
    X, Y, W = _xy(df, cols)
    print(f"{name} n={len(df)} p={len(cols)} pos={Y.mean():.3f} W1={W.mean():.3f}", flush=True)
    auc, v_rf = rf_domain(X, W, seed=seed)
    print(f"  RF-domain AUC={auc:.3f}", flush=True)
    rec = po_risk_fit(X, Y, W, seed=seed)
    print(f"  PO-risk={rec['risk']:.6f}", flush=True)
    po_v = rec["vimp"]
    rf_v = v_rf
    ranked = []
    for j in np.argsort(-po_v):
        ranked.append(
            {
                "rank": len(ranked) + 1,
                "name": cols[j],
                "family": family_of(cols[j]),
                "po_vimp": float(po_v[j]),
                "rf_vimp": float(rf_v[j]),
            }
        )
    fam = sorted({r["family"] for r in ranked})
    po_mass, rf_mass = {}, {}
    for g in fam:
        po_mass[g] = float(sum(r["po_vimp"] for r in ranked if r["family"] == g))
        rf_mass[g] = float(sum(r["rf_vimp"] for r in ranked if r["family"] == g))
    s_po = sum(po_mass.values()) + 1e-12
    s_rf = sum(rf_mass.values()) + 1e-12
    return {
        "name": name,
        "n": int(len(df)),
        "p": int(len(cols)),
        "pos_rate": float(Y.mean()),
        "w1_rate": float(W.mean()),
        "rf_domain_auc": float(auc),
        "po_risk": float(rec["risk"]),
        "po_mass": {k: v / s_po for k, v in po_mass.items()},
        "rf_mass": {k: v / s_rf for k, v in rf_mass.items()},
        "ranked": ranked,
    }


def _md(boards: list[dict]) -> str:
    lines = [
        "# Direct PO on full one-pass X",
        "",
        "Y=`future_cnv`. W=`1{t_end > median}`. 不裁列、不换 Y、不做 LOGO。",
        "",
    ]
    for b in boards:
        lines += [
            f"## `{b['name']}`  n={b['n']} p={b['p']}  pos={b['pos_rate']:.3f}  "
            f"RF-domain **{b['rf_domain_auc']:.3f}**  PO-risk **{b['po_risk']:.6f}**",
            "",
            "family mass | RF-domain (P(X)) | PO-VIMP |",
            "|---|---:|---:|",
        ]
        fams = sorted(b["po_mass"], key=lambda g: -b["po_mass"][g])
        for g in fams:
            lines.append(f"| {g} | {b['rf_mass'][g]:.3f} | {b['po_mass'][g]:.3f} |")
        lines += [
            "",
            "| rank | feat | family | PO-VIMP | RF-VIMP |",
            "|---|---|---|---:|---:|",
        ]
        for r in b["ranked"]:
            lines.append(
                f"| {r['rank']} | `{r['name']}` | {r['family']} | "
                f"{r['po_vimp']:.4f} | {r['rf_vimp']:.4f} |"
            )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    lab = pd.read_parquet(ROOT / "user_feats_selected.parquet")[
        ["user_id", "_seq_t_end", "_label_future_cnv"]
    ]
    user = pd.read_parquet(ROOT / "tables/user.parquet").merge(lab, on="user_id", how="inner")
    onepass_cols = [
        c
        for c in user.columns
        if c not in DROP and np.issubdtype(user[c].dtype, np.number)
    ]
    fs = json.loads((ROOT / "selected_features.json").read_text())
    f150 = pd.read_parquet(ROOT / "user_feats_selected.parquet")
    f150_cols = [x["name"] for x in fs if x["name"] in f150.columns]

    boards = [
        _board("onepass_full", user, onepass_cols, seed=0),
        _board("f150", f150, f150_cols, seed=0),
    ]
    payload = {"clock": "seq_t_end median", "y": "future_cnv", "boards": boards}
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT / "PO_FULL.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (ROOT / "PO_FULL.md").write_text(_md(boards), encoding="utf-8")
    (ROOT / "PO_FULL.txt").write_text(
        "\n".join(
            f"# {b['name']} p={b['p']} AUC={b['rf_domain_auc']:.3f} R={b['po_risk']:.6f}\n"
            + "\n".join(
                f"{r['rank']:3d}  {r['name']:42s}  {r['family']:10s}  "
                f"{r['po_vimp']:.4f}  rf={r['rf_vimp']:.4f}"
                for r in b["ranked"]
            )
            for b in boards
        )
        + "\n",
        encoding="utf-8",
    )
    for b in boards:
        print(f"\n=== {b['name']} top20 ===")
        for r in b["ranked"][:20]:
            print(
                f"{r['rank']:3d}  {r['name']:42s}  {r['family']:10s}  {r['po_vimp']:.4f}"
            )
    print("wrote", ROOT / "PO_FULL.md")


if __name__ == "__main__":
    main()
