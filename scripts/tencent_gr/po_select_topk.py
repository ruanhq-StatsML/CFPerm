#!/usr/bin/env python3
"""PO-VIMP top-k (not F-score) on the 150 matrix, then HGB holdout."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split

sys.path.insert(0, str(Path("scripts/tencent_gr").resolve()))
from po_fs_logo import family_of, po_risk_fit

ROOT = Path("results/tencent_gr_fs150")
K = 32


def main() -> None:
    df = pd.read_parquet(ROOT / "user_feats_selected.parquet")
    fs = json.loads((ROOT / "selected_features.json").read_text())
    names = [x["name"] for x in fs if x["name"] in df.columns]
    X = df[names].fillna(0).to_numpy(np.float64)
    y = df["_label_future_cnv"].to_numpy(np.int64)
    clock = df["_seq_t_end"].to_numpy(np.float64)
    W = (clock > np.median(clock)).astype(int)
    rec = po_risk_fit(X, y, W, seed=0)
    order = np.argsort(-rec["vimp"])
    ranked = [
        {
            "rank": r + 1,
            "name": names[j],
            "family": family_of(names[j]),
            "po_vimp": float(rec["vimp"][j]),
        }
        for r, j in enumerate(order)
    ]
    top = ranked[:K]
    top_names = [x["name"] for x in top]
    (ROOT / "FEATURES_PO32.json").write_text(json.dumps(top, indent=2), encoding="utf-8")
    (ROOT / "FEATURES_PO32.txt").write_text(
        "\n".join(f"{x['rank']:2d}  {x['name']:42s}  {x['family']:10s}  {x['po_vimp']:.4f}" for x in top)
        + "\n",
        encoding="utf-8",
    )

    def holdout(cols):
        Xt = df[cols].fillna(0).to_numpy(np.float64)
        Xtr, Xte, ytr, yte = train_test_split(Xt, y, test_size=0.25, random_state=42, stratify=y)
        m = HistGradientBoostingClassifier(max_depth=4, max_iter=120, learning_rate=0.08, random_state=42)
        m.fit(Xtr, ytr)
        p = m.predict_proba(Xte)[:, 1]
        return {
            "k": len(cols),
            "auc": float(roc_auc_score(yte, p)),
            "ap": float(average_precision_score(yte, p)),
            "acc": float(accuracy_score(yte, (p >= 0.5).astype(int))),
        }

    board = {"po32": holdout(top_names), "f150": holdout(names)}
    (ROOT / "PO32_VS_F150.json").write_text(json.dumps(board, indent=2), encoding="utf-8")
    print(json.dumps(board, indent=2))
    print("top32 families", {f: sum(1 for x in top if x["family"] == f) for f in sorted({x["family"] for x in top})})


if __name__ == "__main__":
    main()
