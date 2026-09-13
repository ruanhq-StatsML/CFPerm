#!/usr/bin/env python3
"""PO-risk tail → DGA alignments → mean-diff subgroup on a planted stream.

  python3 scripts/run_dga_po_localize.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from dga_po_localize import (  # noqa: E402
    localize_one_hop,
    make_planted_subgroup_stream,
    recovery_against_plant,
    run_dga_po_localize,
    top_feature_hit,
)
from msrvtt_multimodal_attribution import write_json  # noqa: E402

OUT = ROOT / "results" / "dga_po_localize"
DOCS = ROOT / "docs" / "method"


def plot_mean_diff(ranked, shift_coords, path):
    import matplotlib.pyplot as plt

    ranked = list(ranked)
    names = [r["name"] for r in ranked]
    d = np.array([r["cohens_d"] for r in ranked], dtype=float)
    planted = {int(j) for j in shift_coords}
    colors = ["#C45C26" if int(r["j"]) in planted else "#5B3A8C" for r in ranked]
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    ypos = np.arange(len(names))[::-1]
    ax.barh(ypos, d, color=colors, height=0.72)
    ax.axvline(0.0, color="#888", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=8.5)
    ax.set_xlabel("Cohen's d (PO tail − complement in $B_{t-1}$)")
    ax.set_title("Post-hoc subgroup: mean-diff localization", loc="left", fontsize=12)
    ax.grid(True, axis="x", color="#eee")
    fig.text(
        0.04,
        0.02,
        "Orange = planted shift coordinates. DGA aligns to the PO tail, not the whole last batch.",
        fontsize=8.0,
        color="#666",
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def write_tex(rec, hop, hit, recov, path):
    top = hop["ranked_features"][:5]
    lines = [
        r"% DGA + PO-risk subset localization. Auto-generated planted stream.",
        r"\begin{table}[ht]\centering",
        r"\caption{PO-risk tail localizes $D_{\mathrm{spe}}=B_{t-1}$, then DGA aligns to that",
        r"pocket. Post-hoc subgroup analysis is the mean-diff / Cohen's $d$ ranking inside",
        r"$B_{t-1}$ (high-PO vs complement), discretized at midpoints. Planted synthetic",
        r"stream; not an Amazon MSE claim.}",
        r"\label{tab:dga-po-subset}",
        r"\small",
        r"\begin{tabular}{@{}lcc@{}}\toprule",
        r"Feature & Cohen's $d$ & planted \\",
        r"\midrule",
    ]
    planted = set(hop.get("shift_coords") or rec["meta"].get("shift_coords") or [])
    for row in top:
        mark = r"yes" if int(row["j"]) in planted else r"no"
        lines.append(r"%s & $%.2f$ & %s \\" % (row["name"], row["cohens_d"], mark))
    lines.extend(
        [
            r"\midrule",
            r"PO-tail precision vs plant & \multicolumn{2}{c}{$%.2f$} \\" % recov["precision"],
            r"PO-tail recall vs plant & \multicolumn{2}{c}{$%.2f$} \\" % recov["recall"],
            r"top-$k$ planted-coord recall & \multicolumn{2}{c}{$%.2f$} \\" % hit["recall"],
            r"online MSE (\texttt{dga\_po\_ridge}) & \multicolumn{2}{c}{$%.3f$} \\"
            % rec["online_mse"],
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def main():
    stream = make_planted_subgroup_stream(n_batches=5, n_per=160, p=12, seed=2026)
    rec = run_dga_po_localize(stream, po_q=0.35, k_features=3)
    spe = int(stream.batch.max()) - 1
    hop = localize_one_hop(stream.X, stream.y, stream.batch, spe, q=0.35, k_features=3)
    hop["shift_coords"] = stream.meta["shift_coords"]
    recov = recovery_against_plant(stream, hop["tail_mask"], spe)
    hit = top_feature_hit(hop["ranked_features"], stream.meta["shift_coords"], k=4)
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "online_mse": rec["online_mse"],
        "online_path": rec["online_path"],
        "last_top_features": rec["last_top_features"],
        "last_rules": rec["last_rules"],
        "recovery": recov,
        "feature_hit": hit,
        "n_spe": hop["n_spe"],
        "n_tail": hop["n_tail"],
        "shift_coords": stream.meta["shift_coords"],
        "planted_frac": stream.meta["planted_frac"],
        "history": rec["history"],
    }
    write_json(OUT / "dga_po_localize.json", payload)
    plot_mean_diff(hop["ranked_features"], stream.meta["shift_coords"], OUT / "dga_po_localize.png")
    write_tex(rec, hop, hit, recov, OUT / "DGA_po_subset.tex")
    write_tex(rec, hop, hit, recov, DOCS / "DGA_po_subset.tex")
    readme = OUT / "README.md"
    readme.write_text(
        "# DGA + PO-risk subset localization\n\n"
        "Planted synthetic stream. PO-risk tail localizes $D_{spe}$, DGA aligns to "
        "that pocket, mean-diff discretizes the subgroup.\n\n"
        "| metric | value |\n|---|---:|\n"
        "| online MSE | %.3f |\n" % rec["online_mse"]
        + "| PO-tail precision | %.3f |\n" % recov["precision"]
        + "| PO-tail recall | %.3f |\n" % recov["recall"]
        + "| planted-coord recall (top-4) | %.3f |\n" % hit["recall"]
        + "| n_spe / n_tail | %d / %d |\n" % (hop["n_spe"], hop["n_tail"])
        + "\nOrange bars in `dga_po_localize.png` are planted coordinates.\n"
    )
    print(json.dumps({k: payload[k] for k in ("online_mse", "recovery", "feature_hit", "n_spe", "n_tail", "last_rules")}, indent=2))


if __name__ == "__main__":
    main()
