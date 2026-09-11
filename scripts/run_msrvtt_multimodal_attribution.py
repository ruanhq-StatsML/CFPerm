#!/usr/bin/env python3
"""Run MSR-VTT multimodal FSDS attribution + resampling inference.

Loads window features from feature_video_audio.zip (array ``s`` with
768 video + 512 audio + 768 text + 1 label). W = early vs late windows
inside each video.

  python3 scripts/run_msrvtt_multimodal_attribution.py
  python3 scripts/run_msrvtt_multimodal_attribution.py --synthetic
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))
sys.path.insert(0, str(ROOT / "vendor" / "fsds"))

from msrvtt_multimodal_attribution import (  # noqa: E402
    find_feature_zip,
    load_window_bundle,
    make_synthetic_bundle,
    run_attribution,
    write_json,
    json_ready,
)
from msrvtt_attribution_plots import write_all_plots  # noqa: E402

OUT = ROOT / "results" / "msrvtt_multimodal_attribution"
DOCS = ROOT / "docs" / "method"


def write_tex(result, path: Path):
    p = result["point"]
    b = result["bootstrap"]
    t = result["tests"]
    lines = [
        "% MSR-VTT multimodal FSDS attribution · tables only\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\\centering\n",
        "\\caption{MSR-VTT sliding-window multimodal attribution. "
        "Columns are video/audio/text VIMP or LOGO shares. "
        "RF-Domain and MMD-LOCO target covariate shift $P(X)$; "
        "PO-risk LOGO targets concept drift through the pseudo-outcome. "
        "Batch $W$ is early vs.\\ late windows inside each video.}\n",
        "\\label{tab:msrvtt-mm-shares}\n\\small\n",
        "\\begin{tabular}{@{}l ccc c@{}}\\toprule\n",
        "Method & Video & Audio & Text & Diagnostic \\\\\n\\midrule\n",
        "RF-Domain VIMP & $%.3f$ & $%.3f$ & $%.3f$ & AUC $%.3f$ \\\\\n"
        % (p["rf_share"]["video"], p["rf_share"]["audio"], p["rf_share"]["text"], p["rf_auc"]),
        "MMD-LOCO & $%.3f$ & $%.3f$ & $%.3f$ & MMD $%.4f$ \\\\\n"
        % (p["mmd_share"]["video"], p["mmd_share"]["audio"], p["mmd_share"]["text"], p["mmd_full"]),
        "PO-risk LOGO & $%.3f$ & $%.3f$ & $%.3f$ & $R_{\\mathrm{PO}}=%.4f$ \\\\\n"
        % (p["po_logo_share"]["video"], p["po_logo_share"]["audio"], p["po_logo_share"]["text"], p["po_risk"]),
        "\\bottomrule\\end{tabular}\\end{table}\n\n",
        "\\begin{table}[ht]\\centering\n",
        "\\caption{Video-clustered bootstrap inference for modality share differences "
        "(percentile 95\\% CI, two-sided bootstrap $p$, Holm within method). "
        "Sampling unit = video.}\n",
        "\\label{tab:msrvtt-mm-bootstrap}\n\\small\n",
        "\\begin{tabular}{@{}l l r r r c@{}}\\toprule\n",
        "Method & Contrast & Mean & 95\\% CI & $p$ & Holm \\\\\n\\midrule\n",
    ]
    for method, lab in (("rf", "RF-Domain"), ("mmd", "MMD-LOCO"), ("po", "PO-risk")):
        for pair, rec in b[method]["pairwise"].items():
            lo, hi = rec["ci95"]
            lines.append(
                "%s & %s & $%.3f$ & $[%.3f, %.3f]$ & $%.4f$ & $%.4f$ \\\\\n"
                % (lab, pair, rec["mean_diff"], lo, hi, rec["p_bootstrap"], rec["p_holm"])
            )
    lines.append("\\bottomrule\\end{tabular}\\end{table}\n\n")

    pv = t["per_video_rf_shares"]
    lines += [
        "\\begin{table}[ht]\\centering\n",
        "\\caption{Per-video RF VIMP shares: Friedman test of equal modality contributions "
        "and pairwise Wilcoxon signed-rank tests (Holm-adjusted).}\n",
        "\\label{tab:msrvtt-mm-pervideo-tests}\n\\small\n",
        "\\begin{tabular}{@{}l r r@{}}\\toprule\n",
        "Test & Statistic & $p$ \\\\\n\\midrule\n",
        "Friedman (3 modalities) & $%.3f$ & $%.4f$ \\\\\n"
        % (pv.get("friedman", {}).get("stat", float("nan")), pv.get("friedman", {}).get("p", float("nan"))),
    ]
    for pair, rec in pv.get("wilcoxon", {}).items():
        lines.append(
            "Wilcoxon %s & $%.3f$ & $%.4f$ (Holm $%.4f$) \\\\\n"
            % (pair, rec.get("stat", float("nan")), rec.get("p", float("nan")), rec.get("p_holm", float("nan")))
        )
    perm = t["within_video_auc_perm"]
    gperm = t["rf_group_label_perm"]
    lines += [
        "Within-video AUC permutation & AUC $%.3f$ & $%.4f$ \\\\\n" % (perm["obs_auc"], perm["p"]),
        "RF group-label permutation (gap) & $%.3f$ & $%.4f$ \\\\\n" % (gperm["obs_gap"], gperm["p"]),
        "\\bottomrule\\end{tabular}\\end{table}\n",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))


def write_readme(result, plot_paths, path: Path):
    p = result["point"]
    lines = [
        "# MSR-VTT multimodal FSDS attribution\n\n",
        "Window matrix `s`: **768 video + 512 audio + 768 text + 1 label**. ",
        "Batch $W$ = first vs second half of sliding windows **inside each video**.\n\n",
        "| Method | Video | Audio | Text | Diagnostic |\n",
        "|--------|------:|------:|-----:|------------|\n",
        "| RF-Domain | %.3f | %.3f | %.3f | AUC %.3f |\n"
        % (p["rf_share"]["video"], p["rf_share"]["audio"], p["rf_share"]["text"], p["rf_auc"]),
        "| MMD-LOCO | %.3f | %.3f | %.3f | MMD %.4f |\n"
        % (p["mmd_share"]["video"], p["mmd_share"]["audio"], p["mmd_share"]["text"], p["mmd_full"]),
        "| PO-risk LOGO | %.3f | %.3f | %.3f | R=%.4f |\n"
        % (p["po_logo_share"]["video"], p["po_logo_share"]["audio"], p["po_logo_share"]["text"], p["po_risk"]),
        "\nInference: video-clustered bootstrap, Holm pairwise tests, ",
        "Friedman/Wilcoxon on per-video shares, within-video permutation of $W$ for AUC, ",
        "and a group-label permutation test that the named 768/512/768 blocks are more ",
        "imbalanced than random partitions of the same sizes.\n\n",
        "n=%d windows, n_videos=%d.\n" % (result["n"], result["n_videos"]),
    ]
    path.write_text("".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", type=str, default=None, help="path to feature_video_audio.zip")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--B", type=int, default=79)
    ap.add_argument("--B-perm", type=int, default=79)
    ap.add_argument("--n-estimators", type=int, default=120)
    ap.add_argument("--bootstrap-estimators", type=int, default=50)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    if args.synthetic:
        bundle = make_synthetic_bundle(seed=args.seed)
        print("synthetic bundle n=%d videos=%d" % (len(bundle.y), len(set(bundle.video_id))), flush=True)
    else:
        zip_path = Path(args.zip) if args.zip else find_feature_zip(ROOT)
        if zip_path is None:
            # last-chance glob including parent dirs
            hits = list(Path("/workspace").rglob("feature_video_audio.zip"))
            zip_path = hits[0] if hits else None
        if zip_path is None:
            raise SystemExit(
                "feature_video_audio.zip not found. Pass --zip PATH or place the archive "
                "at data/msrvtt/feature_video_audio.zip"
            )
        print("loading", zip_path, flush=True)
        bundle = load_window_bundle(zip_path, root=ROOT)
        print(
            "loaded s-layout X=%s n_videos=%d n0=%d n1=%d"
            % (bundle.X.shape, len(set(map(int, bundle.video_id))), (bundle.W == 0).sum(), (bundle.W == 1).sum()),
            flush=True,
        )

    result = run_attribution(
        bundle,
        B=args.B,
        B_perm=args.B_perm,
        seed=args.seed,
        n_estimators=args.n_estimators,
        bootstrap_estimators=args.bootstrap_estimators,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "msrvtt_multimodal_attribution.json", {k: v for k, v in result.items() if k != "arrays"})
    np_dir = OUT / "arrays"
    np_dir.mkdir(exist_ok=True)
    import numpy as np

    for name, arr in result["arrays"].items():
        np.save(np_dir / ("%s.npy" % name), arr)
    plots = write_all_plots(result, OUT)
    write_tex(result, OUT / "MSRVTT_Multimodal_Attribution_tables_only.tex")
    write_tex(result, DOCS / "MSRVTT_Multimodal_Attribution_tables_only.tex")
    write_readme(result, plots, OUT / "README.md")
    print("wrote", OUT, flush=True)
    print("rf_share", result["point"]["rf_share"], "auc", result["point"]["rf_auc"], flush=True)
    print("mmd_share", result["point"]["mmd_share"], flush=True)
    print("po_logo_share", result["point"]["po_logo_share"], flush=True)
    print("bootstrap pairwise rf", result["bootstrap"]["rf"]["pairwise"], flush=True)
    print("auc perm p", result["tests"]["within_video_auc_perm"]["p"], flush=True)


if __name__ == "__main__":
    main()
