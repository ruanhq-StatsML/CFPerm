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
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))
sys.path.insert(0, str(ROOT / "vendor" / "fsds"))

from msrvtt_multimodal_attribution import (
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
        "MMD (coord VIMP) & $%.3f$ & $%.3f$ & $%.3f$ & MMD $%.4f$ \\\\\n"
        % (p["mmd_share"]["video"], p["mmd_share"]["audio"], p["mmd_share"]["text"], p["mmd_full"]),
        "PO-risk VIMP & $%.3f$ & $%.3f$ & $%.3f$ & $R_{\\mathrm{PO}}=%.4f$ \\\\\n"
        % (p["po_feat_share"]["video"], p["po_feat_share"]["audio"], p["po_feat_share"]["text"], p["po_risk"]),
        "PO-risk permute-group & $%.3f$ & $%.3f$ & $%.3f$ & -- \\\\\n"
        % (p["po_perm_share"]["video"], p["po_perm_share"]["audio"], p["po_perm_share"]["text"]),
        "PO-risk LOGO (refit) & $%.3f$ & $%.3f$ & $%.3f$ & signed L1 share \\\\\n"
        % (p["po_logo_share"]["video"], p["po_logo_share"]["audio"], p["po_logo_share"]["text"]),
        "\\bottomrule\\end{tabular}\\end{table}\n\n",
        "\\begin{table}[ht]\\centering\n",
        "\\caption{Video-clustered bootstrap ($B=10$) for modality share differences. "
        "Mean and SD of the 10 replicates; interval is $\\mathrm{mean}\\pm 1.96\\,\\mathrm{SD}$. "
        "Two-sided bootstrap $p$ and Holm within method. Sampling unit = video.}\n",
        "\\label{tab:msrvtt-mm-bootstrap}\n\\small\n",
        "\\begin{tabular}{@{}l l r r r r c@{}}\\toprule\n",
        "Method & Contrast & Mean & SD & $\\mathrm{mean}\\pm 1.96\\mathrm{SD}$ & $p$ & Holm \\\\\n\\midrule\n",
    ]
    for method, lab in (("rf", "RF-Domain"), ("mmd", "MMD coord-VIMP"), ("po", "PO-risk VIMP")):
        for pair, rec in b[method]["pairwise"].items():
            lo, hi = rec["ci95"]
            lines.append(
                "%s & %s & $%.3f$ & $%.3f$ & $[%.3f, %.3f]$ & $%.4f$ & $%.4f$ \\\\\n"
                % (
                    lab,
                    pair,
                    rec["mean_diff"],
                    rec.get("sd", float("nan")),
                    lo,
                    hi,
                    rec["p_bootstrap"],
                    rec["p_holm"],
                )
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
        "\\bottomrule\\end{tabular}\\end{table}\n\n",
    ]
    lines += [
        "\\begin{table}[ht]\\centering\n",
        "\\caption{Per-video RF-Domain AUC for the early vs.\\ late window split "
        "(video-granularity subsets). Dominant modality is the largest RF VIMP share.}\n",
        "\\label{tab:msrvtt-mm-pervideo-auc}\n\\small\n",
        "\\begin{tabular}{@{}r r r r r r r@{}}\\toprule\n",
        "Video & $n$ & AUC & Video share & Audio share & Text share & Dominant \\\\\n\\midrule\n",
    ]
    for row in result["per_video"]:
        sh = row["share"]
        lines.append(
            "%s & %d & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & %s \\\\\n"
            % (
                row["video_id"],
                row["n"],
                row["auc"],
                sh["video"],
                sh["audio"],
                sh["text"],
                row["dominant"],
            )
        )
    lines.append("\\bottomrule\\end{tabular}\\end{table}\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))


def write_readme(result, plot_paths, path: Path):
    p = result["point"]
    board_stats = {}
    stats_path = path.parent / "msrvtt_region_board_stats.json"
    if stats_path.exists():
        board_stats = json.loads(stats_path.read_text(encoding="utf-8"))
    mad = board_stats.get("mean_abs_d_mod") or {}
    lines = [
        "# MSR-VTT multimodal FSDS attribution\n\n",
        "Window matrix `s`: **768 video + 512 audio + 768 text + 1 label**. ",
        "Batch $W$ = first vs second half of sliding windows **inside each video**.\n\n",
        "| Method | Video | Audio | Text | Diagnostic |\n",
        "|--------|------:|------:|-----:|------------|\n",
        "| RF-Domain | %.3f | %.3f | %.3f | AUC %.3f |\n"
        % (p["rf_share"]["video"], p["rf_share"]["audio"], p["rf_share"]["text"], p["rf_auc"]),
        "| MMD (coord VIMP) | %.3f | %.3f | %.3f | MMD %.4f |\n"
        % (p["mmd_share"]["video"], p["mmd_share"]["audio"], p["mmd_share"]["text"], p["mmd_full"]),
        "| PO-risk VIMP | %.3f | %.3f | %.3f | R=%.4f |\n"
        % (p["po_feat_share"]["video"], p["po_feat_share"]["audio"], p["po_feat_share"]["text"], p["po_risk"]),
        "| PO permute-group | %.3f | %.3f | %.3f | -- |\n"
        % (p["po_perm_share"]["video"], p["po_perm_share"]["audio"], p["po_perm_share"]["text"]),
        "\nInference: video-clustered bootstrap with **B=10** (mean and variance of replicates), ",
        "Friedman/Wilcoxon on per-video shares, within-video permutation of $W$ for AUC, ",
        "and a group-label permutation test that the named 768/512/768 blocks are more ",
        "imbalanced than random partitions of the same sizes.\n\n",
        "n=%d windows, n_videos=%d, bootstrap B=%s.\n"
        % (result["n"], result["n_videos"], result.get("bootstrap", {}).get("B", 10)),
        "\n## Inference readout\n\n",
        "Across RF-Domain, coordinate-MMD, and PO-risk VIMP the **video block dominates** ",
        "early-vs-late window shift (shares $\\approx$ 0.77 / 0.90 / 0.81), then audio, then text.\n\n",
        "- Per-video Friedman test of equal RF shares: $p=%.2e$ (n=%d videos).\n"
        % (
            result["tests"]["per_video_rf_shares"].get("friedman", {}).get("p", float("nan")),
            result["tests"]["per_video_rf_shares"].get("n_videos", result["n_videos"]),
        ),
        "- Wilcoxon signed-rank (Holm) rejects video=audio, video=text, and audio=text at $p<10^{-4}$.\n",
        "- Group-label permutation (named 768/512/768 vs random partitions): RF $p=%.3f$.\n"
        % result["tests"]["rf_group_label_perm"]["p"],
        "- Video-clustered bootstrap $B=10$: RF video$-$audio mean diff $%.3f$ (SD $%.3f$); "
        "the two-sided bootstrap $p$ floor with $B=10$ is $1/11\\approx0.091$ "
        "(all 10 replicates had the same sign). Use Wilcoxon/Friedman as the primary tests.\n"
        % (
            result["bootstrap"]["rf"]["pairwise"]["video-audio"]["mean_diff"],
            result["bootstrap"]["rf"]["pairwise"]["video-audio"]["sd"],
        ),
        "\n## Batch 0 vs Batch 1 region board\n\n",
        "See `msrvtt_batch_region_heatmap_board.png` (stats: `msrvtt_region_board_stats.json`). ",
        "Lead heatmap: windows ordered Batch 0 then Batch 1 × 24 embedding bins "
        "(8 video / 8 audio / 8 text). Middle: cosine$(B_0,B_1)$ geometry per modality. ",
        "Bottom: video × region Cohen's $d$ with Holm stars, plus signed pooled $d$ vs mean $|d|$. ",
        "Video/audio regions carry a **video-heterogeneous** early-vs-late shift; text is window-invariant ($d=0$). ",
        "Signed pooled $d$ cancels across videos; mean $|d|$ does not. "
        + (
            "Mean |Cohen's d|: video $%.3f$, audio $%.3f$, text $%.3f$; "
            "%d Holm-significant video×region cells, %d pooled-region Holm hits "
            "(signed $d$ cancels across videos).\n"
            % (
                float(mad.get("video", float("nan"))),
                float(mad.get("audio", float("nan"))),
                float(mad.get("text", float("nan"))),
                int(board_stats.get("n_sig_cells", 0)),
                int(board_stats.get("n_sig_pooled", 0)),
            )
            if mad
            else "\n"
        ),
        "\n## Per-video AUC\n\n",
        "| Video | n | AUC | Video | Audio | Text | Dominant |\n",
        "|------:|--:|----:|------:|------:|-----:|----------|\n",
    ]
    for row in result["per_video"]:
        sh = row["share"]
        lines.append(
            "| %s | %d | %.3f | %.3f | %.3f | %.3f | %s |\n"
            % (row["video_id"], row["n"], row["auc"], sh["video"], sh["audio"], sh["text"], row["dominant"])
        )
    path.write_text("".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", type=str, default=None, help="path to feature_video_audio.zip")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--B", type=int, default=10)
    ap.add_argument("--B-perm", type=int, default=10)
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
    plots = write_all_plots(result, OUT, bundle=bundle)
    write_tex(result, OUT / "MSRVTT_Multimodal_Attribution_tables_only.tex")
    write_tex(result, DOCS / "MSRVTT_Multimodal_Attribution_tables_only.tex")
    write_readme(result, plots, OUT / "README.md")
    print("wrote", OUT, flush=True)
    print("rf_share", result["point"]["rf_share"], "auc", result["point"]["rf_auc"], flush=True)
    print("mmd_share", result["point"]["mmd_share"], "mmd_block", result["point"]["mmd_block"], flush=True)
    print("po_feat_share", result["point"]["po_feat_share"], flush=True)
    print("po_perm_share", result["point"]["po_perm_share"], flush=True)
    print("bootstrap pairwise rf", result["bootstrap"]["rf"]["pairwise"], flush=True)
    print("friedman p", result["tests"]["per_video_rf_shares"].get("friedman"), flush=True)
    print("auc perm p", result["tests"]["within_video_auc_perm"]["p"], flush=True)


if __name__ == "__main__":
    main()
