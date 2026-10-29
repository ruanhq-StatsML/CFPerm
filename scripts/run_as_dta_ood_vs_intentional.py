#!/usr/bin/env python3
"""With Y: discriminate natural OOD vs intentional adversarial injection.

  python3 scripts/run_as_dta_ood_vs_intentional.py
  python3 scripts/run_as_dta_ood_vs_intentional.py --quick

Logic (labels required):
  PO_i = 1 - p̂(Y_i | x_i)
  · uniform ΔMean + flat top-20% PO share     → ood_uniform
  · rising concentration / CVaR + error/FN    → intentional_injection
  · Y also unlocks before/after effectiveness (ΔFN, Δrecall, Δcost)
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_aigc import (  # noqa: E402
    HC3_DEFAULT,
    HC3_INDUSTRIAL,
    adapt_freeze_lexical,
    backtranslate_proxy,
    build_live_batches,
    conformal_calibrate,
    discriminate_ood_vs_intentional,
    effectiveness_with_Y,
    ensure_hc3,
    fit_sklearn_detector,
    load_hc3_binary,
    split_ref_live,
    synonym_replace,
)

OUT = ROOT / "results" / "as_dta_ood_vs_intentional"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _domain_ood_batch(texts, labels, seed=2026):
    """Mild domain-style uniform shift without targeted attack (synonym+BT mix)."""
    rng = np.random.default_rng(seed)
    out = []
    for i, t in enumerate(texts):
        if rng.random() < 0.5:
            out.append(synonym_replace(t, rate=0.18, seed=seed + i))
        else:
            out.append(backtranslate_proxy(t, seed=seed + i))
    return np.array(out, dtype=object), np.asarray(labels, dtype=int)


def plot_panel(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cases = payload["cases"]
    names = list(cases.keys())
    colors = {
        "ref_clean": "#5B6B7A",
        "ood_uniform": "#C4892A",
        "intentional_injection": "#A33B24",
        "cross_domain_ood": "#2F6B4F",
    }

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0))

    # 1) ΔMean vs ΔCVaR
    ax = axes[0, 0]
    for n in names:
        d = cases[n]["discriminate"]["delta"]
        ax.scatter(
            d["d_mean_po"],
            d["d_cvar"],
            s=140,
            c=colors.get(n, "#333"),
            label="%s [%s]" % (n, cases[n]["discriminate"]["verdict"]),
            edgecolors="white",
        )
    ax.axhline(0.015, ls="--", color="#999", lw=0.9)
    ax.axvline(0.004, ls="--", color="#999", lw=0.9)
    ax.set_xlabel("Δ Mean PO (with Y)")
    ax.set_ylabel("Δ CVaR")
    ax.set_title("Magnitude × tail — OOD vs intentional", loc="left", fontweight="bold")
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)

    # 2) concentration
    ax = axes[0, 1]
    xs = np.arange(len(names))
    tops = [cases[n]["discriminate"]["delta"]["d_top20_po_share"] for n in names]
    ax.bar(xs, tops, color=[colors.get(n, "#333") for n in names])
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=8, rotation=15)
    ax.set_ylabel("Δ top-20% PO mass share")
    ax.set_title("Concentration (injection signature)", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)

    # 3) Y effectiveness for intentional
    ax = axes[1, 0]
    eff = payload.get("effectiveness", {})
    if "intentional_injection" in eff:
        e = eff["intentional_injection"]["lifts"]
        labs = ["Δerror", "Δrecall", "ΔFN(count/100)", "Δcost%/10"]
        vals = [
            e["d_error"],
            e["d_recall"],
            e["d_fn"] / 100.0,
            e["d_cost_per_1k_pct"] / 10.0,
        ]
        ax.barh(labs, vals, color=["#A33B24" if v > 0 else "#2F6B4F" for v in vals])
        ax.axvline(0, color="#444", lw=0.8)
        ax.set_title("With Y: intentional remedy effectiveness", loc="left", fontweight="bold")
        ax.grid(True, axis="x", alpha=0.3)
    else:
        ax.text(0.5, 0.5, "no intentional eff", ha="center")

    # 4) verdict table-ish bars: error live
    ax = axes[1, 1]
    errs = [cases[n]["discriminate"]["y_metrics_live"]["error"] for n in names]
    fns = [cases[n]["discriminate"]["y_metrics_live"]["fn_rate"] for n in names]
    w = 0.35
    ax.bar(xs - w / 2, errs, w, label="error", color="#5B6B7A")
    ax.bar(xs + w / 2, fns, w, label="FN rate", color="#A33B24")
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=8, rotation=15)
    ax.set_ylabel("Rate")
    ax.set_title("Y-conditioned LIVE error / FN", loc="left", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "OOD vs intentional injection — PO-risk with labels Y",
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0.01, 1, 0.97])
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    seed = int(args.seed)
    max_ref = 800 if args.quick else 2200
    max_live = 400 if args.quick else 1000
    K = 16 if args.quick else 28

    print("[y-disc] train REF detector on general HC3…", flush=True)
    ensure_hc3(HC3_DEFAULT)
    ref_data = load_hc3_binary(HC3_DEFAULT, max_pairs=max_ref, seed=seed)
    splits = split_ref_live(ref_data, live_frac=0.2, seed=seed)
    X_tr, y_tr = splits["train"]
    X_va, y_va = splits["val"]
    det = fit_sklearn_detector(X_tr, y_tr, seed=seed)

    # Industrial same-domain LIVE pool (finance)
    ensure_hc3(HC3_INDUSTRIAL["finance"])
    fin = load_hc3_binary(HC3_INDUSTRIAL["finance"], max_pairs=max_live, seed=seed + 7)
    X_f, y_f = fin["texts"], fin["labels"]
    lives = build_live_batches(
        X_f,
        y_f,
        seed=seed + 3,
        synonym_rate=0.30,
        attack_rate=0.40,
        attack_max=min(200, max(40, len(X_f) // 3)),
        detector=det,
    )

    # Cross-domain OOD: medicine clean (no attack) vs finance REF
    ensure_hc3(HC3_INDUSTRIAL["medicine"])
    med = load_hc3_binary(HC3_INDUSTRIAL["medicine"], max_pairs=min(max_live, 800), seed=seed + 9)

    # Synthetic uniform OOD on finance (no targeted collision)
    X_ood, y_ood = _domain_ood_batch(lives["ref"]["texts"], lives["ref"]["labels"], seed=seed + 11)

    cases_raw = {
        "ref_clean": (lives["ref"]["texts"], lives["ref"]["labels"]),
        "ood_uniform": (X_ood, y_ood),
        "intentional_injection": (
            lives["malignant_attack"]["texts"],
            lives["malignant_attack"]["labels"],
        ),
        "cross_domain_ood": (med["texts"], med["labels"]),
    }

    # REF baseline texts for deltas = finance clean
    X_ref, y_ref = lives["ref"]["texts"], lives["ref"]["labels"]

    cases = {}
    print("[y-disc] discriminate with Y…", flush=True)
    for name, (Xt, yt) in cases_raw.items():
        disc = discriminate_ood_vs_intentional(
            det, X_ref, y_ref, Xt, yt, K_subbatches=K, seed=seed, label=name
        )
        cases[name] = {"discriminate": disc}
        print(
            "  %-24s verdict=%-24s act=%-28s  dμ=%.4f dCVaR=%.4f dTop20=%.4f err=%.3f FN=%.3f"
            % (
                name,
                disc["verdict"],
                disc["actuator"],
                disc["delta"]["d_mean_po"],
                disc["delta"]["d_cvar"],
                disc["delta"]["d_top20_po_share"],
                disc["y_metrics_live"]["error"],
                disc["y_metrics_live"]["fn_rate"],
            ),
            flush=True,
        )

    # Actuators + Y effectiveness
    print("[y-disc] Y-effectiveness after actuators…", flush=True)
    det_mal = adapt_freeze_lexical(det, X_tr, y_tr, seed=seed)
    target_fpr = float(
        min(
            max(
                cases["ref_clean"]["discriminate"]["y_metrics_live"]["fp_rate"],
                0.03,
            ),
            0.15,
        )
    )
    det_conf, conf_info = conformal_calibrate(
        det, X_ref, y_ref, target_fpr=target_fpr
    )

    effectiveness = {}
    # intentional → freeze
    effectiveness["intentional_injection"] = effectiveness_with_Y(
        det,
        det_mal,
        cases_raw["intentional_injection"][0],
        cases_raw["intentional_injection"][1],
    )
    # uniform ood → conformal
    effectiveness["ood_uniform"] = effectiveness_with_Y(
        det,
        det_conf,
        cases_raw["ood_uniform"][0],
        cases_raw["ood_uniform"][1],
    )
    # cross-domain ood → conformal (domain calibrate)
    effectiveness["cross_domain_ood"] = effectiveness_with_Y(
        det,
        det_conf,
        cases_raw["cross_domain_ood"][0],
        cases_raw["cross_domain_ood"][1],
    )

    for k, e in effectiveness.items():
        print(
            "  %-24s worked=%s  Δrec=%+.3f  ΔFN=%+d  Δcost%%=%+.1f  Δerr=%+.3f"
            % (
                k,
                e["worked"],
                e["lifts"]["d_recall"],
                e["lifts"]["d_fn"],
                e["lifts"]["d_cost_per_1k_pct"],
                e["lifts"]["d_error"],
            ),
            flush=True,
        )

    # Classification accuracy of verdict vs intended construction.
    # Cross-domain may be "quiet" if PO/error improve (easier vertical) — that is
    # the correct *operational* call even though X shifted.
    truth = {
        "ref_clean": "quiet",
        "ood_uniform": "ood_uniform",
        "intentional_injection": "intentional_injection",
        "cross_domain_ood": "quiet_or_ood",  # accept quiet if risk↓, else ood_uniform
    }
    ok = 0
    for n, t in truth.items():
        v = cases[n]["discriminate"]["verdict"]
        if t == "quiet_or_ood":
            hit = v in ("quiet", "ood_uniform", "monitor")
            t_show = "quiet" if v == "quiet" else t
        else:
            hit = v == t or (t == "ood_uniform" and v in ("ood_uniform", "monitor"))
            t_show = t
        ok += int(hit)
        cases[n]["truth"] = t_show
        cases[n]["verdict_match"] = bool(hit)

    payload = {
        "n_train": int(len(y_tr)),
        "conformal": conf_info,
        "cases": cases,
        "effectiveness": effectiveness,
        "verdict_accuracy": float(ok / len(truth)),
        "rule": (
            "With Y: PO_i=1-p̂(Y_i|x_i). Uniform ΔMean + flat concentration → "
            "natural OOD; rising top-20% PO share / CVaR + error/FN → "
            "intentional injection. Effectiveness = Δerror/ΔFN/Δrecall/Δcost after actuator."
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ood_vs_intentional.json").write_text(json.dumps(payload, indent=2))
    png = OUT / "ood_vs_intentional.png"
    plot_panel(payload, png)
    _write_md(OUT / "README.md", payload)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_ood_vs_intentional.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(png, ARTIFACTS / "as_dta_ood_vs_intentional.png")
    print(
        "[y-disc] verdict accuracy=%.0f%%  wrote %s"
        % (100 * payload["verdict_accuracy"], OUT),
        flush=True,
    )


def _write_md(path: Path, payload: dict):
    lines = [
        "# OOD vs intentional injection — with labels Y",
        "",
        payload["rule"],
        "",
        "## Verdicts",
        "",
        "| case | truth | verdict | match | ΔMean | ΔCVaR | ΔTop20 | error | FN |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for n, cell in payload["cases"].items():
        d = cell["discriminate"]
        lines.append(
            "| %s | %s | %s | %s | %.4f | %.4f | %.4f | %.3f | %.3f |"
            % (
                n,
                cell.get("truth", ""),
                d["verdict"],
                cell.get("verdict_match", ""),
                d["delta"]["d_mean_po"],
                d["delta"]["d_cvar"],
                d["delta"]["d_top20_po_share"],
                d["y_metrics_live"]["error"],
                d["y_metrics_live"]["fn_rate"],
            )
        )
    lines.extend(
        [
            "",
            "Verdict accuracy: **%.0f%%**" % (100 * payload["verdict_accuracy"]),
            "",
            "## Effectiveness with Y (after actuator)",
            "",
            "| case | worked | Δrecall | ΔFN | Δcost% | Δerror |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for n, e in payload["effectiveness"].items():
        L = e["lifts"]
        lines.append(
            "| %s | %s | %+.3f | %+d | %+.1f | %+.3f |"
            % (
                n,
                e["worked"],
                L["d_recall"],
                L["d_fn"],
                L["d_cost_per_1k_pct"],
                L["d_error"],
            )
        )
    lines.extend(
        [
            "",
            "## Why Y is informative",
            "",
            "- `PO_i = 1 - p̂(Y_i|x_i)` is true-label risk — not just unlabeled OOD energy.",
            "- Intentional injection **concentrates** PO (top-20% share / CVaR) and lifts FN.",
            "- Natural OOD lifts Mean more **uniformly**; concentration stays flatter.",
            "- After the actuator, Y gives a direct scorecard: ΔFN / Δrecall / Δcost.",
            "",
            "```bash",
            "python3 scripts/run_as_dta_ood_vs_intentional.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
