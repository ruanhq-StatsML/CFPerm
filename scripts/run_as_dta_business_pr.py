#!/usr/bin/env python3
"""AS-DTA business metrics: industrial HC3 verticals + precision-recall + ROI.

Train on general HC3 REF, evaluate LIVE shifts on **finance** / **medicine**
(industrial integrity verticals). Reports:

  - drift-shape detection (quiet / benign / malignant) from Mean×Peak
  - precision / recall / FPR before vs after targeted adaptation
  - expected cost = c_FP·FP + c_FN·FN  (default $50 / $20)

  python3 scripts/run_as_dta_business_pr.py
  python3 scripts/run_as_dta_business_pr.py --quick
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
    COST_FN_DEFAULT,
    COST_FP_DEFAULT,
    HC3_DEFAULT,
    HC3_INDUSTRIAL,
    adapt_freeze_lexical,
    build_live_batches,
    business_pr_metrics,
    conformal_calibrate,
    diagnose_batch,
    ensure_hc3,
    fit_sklearn_detector,
    load_hc3_binary,
    route_action,
    shape_detection_score,
    split_ref_live,
)

OUT = ROOT / "results" / "as_dta_business_pr"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _pr_pair(det, texts, labels, cost_fp, cost_fn):
    return business_pr_metrics(det, texts, labels, cost_fp=cost_fp, cost_fn=cost_fn)


def plot_business_board(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    verticals = payload["verticals"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))

    # --- PR curves malignant before/after (finance) ---
    ax = axes[0, 0]
    for vert, style in (("finance", "-"), ("medicine", "--")):
        cell = verticals[vert]["malignant_attack"]
        b = cell["before"]["pr_curve"]
        a = cell["after"]["pr_curve"]
        ax.plot(b["recall"], b["precision"], style, color="#A33B24", lw=1.6, label="%s before" % vert)
        ax.plot(a["recall"], a["precision"], style, color="#2F6B4F", lw=1.6, label="%s after" % vert)
    ax.set_xlabel("Recall (catch AI)")
    ax.set_ylabel("Precision (flagged ∩ AI)")
    ax.set_title("PR curve — malignant LIVE (before/after freeze+aug)", loc="left", fontweight="bold")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0.5, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.5, ncol=2)

    # --- Recall / Precision bars before→after ---
    ax = axes[0, 1]
    rows = []
    labels = []
    for vert in verticals:
        for kind in ("benign_synonym", "malignant_attack"):
            cell = verticals[vert][kind]
            labels.append("%s\n%s" % (vert[:3], kind.split("_")[-1][:6]))
            rows.append(
                (
                    cell["before"]["recall"],
                    cell["after"]["recall"],
                    cell["before"]["precision"],
                    cell["after"]["precision"],
                )
            )
    x = np.arange(len(rows))
    w = 0.2
    ax.bar(x - 1.5 * w, [r[0] for r in rows], w, label="rec before", color="#C4892A")
    ax.bar(x - 0.5 * w, [r[1] for r in rows], w, label="rec after", color="#2F6B4F")
    ax.bar(x + 0.5 * w, [r[2] for r in rows], w, label="prec before", color="#A33B24", alpha=0.7)
    ax.bar(x + 1.5 * w, [r[3] for r in rows], w, label="prec after", color="#1F6F8B", alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylim(0.5, 1.05)
    ax.set_title("Precision & Recall — adaptation effect", loc="left", fontweight="bold")
    ax.legend(fontsize=7, ncol=2, loc="lower left")
    ax.grid(True, axis="y", alpha=0.3)

    # --- Cost per 1k ---
    ax = axes[1, 0]
    cost_labels, cost_b, cost_a = [], [], []
    for vert in verticals:
        for kind in ("benign_synonym", "malignant_attack"):
            cell = verticals[vert][kind]
            cost_labels.append("%s/%s" % (vert[:3], kind.split("_")[-1][:6]))
            cost_b.append(cell["before"]["cost_per_1k_usd"])
            cost_a.append(cell["after"]["cost_per_1k_usd"])
    x = np.arange(len(cost_labels))
    ax.bar(x - 0.18, cost_b, 0.36, label="before", color="#A33B24")
    ax.bar(x + 0.18, cost_a, 0.36, label="after", color="#2F6B4F")
    ax.set_xticks(x)
    ax.set_xticklabels(cost_labels, fontsize=8)
    ax.set_ylabel("Expected cost / 1k items (USD)")
    ax.set_title(
        "Business cost  ($%g·FP + $%g·FN)"
        % (payload["cost_fp"], payload["cost_fn"]),
        loc="left",
        fontweight="bold",
    )
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    # --- Shape detection ---
    ax = axes[1, 1]
    shape = payload["shape_detection"]
    mats = shape["confusion"]
    keys = list(mats.keys())
    mat = np.array([[mats[t].get(p, 0) for p in keys] for t in keys], dtype=float)
    im = ax.imshow(mat, cmap="YlOrBr")
    ax.set_xticks(range(len(keys)))
    ax.set_yticks(range(len(keys)))
    ax.set_xticklabels(keys, fontsize=8)
    ax.set_yticklabels(keys, fontsize=8)
    ax.set_xlabel("Predicted shape (Mean×Peak route)")
    ax.set_ylabel("True LIVE kind")
    ax.set_title(
        "Drift-shape detection  acc=%.0f%%" % (100 * shape["accuracy"]),
        loc="left",
        fontweight="bold",
    )
    for i in range(len(keys)):
        for j in range(len(keys)):
            ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046)

    fig.suptitle(
        "AS-DTA business board — HC3 finance/medicine industrial LIVE",
        fontsize=12.5,
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0.01, 1, 0.97])
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def run_business(args):
    cost_fp = float(args.cost_fp)
    cost_fn = float(args.cost_fn)
    max_ref = 1500 if args.quick else int(args.max_ref)
    max_live = 600 if args.quick else int(args.max_live)
    K = 16 if args.quick else int(args.K)
    seed = int(args.seed)

    print("[biz] fine-tune REF on general HC3…", flush=True)
    ensure_hc3(HC3_DEFAULT)
    ref_data = load_hc3_binary(HC3_DEFAULT, max_pairs=max_ref, seed=seed)
    # use most of general as train/val; keep a small general live only for shape baseline
    splits = split_ref_live(ref_data, live_frac=0.15, seed=seed)
    X_tr, y_tr = splits["train"]
    X_va, y_va = splits["val"]
    det = fit_sklearn_detector(X_tr, y_tr, seed=seed)
    print("  train n=%d  val n=%d" % (len(X_tr), len(X_va)), flush=True)

    det_mal = adapt_freeze_lexical(det, X_tr, y_tr, seed=seed)

    verticals = {}
    shape_true, shape_pred = [], []
    true_map = {
        "ref": "quiet",
        "benign_backtranslate": "quiet",  # expect near-hold
        "benign_synonym": "benign",
        "malignant_attack": "malignant",
    }
    conf_info = {"note": "per-vertical conformal on industrial REF humans"}

    for vert in ("finance", "medicine"):
        path = HC3_INDUSTRIAL[vert]
        ensure_hc3(path)
        print("[biz] LIVE vertical=%s" % vert, flush=True)
        live_data = load_hc3_binary(path, max_pairs=max_live, seed=seed + 11)
        X_live, y_live = live_data["texts"], live_data["labels"]
        lives = build_live_batches(
            X_live,
            y_live,
            seed=seed + 3,
            synonym_rate=0.30,
            attack_rate=0.40,
            attack_max=min(220, max(40, len(X_live) // 3)),
            detector=det,
        )
        # Domain-calibrated conformal: match FPR of this vertical's clean REF
        # (absolute 5% on general-val over-flags under industrial covariate shift).
        ref_pr = _pr_pair(det, lives["ref"]["texts"], lives["ref"]["labels"], cost_fp, cost_fn)
        target_fpr_vert = float(min(max(ref_pr["fpr"], 0.03), 0.15))
        det_conf, conf_v = conformal_calibrate(
            det, lives["ref"]["texts"], lives["ref"]["labels"], target_fpr=target_fpr_vert
        )
        conf_info[vert] = {**conf_v, "target_fpr_matched_to_ref": target_fpr_vert}

        # REF diagnosis on clean industrial LIVE for relative routing
        ref_diag = diagnose_batch(
            lives["ref"]["texts"], lives["ref"]["labels"], det, K_subbatches=K, seed=seed + 7, kind="ref"
        )
        vert_out = {"ref_diag": ref_diag.to_dict(), "batches": {}}
        for kind, batch in lives.items():
            diag = diagnose_batch(
                batch["texts"], batch["labels"], det, K_subbatches=K, seed=seed + 7, kind=kind
            )
            shape = shape_detection_score(diag, ref_diag)
            shape_true.append(true_map[kind])
            # collapse predicted quiet/monitor using continuous malignant_score
            pred = shape["shape"]
            if pred == "monitor":
                pred = "malignant" if shape["malignant_score"] >= 0.08 else (
                    "benign" if shape["d_mean"] > 0.004 else "quiet"
                )
            # hard mean shock override (industrial mixture attacks)
            if shape["d_mean"] >= 0.04 or shape["d_cvar"] >= 0.02:
                pred = "malignant"
                shape["shape"] = "malignant"
            shape_pred.append(pred)

            before = _pr_pair(det, batch["texts"], batch["labels"], cost_fp, cost_fn)

            def _apply(actuator_name):
                if actuator_name == "freeze_lexical_adv_aug":
                    return _pr_pair(det_mal, batch["texts"], batch["labels"], cost_fp, cost_fn)
                if actuator_name == "conformal_domain_calibrated":
                    return _pr_pair(det_conf, batch["texts"], batch["labels"], cost_fp, cost_fn)
                return before

            # Oracle medicine = by true LIVE kind (upper bound on actuator value)
            if kind == "malignant_attack":
                oracle_act = "freeze_lexical_adv_aug"
            elif kind == "benign_synonym":
                oracle_act = "conformal_domain_calibrated"
            else:
                oracle_act = "hold"
            # Routed medicine = by detected shape (end-to-end AS-DTA)
            if pred == "malignant":
                routed_act = "freeze_lexical_adv_aug"
            elif pred == "benign":
                routed_act = "conformal_domain_calibrated"
            else:
                routed_act = "hold"

            after_oracle = _apply(oracle_act)
            after_routed = _apply(routed_act)
            # Primary "after" for PR board = oracle (answers "does the right medicine work?")
            after = after_oracle
            actuator = oracle_act

            vert_out[kind] = {
                "diagnosis": diag.to_dict(),
                "shape": shape,
                "before": before,
                "after": after,
                "after_routed": after_routed,
                "actuator": actuator,
                "actuator_routed": routed_act,
                "route": route_action(diag, ref=ref_diag),
                "detected_shape": pred,
            }
            print(
                "  %-22s detect=%-10s oracle=%-28s  rec=%.3f→%.3f  prec=%.3f→%.3f  cost/1k=$%.0f→$%.0f  (routed cost→$%.0f)"
                % (
                    kind,
                    pred,
                    oracle_act,
                    before["recall"],
                    after["recall"],
                    before["precision"],
                    after["precision"],
                    before["cost_per_1k_usd"],
                    after["cost_per_1k_usd"],
                    after_routed["cost_per_1k_usd"],
                ),
                flush=True,
            )
        verticals[vert] = vert_out

    # shape detection confusion
    labels_u = ["quiet", "benign", "malignant"]
    conf = {t: {p: 0 for p in labels_u} for t in labels_u}
    for t, p in zip(shape_true, shape_pred):
        if t not in conf:
            continue
        if p not in conf[t]:
            p = "benign"
        conf[t][p] += 1
    acc = float(np.mean([t == p for t, p in zip(shape_true, shape_pred)]))

    # ROI rollup: oracle medicine vs static; also routed end-to-end
    cost_saved = 0.0
    cost_saved_routed = 0.0
    recall_lifts = []
    for vert, vo in verticals.items():
        for kind in ("benign_synonym", "malignant_attack"):
            b = vo[kind]["before"]
            a = vo[kind]["after"]
            ar = vo[kind]["after_routed"]
            cost_saved += b["expected_cost_usd"] - a["expected_cost_usd"]
            cost_saved_routed += b["expected_cost_usd"] - ar["expected_cost_usd"]
            if kind == "malignant_attack":
                recall_lifts.append(a["recall"] - b["recall"])

    payload = {
        "cost_fp": cost_fp,
        "cost_fn": cost_fn,
        "conformal": conf_info,
        "ref_train_n": int(len(X_tr)),
        "verticals": {
            vert: {
                k: {
                    "route": verticals[vert][k]["route"],
                    "actuator": verticals[vert][k]["actuator"],
                    "shape": verticals[vert][k]["shape"],
                    "detected_shape": verticals[vert][k]["detected_shape"],
                    "actuator_routed": verticals[vert][k]["actuator_routed"],
                    "before": {
                        kk: verticals[vert][k]["before"][kk]
                        for kk in (
                            "precision",
                            "recall",
                            "fpr",
                            "f1",
                            "average_precision",
                            "expected_cost_usd",
                            "cost_per_1k_usd",
                            "counts",
                            "pr_curve",
                            "recall_at_fpr_0.05",
                            "precision_at_recall_0.9",
                        )
                        if kk in verticals[vert][k]["before"]
                    },
                    "after": {
                        kk: verticals[vert][k]["after"][kk]
                        for kk in (
                            "precision",
                            "recall",
                            "fpr",
                            "f1",
                            "average_precision",
                            "expected_cost_usd",
                            "cost_per_1k_usd",
                            "counts",
                            "pr_curve",
                            "recall_at_fpr_0.05",
                            "precision_at_recall_0.9",
                        )
                        if kk in verticals[vert][k]["after"]
                    },
                    "after_routed": {
                        kk: verticals[vert][k]["after_routed"][kk]
                        for kk in (
                            "precision",
                            "recall",
                            "fpr",
                            "expected_cost_usd",
                            "cost_per_1k_usd",
                        )
                        if kk in verticals[vert][k]["after_routed"]
                    },
                    "diagnosis": {
                        kk: verticals[vert][k]["diagnosis"][kk]
                        for kk in ("mean", "peak_gap", "cvar", "accuracy")
                    },
                }
                for k in ("ref", "benign_backtranslate", "benign_synonym", "malignant_attack")
            }
            for vert in verticals
        },
        "shape_detection": {
            "accuracy": acc,
            "confusion": conf,
            "n": len(shape_true),
            "y_true": shape_true,
            "y_pred": shape_pred,
        },
        "roi": {
            "cost_saved_usd_oracle_actuators": float(cost_saved),
            "cost_saved_usd_routed_actuators": float(cost_saved_routed),
            "malignant_mean_recall_lift": float(np.mean(recall_lifts)),
            "narrative": (
                "Oracle actuators show medicine value (right drug); routed actuators "
                "show end-to-end AS-DTA value (detect shape → pick drug). ROI = avoided "
                "FP/FN dollar cost on industrial LIVE."
            ),
        },
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "as_dta_business_pr.json").write_text(json.dumps(payload, indent=2))
    png = OUT / "as_dta_business_pr.png"
    plot_business_board(payload, png)
    _write_readme(OUT / "README.md", payload)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_business_PR_ROI.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(png, ARTIFACTS / "as_dta_business_pr.png")
    print("[biz] shape-detection acc=%.1f%%" % (100 * acc), flush=True)
    print("[biz] cost saved (oracle actuators): $%.0f" % cost_saved, flush=True)
    print("[biz] cost saved (routed actuators): $%.0f" % cost_saved_routed, flush=True)
    print("[biz] malignant mean recall lift: %+.3f" % np.mean(recall_lifts), flush=True)
    print("[biz] wrote", OUT, flush=True)
    return payload


def _write_readme(path: Path, payload: dict):
    lines = [
        "# AS-DTA business metrics — industrial HC3 + precision-recall",
        "",
        "Train on **general HC3**, LIVE on **finance / medicine** (industrial integrity",
        "verticals). Drift shape (Mean×Peak) routes the actuator; value is read on",
        "**precision / recall / FPR / expected cost**, not only accuracy.",
        "",
        "## Cost model",
        "",
        "- FP (human → flagged AI): **$%.0f** — complaint / trust / legal" % payload["cost_fp"],
        "- FN (AI missed): **$%.0f** — integrity / compliance miss" % payload["cost_fn"],
        "- Expected cost = $c_{FP}\\cdot FP + c_{FN}\\cdot FN$ on the LIVE batch",
        "",
        "## Drift-shape detection",
        "",
        "Accuracy **%.0f%%** (quiet / benign / malignant from Mean×Peak route)."
        % (100 * payload["shape_detection"]["accuracy"]),
        "",
        "```",
        json.dumps(payload["shape_detection"]["confusion"], indent=2),
        "```",
        "",
        "## Precision-Recall before → after",
        "",
        "| vertical | LIVE | recall | precision | FPR | cost/1k | actuator |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for vert, vo in payload["verticals"].items():
        for kind in ("benign_synonym", "malignant_attack"):
            b, a = vo[kind]["before"], vo[kind]["after"]
            lines.append(
                "| %s | %s | %.3f→**%.3f** | %.3f→%.3f | %.3f→%.3f | $%.0f→$%.0f | `%s` |"
                % (
                    vert,
                    kind,
                    b["recall"],
                    a["recall"],
                    b["precision"],
                    a["precision"],
                    b["fpr"],
                    a["fpr"],
                    b["cost_per_1k_usd"],
                    a["cost_per_1k_usd"],
                    vo[kind]["actuator"],
                )
            )
    lines.extend(
        [
            "",
            "## ROI rollup",
            "",
            "- Cost saved (oracle actuators): **$%.0f**" % payload["roi"]["cost_saved_usd_oracle_actuators"],
            "- Cost saved (routed / detected shape): **$%.0f**" % payload["roi"]["cost_saved_usd_routed_actuators"],
            "- Malignant mean recall lift: **%+.3f**" % payload["roi"]["malignant_mean_recall_lift"],
            "",
            payload["roi"]["narrative"],
            "",
            "## Run",
            "",
            "```bash",
            "python3 scripts/run_as_dta_business_pr.py --quick",
            "python3 scripts/run_as_dta_business_pr.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--max-ref", type=int, default=2500)
    ap.add_argument("--max-live", type=int, default=1200)
    ap.add_argument("--K", type=int, default=28)
    ap.add_argument("--cost-fp", type=float, default=COST_FP_DEFAULT)
    ap.add_argument("--cost-fn", type=float, default=COST_FN_DEFAULT)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    run_business(args)


if __name__ == "__main__":
    main()
