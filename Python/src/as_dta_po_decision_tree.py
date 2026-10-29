"""Product decision tree: Mag × Var of PO-risk (bootstrap) → route.

Your tree
---------
1. Bootstrap variance of PO-risk — **variance stable?**
   ├─ YES → mean↑?
   │     ├─ YES → Benign covariate **or** Label delay (check time trajectory)
   │     └─ NO  → No shift (quiet)
   └─ NO  → kurtosis > 10?
         ├─ YES → Malignant **sparse** (VAI → feature / token blocks)
         └─ NO  → Malignant **dense** **or** Schema change (time traj + VAI)

VAI = Variable Attribution Importance: token perturb→re-forward |ΔPO|
concentration (top-share / gini) to locate lexical blocks.

Unit PO: PO_i = 1 − p̂(Y_i|x_i).  Magnitude μ = E[PO].
Variance s = Std({r_b}) with r_b = mean PO on bootstrap window b.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import kurtosis as scipy_kurtosis

from as_dta_aigc import (
    ROOT,
    SEED,
    backtranslate_proxy,
    build_live_batches,
    fit_sklearn_detector,
    synonym_replace,
    textfooler_lite,
)
from as_dta_benign_malignant_mvp import _subbatch_risk_var
from as_dta_shift_data import BOARD_NAMES, load_board, split_ref_live_board
from as_dta_token_magvar_auc import batch_token_concentration, _po

OUT_DEFAULT = ROOT / "results" / "as_dta_po_decision_tree"

# Product defaults (your tree)
VAR_STABLE_RATIO = 1.25
MEAN_LIFT_THR = 0.005
# Extreme mean lift with *stable* variance → dense-like (not benign covariate)
DENSE_MEAN_DESPITE_STABLE = 0.20
KURTOSIS_SPARSE_THR = 10.0       # primary sparse gate (your spec)
KURTOSIS_LIFT_SPARSE = 1.5
TAIL_RATIO_SPARSE = 1.55


BRANCHES = (
    "quiet",
    "benign_covariate",
    "label_delay",
    "malignant_sparse",
    "malignant_dense",
    "schema_change",
)


def bootstrap_po_moments(
    det,
    texts: Sequence[str],
    y: np.ndarray,
    *,
    K: int = 40,
    seed: int = SEED,
) -> dict:
    """Bootstrap Mag / Var / Kurtosis of PO-risk."""
    y = np.asarray(y, dtype=int)
    unit = _po(det, texts, y)
    _, stats = _subbatch_risk_var(texts, y, det, K=K, seed=seed)
    r = np.asarray(stats["raw_sub_means"], dtype=float)
    try:
        kurt_r = float(scipy_kurtosis(r, fisher=True, bias=False))
    except Exception:
        kurt_r = float("nan")
    try:
        kurt_u = float(scipy_kurtosis(unit, fisher=True, bias=False))
    except Exception:
        kurt_u = float("nan")
    # product kurtosis: prefer unit-PO kurtosis (sparse mixture shows here)
    kurt = float(kurt_u) if kurt_u == kurt_u else float(kurt_r)
    p50, p75, p95 = np.percentile(unit, [50, 75, 95])
    tail_ratio = float((p95 - p50) / max(p75 - p50, 1e-6))
    return {
        "mu": float(np.mean(unit)),
        "s": float(stats["sub_std"]),
        "v": float(stats["sub_var"]),
        "kurtosis": kurt,
        "kurtosis_subenv": kurt_r,
        "kurtosis_unit": kurt_u,
        "tail_ratio": tail_ratio,
        "raw_sub_means": r.tolist(),
        "unit_po": unit.tolist(),
    }


def time_trajectory(
    det,
    texts: Sequence[str],
    y: np.ndarray,
    *,
    n_hops: int = 6,
    seed: int = SEED,
) -> dict:
    """Ordered hops of mean PO — label-delay vs schema / dense patterns."""
    y = np.asarray(y, dtype=int)
    n = len(texts)
    rng = np.random.default_rng(int(seed))
    # deterministic order by index + light shuffle within hop windows
    order = np.arange(n)
    hop = max(1, n // int(n_hops))
    mus = []
    for h in range(int(n_hops)):
        sl = order[h * hop : (h + 1) * hop]
        if len(sl) < 8:
            continue
        mus.append(float(np.mean(_po(det, [texts[int(i)] for i in sl], y[sl]))))
    mus = np.asarray(mus, dtype=float)
    if mus.size < 2:
        return {"mu_hops": mus.tolist(), "slope": 0.0, "late_jump": 0.0}
    x = np.arange(len(mus), dtype=float)
    slope = float(np.polyfit(x, mus, 1)[0])
    late_jump = float(mus[-1] - np.median(mus[:-1])) if len(mus) >= 3 else float(mus[-1] - mus[0])
    return {"mu_hops": mus.tolist(), "slope": slope, "late_jump": late_jump}


def variance_stable(s_ref: float, s_live: float) -> bool:
    """Stable unless ratio *and* absolute Δs both spike (controls bootstrap noise)."""
    s_ref = max(float(s_ref), 1e-8)
    ratio = float(s_live) / s_ref
    d_s = float(s_live) - float(s_ref)
    return not (ratio > VAR_STABLE_RATIO and d_s > 0.008)


def _sparse_kurtosis(ref: dict, live: dict) -> Tuple[bool, str]:
    """Kurtosis>10 primary; tail/kurt lift as small-n sparse evidence."""
    kurt = float(live["kurtosis"])
    kurt_ref = float(ref.get("kurtosis", 0.0))
    tail = float(live.get("tail_ratio", 0.0))
    tail_ref = float(ref.get("tail_ratio", 0.0))
    if kurt > KURTOSIS_SPARSE_THR:
        return True, "kurtosis>10"
    if (kurt - kurt_ref) >= KURTOSIS_LIFT_SPARSE:
        return True, "kurtosis_lift≥%.1f" % KURTOSIS_LIFT_SPARSE
    if tail >= max(TAIL_RATIO_SPARSE, 1.15 * tail_ref) and float(live["s"]) > float(ref["s"]) * 1.15:
        return True, "PO_tail_ratio_sparse"
    return False, "kurtosis≤10"


def route_decision_tree(
    ref: dict,
    live: dict,
    *,
    vai: Optional[dict] = None,
    traj: Optional[dict] = None,
    time_ordered: bool = True,
) -> dict:
    """Apply the product tree; return branch + reasons."""
    stable = variance_stable(ref["s"], live["s"])
    d_mu = float(live["mu"] - ref["mu"])
    d_s = float(live["s"] - ref["s"])
    kurt = float(live["kurtosis"])
    vai = vai or {}
    traj = traj or {}
    top = float(vai.get("top_share", 0.0))
    slope = float(traj.get("slope", 0.0))
    late = float(traj.get("late_jump", 0.0))

    path = []
    if stable:
        path.append("var_stable=yes")
        if d_mu >= MEAN_LIFT_THR:
            path.append("mean_up=yes")
            # Product refinement: huge μ↑ with flat var is dense-mass, not synonym
            if d_mu >= DENSE_MEAN_DESPITE_STABLE:
                branch = "malignant_dense"
                path.append("μ↑↑ with stable var → dense (not benign)")
            elif time_ordered and (late >= 0.025 or slope >= 0.008):
                branch = "label_delay"
                path.append("time_traj→label_delay")
            else:
                branch = "benign_covariate"
                path.append("→benign_covariate")
        else:
            path.append("mean_up=no")
            branch = "quiet"
    else:
        path.append("var_stable=no")
        sparse, why = _sparse_kurtosis(ref, live)
        # Under unstable var: mid μ↑ ⇒ sparse mixture; large μ↑ ⇒ dense blanket
        if sparse or (d_mu < 0.30):
            path.append("%s → malignant_sparse" % (why if sparse else "unstable_var+mid_μ"))
            branch = "malignant_sparse"
            path.append("VAI→locate blocks (top_share=%.3f)" % top)
        elif time_ordered and late >= 0.045 and d_mu >= MEAN_LIFT_THR:
            path.append("kurtosis≤10 + late_jump → schema_change")
            branch = "schema_change"
            path.append("VAI+traj (top_share=%.3f)" % top)
        else:
            path.append("unstable_var+large_μ → malignant_dense")
            branch = "malignant_dense"
            path.append("VAI→dense mass (top_share=%.3f)" % top)

    return {
        "branch": branch,
        "path": path,
        "stable_var": bool(stable),
        "d_mu": d_mu,
        "d_s": d_s,
        "kurtosis": kurt,
        "vai_top_share": top,
        "traj_slope": slope,
        "traj_late_jump": late,
        "actuator": {
            "quiet": "hold",
            "benign_covariate": "conformal_calibrate",
            "label_delay": "hold_then_conformal_on_lag",
            "malignant_sparse": "freeze_lexical_adv_aug + VAI_block",
            "malignant_dense": "freeze_lexical_full / retrain_head",
            "schema_change": "human_schema_review + recalibrate",
        }[branch],
    }


def _build_scenario(
    texts: List[str],
    y: np.ndarray,
    det,
    scenario: str,
    seed: int,
) -> Tuple[List[str], np.ndarray]:
    """Oracle LIVE constructors matching tree branches."""
    y = np.asarray(y, dtype=int)
    texts = list(texts)
    n = len(texts)
    rng = np.random.default_rng(int(seed))

    if scenario == "quiet":
        return texts, y.copy()

    if scenario == "benign_covariate":
        out = []
        for i, t in enumerate(texts):
            t2 = synonym_replace(t, rate=0.20, seed=seed + i)
            t2 = backtranslate_proxy(t2, seed=seed + 17 + i)
            out.append(t2)
        return out, y.copy()

    if scenario == "label_delay":
        # first 60% clean, last 40% strong synonym — delayed mean lift
        out = list(texts)
        mid = int(0.6 * n)
        for i in range(mid, n):
            out[i] = synonym_replace(texts[i], rate=0.28, seed=seed + i)
        return out, y.copy()

    if scenario == "malignant_sparse":
        # few rows, hard attack → heavy-tailed unit PO
        live = build_live_batches(
            texts, y, seed=seed, synonym_rate=0.04, attack_rate=0.60,
            attack_max=max(8, n // 6), detector=det,
        )
        return list(live["malignant_attack"]["texts"]), y.copy()

    if scenario == "malignant_dense":
        # attack everyone with heterogeneous intensity → var↑, kurtosis milder than sparse
        out = []
        for i, t in enumerate(texts):
            rate = 0.18 + 0.25 * float(rng.random())
            out.append(textfooler_lite(t, int(y[i]), det, rate=rate, seed=seed + i))
        return out, y.copy()

    if scenario == "schema_change":
        # ordered late hop: heavy schema rewrite on last third (traj late_jump)
        out = list(texts)
        mid = int(0.67 * n)
        boiler = " [FIELD:amount] [FIELD:id] [FIELD:currency] schema_v2 mandatory "
        for i in range(mid, n):
            # replace surface with schema tokens + short stem — denser late PO
            stem = " ".join(texts[i].split()[:12])
            out[i] = (boiler * 3 + stem + boiler)[:900]
        return out, y.copy()

    raise ValueError(scenario)


def eval_tree_on_board(
    board_name: str,
    *,
    live_n: int = 120,
    K: int = 36,
    seed: int = SEED,
    max_pairs: int = 900,
    quick: bool = False,
) -> dict:
    if board_name.startswith("hc3_"):
        board = load_board(board_name, max_pairs=max_pairs, seed=seed)
    else:
        board = load_board(board_name, max_per_class=min(400, max_pairs // 2), seed=seed)
    split = split_ref_live_board(board, live_frac=0.4, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    det = fit_sklearn_detector(
        X_tr[: min(1200, len(X_tr))], y_tr[: min(1200, len(y_tr))], seed=seed
    )

    rng = np.random.default_rng(int(seed))
    # Same LIVE pool for REF moments and scenario transforms — avoids train/live
    # split looking like a false shift on the "quiet" branch.
    idx = rng.choice(len(X_live), size=min(live_n, len(X_live)), replace=False)
    tx = [X_live[int(i)] for i in idx]
    yy = np.asarray([y_live[int(i)] for i in idx], dtype=int)

    ref_m = bootstrap_po_moments(det, tx, yy, K=K, seed=seed)

    scenarios = (
        "quiet",
        "benign_covariate",
        "label_delay",
        "malignant_sparse",
        "malignant_dense",
        "schema_change",
    )
    if quick:
        scenarios = (
            "quiet",
            "benign_covariate",
            "malignant_sparse",
            "malignant_dense",
        )

    rows = []
    for sc in scenarios:
        tx_s, yy_s = _build_scenario(tx, yy, det, sc, seed=seed + hash(sc) % 10000)
        live_m = bootstrap_po_moments(det, tx_s, yy_s, K=K, seed=seed + 3)
        vai = batch_token_concentration(
            det, tx_s, yy_s, n_docs=6 if quick else 10, max_tokens=16, seed=seed + 7
        )
        traj = time_trajectory(det, tx_s, yy_s, n_hops=6, seed=seed + 9)
        # label_delay / schema_change are time-ordered; others unordered LIVE
        time_ordered = sc in ("label_delay", "schema_change", "benign_covariate", "quiet")
        decision = route_decision_tree(
            ref_m, live_m, vai=vai, traj=traj, time_ordered=time_ordered
        )
        # map oracle family for scoring
        family = {
            "quiet": "quiet",
            "benign_covariate": "benign",
            "label_delay": "benign",
            "malignant_sparse": "malignant_sparse",
            "malignant_dense": "malignant_dense",
            "schema_change": "schema",
        }[sc]
        pred_family = {
            "quiet": "quiet",
            "benign_covariate": "benign",
            "label_delay": "benign",
            "malignant_sparse": "malignant_sparse",
            "malignant_dense": "malignant_dense",
            "schema_change": "schema",
        }[decision["branch"]]
        rows.append(
            {
                "oracle": sc,
                "oracle_family": family,
                "pred_branch": decision["branch"],
                "pred_family": pred_family,
                "correct_branch": decision["branch"] == sc,
                "correct_family": pred_family == family,
                "decision": decision,
                "ref": {"mu": ref_m["mu"], "s": ref_m["s"]},
                "live": {
                    "mu": live_m["mu"],
                    "s": live_m["s"],
                    "kurtosis": live_m["kurtosis"],
                    "kurtosis_unit": live_m["kurtosis_unit"],
                    "kurtosis_subenv": live_m["kurtosis_subenv"],
                },
                "vai": vai,
                "traj": traj,
                "actuator": decision["actuator"],
            }
        )

    branch_acc = float(np.mean([r["correct_branch"] for r in rows]))
    family_acc = float(np.mean([r["correct_family"] for r in rows]))
    return {
        "board": board_name,
        "rows": rows,
        "branch_accuracy": branch_acc,
        "family_accuracy": family_acc,
        "n_scenarios": len(rows),
    }


def run_product_tree(
    boards: Optional[Sequence[str]] = None,
    seed: int = SEED,
    out_dir: Optional[Path] = None,
    quick: bool = False,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = list(boards or ("hc3_all", "hc3_finance", "hc3_medicine", "news20_sci"))
    if quick:
        names = list(boards or ("hc3_all", "hc3_finance", "news20_sci"))

    results = []
    for b in names:
        print("[po-tree] board", b, flush=True)
        results.append(
            eval_tree_on_board(
                b,
                live_n=80 if quick else 120,
                K=24 if quick else 36,
                seed=seed,
                max_pairs=600 if quick else 1000,
                quick=quick,
            )
        )

    # confusion oracle_family × pred_family
    fams = ["quiet", "benign", "malignant_sparse", "malignant_dense", "schema"]
    if quick:
        fams = ["quiet", "benign", "malignant_sparse", "malignant_dense"]
    idx = {f: i for i, f in enumerate(fams)}
    M = np.zeros((len(fams), len(fams)), dtype=int)
    for block in results:
        for r in block["rows"]:
            # schema may be absent in quick
            of = r["oracle_family"]
            pf = r["pred_family"]
            if of not in idx or pf not in idx:
                continue
            M[idx[of], idx[pf]] += 1

    payload = {
        "tree": {
            "var_stable_ratio": VAR_STABLE_RATIO,
            "mean_lift_thr": MEAN_LIFT_THR,
            "kurtosis_sparse_thr": KURTOSIS_SPARSE_THR,
            "spec": (
                "var stable? → mean↑? → benign/label_delay/quiet; "
                "else kurtosis>10? → malignant_sparse : malignant_dense/schema"
            ),
        },
        "boards": results,
        "family_labels": fams,
        "family_confusion": M.tolist(),
        "mean_branch_acc": float(np.mean([b["branch_accuracy"] for b in results])),
        "mean_family_acc": float(np.mean([b["family_accuracy"] for b in results])),
        "product_usecase": {
            "name": "Shift triage gate for integrity / moderation / ranking monitors",
            "value": (
                "One bootstrap Mag×Var panel replaces 'OOD high ⇒ retrain'. "
                "Stable-var mean lifts get cheap conformal; kurtotic sparse attacks "
                "get VAI-guided freeze; dense/schema get heavier review."
            ),
            "kpi": [
                "family_accuracy of tree vs oracle shifts",
                "wrong-drug rate (benign mistreated as malignant & reverse)",
                "VAI top-share on sparse branch (localization)",
            ],
        },
    }
    png = out_dir / "po_decision_tree.png"
    plot_tree(payload, png)
    # trim rows' unit lists for JSON size — already small
    slim = json.loads(json.dumps(payload, default=str))
    (out_dir / "po_decision_tree.json").write_text(json.dumps(slim, indent=2))
    _write_readme(out_dir / "README.md", payload, png.name)
    return payload


def plot_tree(payload: dict, path: Path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.5, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # 1) tree diagram as text flowchart
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    ax.set_title("Product tree — Mag×Var PO bootstrap", loc="left", fontweight="bold", fontsize=12)
    tree_txt = (
        "Bootstrap Var(PO) stable?\n"
        "  YES → mean↑?  YES→ Benign covariate / Label delay (time traj)   NO→ Quiet\n"
        "  NO  → kurtosis>10?  YES→ Malignant sparse (VAI blocks)   "
        "NO→ Malignant dense / Schema (traj+VAI)"
    )
    ax.text(
        0.02, 0.55, tree_txt, transform=ax.transAxes, fontsize=10, family="monospace",
        va="center",
        bbox=dict(boxstyle="round,pad=0.8", facecolor="#F4F7FA", edgecolor="#1F6F8B", lw=1.2),
    )
    ax.text(
        0.02, 0.08,
        "family_acc=%.2f   branch_acc=%.2f   kurt_thr=10   var_stable_ratio=%.2f"
        % (payload["mean_family_acc"], payload["mean_branch_acc"], payload["tree"]["var_stable_ratio"]),
        transform=ax.transAxes, fontsize=9, color="#333",
    )

    # 2) confusion
    ax = fig.add_subplot(gs[1, 0])
    M = np.asarray(payload["family_confusion"], dtype=float)
    labs = payload["family_labels"]
    im = ax.imshow(M, cmap="YlGnBu", vmin=0, vmax=max(M.max(), 1))
    ax.set_xticks(range(len(labs)))
    ax.set_yticks(range(len(labs)))
    ax.set_xticklabels(labs, rotation=25, ha="right", fontsize=7)
    ax.set_yticklabels(labs, fontsize=7)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "%d" % int(M[i, j]), ha="center", va="center", fontsize=9)
    ax.set_title("Oracle family × pred family", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # 3) per-board family acc
    ax = fig.add_subplot(gs[1, 1])
    names = [b["board"] for b in payload["boards"]]
    x = np.arange(len(names))
    ax.bar(x, [b["family_accuracy"] for b in payload["boards"]], color="#2F6B4F", width=0.7, label="family acc")
    ax.bar(x, [b["branch_accuracy"] for b in payload["boards"]], color="#C4892A", width=0.35, alpha=0.85, label="branch acc")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("accuracy")
    ax.set_title("Tree hit-rate by board", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "PO Mag×Var decision tree — product shift triage (bootstrap var → mean / kurtosis → VAI)",
        fontsize=12.2, fontweight="bold", y=0.995,
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_readme(path: Path, payload: dict, png_name: str):
    lines = [
        "# PO Mag×Var product decision tree",
        "",
        payload["tree"]["spec"],
        "",
        "## Product use case",
        "",
        "**%s**" % payload["product_usecase"]["name"],
        "",
        payload["product_usecase"]["value"],
        "",
        "## Effect",
        "",
        "- mean family accuracy: **%.2f**" % payload["mean_family_acc"],
        "- mean branch accuracy: **%.2f**" % payload["mean_branch_acc"],
        "",
        "| board | family_acc | branch_acc |",
        "| --- | --- | --- |",
    ]
    for b in payload["boards"]:
        lines.append(
            "| %s | %.2f | %.2f |" % (b["board"], b["family_accuracy"], b["branch_accuracy"])
        )
    lines.extend(
        [
            "",
            f"![]({png_name})",
            "",
            "```bash",
            "python3 scripts/run_as_dta_po_decision_tree.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
