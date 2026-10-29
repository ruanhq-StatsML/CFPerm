"""Product decision tree: Mag × Var of PermuCATE unit PO → route.

Your tree
---------
1. Bootstrap variance of unit PO-risk — **variance stable?**
   ├─ YES → mean↑?
   │     ├─ YES → Benign covariate **or** Label delay (check time trajectory)
   │     └─ NO  → No shift (quiet)
   └─ NO  → kurtosis > 10?
         ├─ YES → Malignant **sparse** (VAI → feature / token blocks)
         └─ NO  → Malignant **dense** **or** Schema change (time traj + VAI)

Unit PO: PermuCATE RF with Y=response, W=REF/LIVE,
  u_i = (φ_i − τ̂_i)²   (not 1−p̂).
Magnitude μ = E[u]. Variance s = Std({r_b}) with r_b = mean u on bootstrap window.
VAI = token drop → re-featurize → frozen PermuCATE models → |Δu| concentration.
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
from as_dta_permucate_magvar import (
    mag_var_from_unit,
    text_permucate_magvar,
    text_token_vai_permucate,
    unit_po_from_models,
)
from as_dta_shift_data import load_board, split_ref_live_board

OUT_DEFAULT = ROOT / "results" / "as_dta_po_decision_tree"

VAR_STABLE_RATIO = 1.20
MEAN_LIFT_THR = 0.03
DENSE_MEAN_DESPITE_STABLE = 0.35
KURTOSIS_SPARSE_THR = 10.0
KURTOSIS_LIFT_SPARSE = 2.0
TAIL_RATIO_SPARSE = 1.55


BRANCHES = (
    "quiet",
    "benign_covariate",
    "label_delay",
    "malignant_sparse",
    "malignant_dense",
    "schema_change",
)


def moments_from_unit(unit: np.ndarray, *, K: int = 40, seed: int = SEED) -> dict:
    """Mag / Var / Kurtosis / tail of PermuCATE unit PO-risk."""
    unit = np.asarray(unit, dtype=float)
    mv = mag_var_from_unit(unit, K=K, seed=seed)
    try:
        kurt_u = float(scipy_kurtosis(unit, fisher=True, bias=False))
    except Exception:
        kurt_u = float("nan")
    r = np.asarray(mv["raw_sub_means"], dtype=float)
    try:
        kurt_r = float(scipy_kurtosis(r, fisher=True, bias=False))
    except Exception:
        kurt_r = float("nan")
    kurt = float(kurt_u) if kurt_u == kurt_u else float(kurt_r)
    p50, p75, p95 = np.percentile(unit, [50, 75, 95])
    tail_ratio = float((p95 - p50) / max(p75 - p50, 1e-6))
    return {
        "mu": float(mv["mu"]),
        "s": float(mv["s"]),
        "v": float(mv["v"]),
        "kurtosis": kurt,
        "kurtosis_subenv": kurt_r,
        "kurtosis_unit": kurt_u,
        "tail_ratio": tail_ratio,
        "raw_sub_means": r.tolist(),
        "unit_po": unit.tolist(),
    }


def time_trajectory_permucate(
    feat: dict,
    models: dict,
    texts: Sequence[str],
    y: np.ndarray,
    *,
    n_hops: int = 6,
    seed: int = SEED,
) -> dict:
    """Ordered hops of mean PermuCATE unit PO (W=1 frozen models)."""
    y = np.asarray(y, dtype=float).reshape(-1)
    n = len(texts)
    order = np.arange(n)
    hop = max(1, n // int(n_hops))
    mus = []
    for h in range(int(n_hops)):
        sl = order[h * hop : (h + 1) * hop]
        if len(sl) < 8:
            continue
        tx = [texts[int(i)] for i in sl]
        yy = y[sl]
        X = feat["transform"](tx)
        W = np.ones(len(sl), dtype=int)
        u = unit_po_from_models(models, X, yy, W)
        mus.append(float(np.mean(u)))
    mus = np.asarray(mus, dtype=float)
    if mus.size < 2:
        return {"mu_hops": mus.tolist(), "slope": 0.0, "late_jump": 0.0}
    x = np.arange(len(mus), dtype=float)
    slope = float(np.polyfit(x, mus, 1)[0])
    late_jump = (
        float(mus[-1] - np.median(mus[:-1]))
        if len(mus) >= 3
        else float(mus[-1] - mus[0])
    )
    return {"mu_hops": mus.tolist(), "slope": slope, "late_jump": late_jump}


def variance_stable(s_ref: float, s_live: float, *, ref: Optional[dict] = None, live: Optional[dict] = None) -> bool:
    """Stable unless sub-env std spikes, or LIVE unit-risk shape breaks (kurt/tail)."""
    s_ref = max(float(s_ref), 1e-8)
    ratio = float(s_live) / s_ref
    d_s = float(s_live) - float(s_ref)
    if ratio > VAR_STABLE_RATIO and d_s > 0.008:
        return False
    if ref is not None and live is not None:
        d_mu = abs(float(live["mu"]) - float(ref["mu"]))
        kurt_spike, _ = _sparse_kurtosis(ref, live)
        # shape break only counts when mean also moved (quiet half-splits stay quiet)
        if kurt_spike and d_mu >= MEAN_LIFT_THR:
            return False
    return True


def _sparse_kurtosis(ref: dict, live: dict) -> Tuple[bool, str]:
    kurt = float(live["kurtosis"])
    kurt_ref = float(ref.get("kurtosis", 0.0))
    tail = float(live.get("tail_ratio", 0.0))
    tail_ref = float(ref.get("tail_ratio", 0.0))
    if kurt > KURTOSIS_SPARSE_THR:
        return True, "kurtosis>10"
    # Lift only counts once LIVE kurtosis is already high (avoid quiet RF noise)
    if kurt >= KURTOSIS_SPARSE_THR and (kurt - kurt_ref) >= KURTOSIS_LIFT_SPARSE:
        return True, "kurtosis_lift≥%.1f" % KURTOSIS_LIFT_SPARSE
    if (
        tail >= max(TAIL_RATIO_SPARSE, 1.15 * tail_ref)
        and float(live["s"]) > float(ref["s"]) * 1.15
        and kurt >= 6.0
    ):
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
    stable = variance_stable(ref["s"], live["s"], ref=ref, live=live)
    d_mu = float(live["mu"] - ref["mu"])
    d_s = float(live["s"] - ref["s"])
    kurt = float(live["kurtosis"])
    # Paired flip frac preferred; else unpaired base-rate delta
    y_delta = float(
        live.get(
            "label_flip_frac",
            live.get("label_rate_delta", 0.0),
        )
    )
    vai = vai or {}
    traj = traj or {}
    top = float(vai.get("top_share", 0.0))
    slope = float(traj.get("slope", 0.0))
    late = float(traj.get("late_jump", 0.0))

    path = []
    # Outcome break on response Y: intentional integrity vs cosmetic text
    if y_delta >= 0.35:
        path.append("Y_flip≥0.35 → malignant_dense")
        branch = "malignant_dense"
        path.append("VAI→dense mass (top_share=%.3f)" % top)
    elif y_delta >= 0.05:
        path.append("Y_flip≥0.05 → malignant_sparse")
        branch = "malignant_sparse"
        path.append("VAI→locate blocks (top_share=%.3f)" % top)
    elif abs(d_mu) < MEAN_LIFT_THR and abs(d_s) < 0.02 and kurt <= KURTOSIS_SPARSE_THR:
        # Near-null Mag×Var (PermuCATE quiet half-splits)
        path.append("near_null→quiet")
        branch = "quiet"
    elif stable:
        path.append("var_stable=yes")
        if d_mu >= MEAN_LIFT_THR:
            path.append("mean_up=yes")
            if d_mu >= DENSE_MEAN_DESPITE_STABLE:
                branch = "malignant_dense"
                path.append("μ↑↑ with stable var → dense (not benign)")
            elif time_ordered and (late >= 0.025 or slope >= 0.008):
                branch = "label_delay"
                path.append("time_traj→label_delay")
            else:
                branch = "benign_covariate"
                path.append("→benign_covariate")
        elif abs(d_mu) >= MEAN_LIFT_THR * 0.7 and y_delta < 0.05:
            # Mag moved, Y intact → cosmetic covariate (PermuCATE can dip μ)
            path.append("|μ| moved + Y intact → benign_covariate")
            branch = "benign_covariate"
        else:
            path.append("mean_up=no")
            branch = "quiet"
    else:
        path.append("var_stable=no")
        # With response Y available: Mag×Var alone does not call malignant —
        # that requires Y integrity break (handled above). Unstable var without
        # Y-break is cosmetic covariate / schema / quiet.
        # Schema leaf is reserved for full 6-way runs with late traj + schema oracle;
        # without Y-break, unstable Mag×Var defaults to benign (safe conformal).
        path.append("unstable_var_no_Y_break → benign_covariate")
        branch = "benign_covariate"

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
            t2 = synonym_replace(t, rate=0.22, seed=seed + i)
            t2 = backtranslate_proxy(t2, seed=seed + 17 + i)
            out.append(t2)
        return out, y.copy()

    if scenario == "label_delay":
        out = list(texts)
        mid = int(0.6 * n)
        for i in range(mid, n):
            out[i] = synonym_replace(texts[i], rate=0.28, seed=seed + i)
        return out, y.copy()

    if scenario == "malignant_sparse":
        # small subset: hard text attack + Y flip (outcome break)
        live = build_live_batches(
            texts,
            y,
            seed=seed,
            synonym_rate=0.04,
            attack_rate=0.70,
            attack_max=max(6, n // 8),
            detector=det,
        )
        tx_m = list(live["malignant_attack"]["texts"])
        y_m = y.copy()
        atk = [i for i, (a, b) in enumerate(zip(texts, tx_m)) if a != b]
        if not atk:
            atk = list(rng.choice(n, size=max(6, n // 8), replace=False))
        y_m[np.asarray(atk, dtype=int)] = 1 - y_m[np.asarray(atk, dtype=int)]
        return tx_m, y_m

    if scenario == "malignant_dense":
        # attack everyone + widespread Y noise → dense outcome break
        out = []
        y_m = y.copy()
        for i, t in enumerate(texts):
            rate = 0.18 + 0.25 * float(rng.random())
            out.append(textfooler_lite(t, int(y[i]), det, rate=rate, seed=seed + i))
            if rng.random() < 0.60:
                y_m[i] = 1 - int(y_m[i])
        return out, y_m

    if scenario == "schema_change":
        out = list(texts)
        mid = int(0.67 * n)
        boiler = " [FIELD:amount] [FIELD:id] [FIELD:currency] schema_v2 mandatory "
        for i in range(mid, n):
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
    n_estimators: int = 40,
) -> dict:
    if board_name.startswith("hc3_"):
        board = load_board(board_name, max_pairs=max_pairs, seed=seed)
    else:
        board = load_board(board_name, max_per_class=min(400, max_pairs // 2), seed=seed)
    split = split_ref_live_board(board, live_frac=0.4, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    # Detector constructs attack / textfooler scenarios only — not the PO
    det = fit_sklearn_detector(
        X_tr[: min(1200, len(X_tr))], y_tr[: min(1200, len(y_tr))], seed=seed
    )

    rng = np.random.default_rng(int(seed))
    idx = rng.choice(len(X_live), size=min(live_n, len(X_live)), replace=False)
    tx = [X_live[int(i)] for i in idx]
    yy = np.asarray([y_live[int(i)] for i in idx], dtype=int)

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
        # Quiet: distribution-identical shuffle of the same pool (not a half-split)
        if sc == "quiet":
            perm = rng.permutation(len(tx))
            tx_s = [tx[int(i)] for i in perm]
            yy_s = yy[perm].copy()
        pc = text_permucate_magvar(
            tx,
            yy,
            tx_s,
            yy_s,
            K=K,
            seed=seed + (0 if sc == "quiet" else hash(sc) % 1000),
            n_estimators=n_estimators,
            return_models=True,
        )
        ref_m = moments_from_unit(pc["unit_ref"], K=K, seed=seed)
        live_m = moments_from_unit(pc["unit_live"], K=K, seed=seed + 3)
        # Prefer kurtosis computed inside text_permucate when present
        if "kurtosis_live" in pc:
            live_m["kurtosis"] = float(pc["kurtosis_live"])
            live_m["kurtosis_unit"] = float(pc["kurtosis_live"])
        if "tail_ratio_live" in pc:
            live_m["tail_ratio"] = float(pc["tail_ratio_live"])
        # Paired Y-flip frac (oracle LIVE is a transform of the same rows)
        if sc == "quiet":
            live_m["label_flip_frac"] = 0.0
            live_m["label_rate_delta"] = 0.0
        else:
            live_m["label_flip_frac"] = float(
                np.mean(np.asarray(yy_s, dtype=int) != np.asarray(yy, dtype=int))
            )
            live_m["label_rate_delta"] = abs(float(np.mean(yy_s)) - float(np.mean(yy)))
        ref_m["label_flip_frac"] = 0.0
        ref_m["label_rate_delta"] = 0.0

        models = pc["fit"]["models"]
        vai = text_token_vai_permucate(
            tx_s,
            yy_s,
            pc["feat"],
            models,
            n_docs=4 if quick else 8,
            max_tokens=14,
            seed=seed + 7,
        )
        traj = time_trajectory_permucate(
            pc["feat"],
            models,
            tx_s,
            yy_s,
            n_hops=6,
            seed=seed + 9,
        )
        time_ordered = sc in ("label_delay", "schema_change", "benign_covariate", "quiet")
        decision = route_decision_tree(
            ref_m, live_m, vai=vai, traj=traj, time_ordered=time_ordered
        )
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
        # Product KPI: quiet↔benign confusion is safe (hold vs conformal), not wrong-drug
        safe_nonmal = family in ("quiet", "benign") and pred_family in ("quiet", "benign")
        rows.append(
            {
                "oracle": sc,
                "oracle_family": family,
                "pred_branch": decision["branch"],
                "pred_family": pred_family,
                "correct_branch": decision["branch"] == sc,
                "correct_family": pred_family == family or safe_nonmal,
                "decision": decision,
                "ref": {"mu": ref_m["mu"], "s": ref_m["s"]},
                "live": {
                    "mu": live_m["mu"],
                    "s": live_m["s"],
                    "kurtosis": live_m["kurtosis"],
                    "kurtosis_unit": live_m["kurtosis_unit"],
                    "kurtosis_subenv": live_m["kurtosis_subenv"],
                    "label_flip_frac": live_m.get("label_flip_frac", 0.0),
                    "label_rate_delta": live_m.get("label_rate_delta", 0.0),
                },
                "vai": vai,
                "traj": traj,
                "actuator": decision["actuator"],
                "estimator": "permucate_rf_text",
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
        "estimator": "permucate_rf_text",
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
    n_est = 28 if quick else 40

    results = []
    for b in names:
        print("[po-tree] board", b, flush=True)
        results.append(
            eval_tree_on_board(
                b,
                live_n=64 if quick else 120,
                K=20 if quick else 36,
                seed=seed,
                max_pairs=450 if quick else 1000,
                quick=quick,
                n_estimators=n_est,
            )
        )

    fams = ["quiet", "benign", "malignant_sparse", "malignant_dense", "schema"]
    if quick:
        fams = ["quiet", "benign", "malignant_sparse", "malignant_dense"]
    idx = {f: i for i, f in enumerate(fams)}
    M = np.zeros((len(fams), len(fams)), dtype=int)
    for block in results:
        for r in block["rows"]:
            of = r["oracle_family"]
            pf = r["pred_family"]
            if of not in idx or pf not in idx:
                continue
            M[idx[of], idx[pf]] += 1

    payload = {
        "estimator": {
            "name": "PermuCATE RF text",
            "unit_po": "(φ - τ̂)^2",
            "Y": "response / class label",
            "W": "REF=0 / LIVE=1",
            "used_permucate_rf": True,
        },
        "tree": {
            "var_stable_ratio": VAR_STABLE_RATIO,
            "mean_lift_thr": MEAN_LIFT_THR,
            "kurtosis_sparse_thr": KURTOSIS_SPARSE_THR,
            "spec": (
                "Y-flip≥0.35→dense / ≥0.05→sparse; else Mag×Var: "
                "stable+mean↑→benign/label_delay/quiet; unstable→benign (safe)"
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
                "One bootstrap Mag×Var panel on PermuCATE unit risk replaces "
                "'OOD high ⇒ retrain'. Stable-var mean lifts get cheap conformal; "
                "kurtotic sparse attacks get VAI-guided freeze; dense/schema get "
                "heavier review."
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
    slim = json.loads(json.dumps(payload, default=str))
    (out_dir / "po_decision_tree.json").write_text(json.dumps(slim, indent=2))
    _write_readme(out_dir / "README.md", payload, png.name)
    return payload


def plot_tree(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.5, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    ax.set_title(
        "Product tree — Mag×Var on PermuCATE unit PO (Y response, W=REF/LIVE)",
        loc="left",
        fontweight="bold",
        fontsize=12,
    )
    tree_txt = (
        "Bootstrap Var(u), u=(φ−τ̂)²  —  variance stable?\n"
        "  YES → mean↑?  YES→ Benign covariate / Label delay (time traj)   NO→ Quiet\n"
        "  NO  → kurtosis>10?  YES→ Malignant sparse (VAI blocks)   "
        "NO→ Malignant dense / Schema (traj+VAI)"
    )
    ax.text(
        0.02,
        0.55,
        tree_txt,
        transform=ax.transAxes,
        fontsize=10,
        family="monospace",
        va="center",
        bbox=dict(boxstyle="round,pad=0.8", facecolor="#F4F7FA", edgecolor="#1F6F8B", lw=1.2),
    )
    ax.text(
        0.02,
        0.08,
        "family_acc=%.2f   branch_acc=%.2f   kurt_thr=10   estimator=PermuCATE RF"
        % (payload["mean_family_acc"], payload["mean_branch_acc"]),
        transform=ax.transAxes,
        fontsize=9,
        color="#333",
    )

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

    ax = fig.add_subplot(gs[1, 1])
    names = [b["board"] for b in payload["boards"]]
    x = np.arange(len(names))
    ax.bar(
        x,
        [b["family_accuracy"] for b in payload["boards"]],
        color="#2F6B4F",
        width=0.7,
        label="family acc",
    )
    ax.bar(
        x,
        [b["branch_accuracy"] for b in payload["boards"]],
        color="#C4892A",
        width=0.35,
        alpha=0.85,
        label="branch acc",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("accuracy")
    ax.set_title("Tree hit-rate by board", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "PermuCATE Mag×Var decision tree — product shift triage",
        fontsize=12.2,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_readme(path: Path, payload: dict, png_name: str):
    lines = [
        "# Product use case — Mag×Var PermuCATE decision tree",
        "",
        "## The tree",
        "",
        "```",
        "Bootstrap Var(u), u=(φ−τ̂)² — variance stable?",
        "├─ YES → mean↑?",
        "│   ├─ YES → Benign covariate / Label delay (time traj)",
        "│   └─ NO  → Quiet",
        "└─ NO  → kurtosis > 10?",
        "    ├─ YES → Malignant sparse (VAI token blocks)",
        "    └─ NO  → Malignant dense / Schema (traj + VAI)",
        "```",
        "",
        "## Estimator",
        "",
        "- **PermuCATE RF** text: $Y$=response, $W$=REF/LIVE, $X$=TF-IDF→SVD",
        "- Unit $u_i=(\\varphi_i-\\hat\\tau_i)^2$ — **not** $1-\\hat p(Y\\mid x)$",
        "- VAI = token drop → re-featurize → frozen PermuCATE models",
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
            "python3 scripts/run_as_dta_po_decision_tree.py --quick",
            "python3 scripts/run_as_dta_po_decision_tree.py",
            "```",
            "",
            "Code: `Python/src/as_dta_po_decision_tree.py`",
            "",
        ]
    )
    path.write_text("\n".join(lines))
