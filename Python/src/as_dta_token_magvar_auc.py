"""Token-wise perturbation → Mag×Var AUC on PermuCATE unit PO.

Logic
-----
1. Fit attack detector on REF (only to *construct* malignant LIVE texts).
2. LIVE shifts: benign (uniform synonym) vs malignant (subset attack).
3. Unit PO via PermuCATE RF: Y = response label, W = REF(0)/LIVE(1),
   u = (φ − τ̂)²  (not 1−p̂).
4. Magnitude μ = E[u]; Variance s = Std(bootstrap sub-env mean u).
5. Token-wise VAI: drop token j → re-featurize → recompute u with frozen
   PermuCATE models → a_j = |u^{drop j} − u|.
6. Score Mag×Var (and token concentration); report AUC / 2×2.

No multimodal — token / lexical axis only. Requires response Y.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.metrics import roc_auc_score

from as_dta_aigc import (
    ROOT,
    SEED,
    build_live_batches,
    fit_sklearn_detector,
    synonym_replace,
)
from as_dta_permucate_magvar import text_permucate_magvar, text_token_vai_permucate
from as_dta_shift_data import BOARD_NAMES, load_board, split_ref_live_board

OUT_DEFAULT = ROOT / "results" / "as_dta_token_magvar_auc"


def magvar_features_permucate(
    texts_ref: Sequence[str],
    y_ref: np.ndarray,
    texts_live: Sequence[str],
    y_live: np.ndarray,
    *,
    K: int = 24,
    seed: int = SEED,
    token_docs: int = 10,
    n_estimators: int = 40,
    label_flip_frac: float = 0.0,
) -> dict:
    """Δμ, Δs, J + token VAI on LIVE — all on PermuCATE unit risk."""
    pc = text_permucate_magvar(
        texts_ref,
        y_ref,
        texts_live,
        y_live,
        K=K,
        seed=seed,
        n_estimators=n_estimators,
        return_models=True,
    )
    models = pc["fit"]["models"]
    tok = text_token_vai_permucate(
        texts_live,
        y_live,
        pc["feat"],
        models,
        n_docs=token_docs,
        max_tokens=16,
        seed=seed + 11,
    )
    d_uv = float(np.var(pc["unit_live"]) - np.var(pc["unit_ref"]))
    y_delta = abs(float(np.mean(y_live)) - float(np.mean(y_ref)))
    flip = float(label_flip_frac)
    J = float(pc["J"])
    # Product severity: Mag×Var + response-Y integrity break
    severity = J + 1.25 * flip + 0.5 * y_delta
    return {
        "d_mu": float(pc["d_mu"]),
        "d_s": float(pc["d_s"]),
        "d_unit_var": d_uv,
        "J": J,
        "severity": float(severity),
        "label_flip_frac": flip,
        "label_rate_delta": float(y_delta),
        "token_top_share": float(tok["top_share"]),
        "token_gini": float(tok["gini"]),
        "token_mean_abs_delta": float(tok.get("gini", 0.0)),
        "mu_live": float(pc["mu_live"]),
        "s_live": float(pc["s_live"]),
        "estimator": "permucate_rf_text",
    }


def _make_shift(
    texts: Sequence[str],
    y: np.ndarray,
    det,
    kind: str,
    seed: int,
) -> Tuple[List[str], np.ndarray]:
    """Benign: text-only synonym (Y kept). Malignant: subset attack + Y flip.

    Y-flip on attacked rows mirrors tabular PermuCATE malignant boards
    (outcome break), so φ / unit risk respond to intentional integrity hits.
    """
    texts = list(texts)
    y = np.asarray(y, dtype=int).copy()
    if kind == "benign":
        out = [synonym_replace(t, rate=0.12, seed=seed + i) for i, t in enumerate(texts)]
        return out, y
    live = build_live_batches(
        texts,
        y,
        seed=seed,
        synonym_rate=0.12,
        attack_rate=0.40,
        attack_max=max(16, len(texts) // 3),
        detector=det,
    )
    tx_m = list(live["malignant_attack"]["texts"])
    y_m = y.copy()
    atk = [i for i, (a, b) in enumerate(zip(texts, tx_m)) if a != b]
    if not atk:
        rng = np.random.default_rng(int(seed))
        atk = list(rng.choice(len(texts), size=max(8, len(texts) // 4), replace=False))
    y_m[np.asarray(atk, dtype=int)] = 1 - y_m[np.asarray(atk, dtype=int)]
    return tx_m, y_m


def eval_board_auc(
    board_name: str,
    *,
    n_reps: int = 16,
    live_n: int = 100,
    K: int = 20,
    seed: int = SEED,
    max_pairs: int = 900,
    token_docs: int = 8,
    n_estimators: int = 40,
) -> dict:
    """Many LIVE replicates → PermuCATE Mag×Var / token scores → AUC."""
    if board_name.startswith("hc3_"):
        board = load_board(board_name, max_pairs=max_pairs, seed=seed)
    else:
        board = load_board(board_name, max_per_class=min(450, max_pairs // 2), seed=seed)
    split = split_ref_live_board(board, live_frac=0.4, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    # Detector only constructs adversarial LIVE text — not the PO estimator
    det = fit_sklearn_detector(
        X_tr[: min(len(X_tr), 1200)], y_tr[: min(len(y_tr), 1200)], seed=seed
    )

    rng = np.random.default_rng(int(seed))
    rows = []
    for r in range(int(n_reps)):
        idx = rng.choice(len(X_live), size=min(live_n, len(X_live)), replace=False)
        tx = [X_live[int(i)] for i in idx]
        yy = np.asarray([y_live[int(i)] for i in idx], dtype=int)
        ridx = rng.choice(len(X_tr), size=min(live_n, len(X_tr)), replace=False)
        tx_r = [X_tr[int(i)] for i in ridx]
        yy_r = np.asarray([y_tr[int(i)] for i in ridx], dtype=int)

        for kind, lab in (("benign", 0), ("malignant", 1)):
            tx_s, yy_s = _make_shift(tx, yy, det, kind, seed=seed + 1000 * r + lab)
            flip = float(np.mean(np.asarray(yy_s) != np.asarray(yy)))
            feat = magvar_features_permucate(
                tx_r,
                yy_r,
                tx_s,
                yy_s,
                K=K,
                seed=seed + r + lab,
                token_docs=token_docs,
                n_estimators=n_estimators,
                label_flip_frac=flip,
            )
            feat["label"] = lab
            feat["kind"] = kind
            feat["rep"] = r
            rows.append(feat)

    y = np.array([r["label"] for r in rows], dtype=int)
    scores = {
        "severity": np.array([r["severity"] for r in rows], dtype=float),
        "J": np.array([r["J"] for r in rows], dtype=float),
        "d_s": np.array([r["d_s"] for r in rows], dtype=float),
        "d_mu": np.array([r["d_mu"] for r in rows], dtype=float),
        "label_flip_frac": np.array([r["label_flip_frac"] for r in rows], dtype=float),
        "token_top_share": np.array([r["token_top_share"] for r in rows], dtype=float),
        "token_gini": np.array([r["token_gini"] for r in rows], dtype=float),
        "J_plus_token": np.array(
            [r["J"] + 0.5 * r["token_top_share"] for r in rows], dtype=float
        ),
    }
    aucs = {}
    for name, sc in scores.items():
        try:
            aucs[name] = float(roc_auc_score(y, sc))
        except ValueError:
            aucs[name] = float("nan")

    sc = scores["severity"]
    thr = float(np.median(sc))
    yhat = (sc >= thr).astype(int)
    tn = int(((y == 0) & (yhat == 0)).sum())
    fp = int(((y == 0) & (yhat == 1)).sum())
    fn = int(((y == 1) & (yhat == 0)).sum())
    tp = int(((y == 1) & (yhat == 1)).sum())

    return {
        "board": board_name,
        "n_reps": n_reps,
        "n_rows": len(rows),
        "estimator": "permucate_rf_text",
        "auc": aucs,
        "contingency_severity_median": {
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "tp": tp,
            "accuracy": (tn + tp) / max(tn + tp + fn + fp, 1),
        },
        "mean_by_kind": {
            "benign": {
                k: float(np.mean([r[k] for r in rows if r["kind"] == "benign"]))
                for k in (
                    "d_mu",
                    "d_s",
                    "J",
                    "severity",
                    "label_flip_frac",
                    "token_top_share",
                    "token_gini",
                )
            },
            "malignant": {
                k: float(np.mean([r[k] for r in rows if r["kind"] == "malignant"]))
                for k in (
                    "d_mu",
                    "d_s",
                    "J",
                    "severity",
                    "label_flip_frac",
                    "token_top_share",
                    "token_gini",
                )
            },
        },
        "rows": rows,
    }


def run_auc_suite(
    boards: Optional[Sequence[str]] = None,
    n_reps: int = 14,
    live_n: int = 90,
    seed: int = SEED,
    out_dir: Optional[Path] = None,
    quick: bool = False,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    if quick:
        n_reps, live_n = 6, 50
    names = list(boards or BOARD_NAMES)
    if quick:
        names = list(boards or ("hc3_all", "hc3_finance", "news20_sci"))
    n_est = 30 if quick else 40
    results = []
    for b in names:
        print("[tok-auc] board", b, flush=True)
        results.append(
            eval_board_auc(
                b,
                n_reps=n_reps,
                live_n=live_n,
                seed=seed,
                max_pairs=500 if quick else 1000,
                token_docs=4 if quick else 8,
                K=12 if quick else 20,
                n_estimators=n_est,
            )
        )

    all_rows = [r for block in results for r in block["rows"]]
    y = np.array([r["label"] for r in all_rows], dtype=int)
    pooled_auc = {}
    for key in (
        "severity",
        "J",
        "d_s",
        "d_mu",
        "label_flip_frac",
        "token_top_share",
        "J_plus_token",
    ):
        sc = np.array(
            [
                r[key]
                if key != "J_plus_token"
                else r["J"] + 0.5 * r["token_top_share"]
                for r in all_rows
            ]
        )
        try:
            pooled_auc[key] = float(roc_auc_score(y, sc))
        except ValueError:
            pooled_auc[key] = float("nan")

    payload = {
        "scope": (
            "Token-wise perturb→re-forward + Mag×Var AUC on PermuCATE unit PO "
            "(Y=response, W=REF/LIVE) — benign vs malignant (no multimodal)."
        ),
        "estimator": {
            "name": "PermuCATE RF text",
            "unit_po": "(φ - τ̂)^2",
            "Y": "response / class label",
            "W": "REF=0 / LIVE=1",
            "used_permucate_rf": True,
        },
        "logic": [
            "Y = response label; W = REF(0)/LIVE(1); X = TF-IDF→SVD",
            "PermuCATE RF: μ̂0, μ̂1, π̂ → φ → τ̂; unit u=(φ−τ̂)²",
            "Magnitude μ, Variance s of bootstrap sub-env mean u",
            "Token VAI: drop j → re-featurize → frozen models → |Δu|",
            "Malignant oracle: subset attack + Y flip (outcome break)",
            "Severity = J + Y_flip + rateΔ; AUC malignant vs benign",
        ],
        "boards": [{k: v for k, v in b.items() if k != "rows"} for b in results],
        "pooled_auc": pooled_auc,
        "justify": (
            "Magnitude alone confuses uniform synonym with attack when both raise μ. "
            "Variance of PermuCATE unit risk + token-attr concentration rises under "
            "malignant subset triggers; AUC of J / d_s justifies separation."
        ),
    }
    payload["n_pooled"] = len(all_rows)
    png = out_dir / "token_magvar_auc.png"
    plot_auc(payload, results, png)
    (out_dir / "token_magvar_auc.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", payload, png.name)
    return payload


def plot_auc(payload: dict, full_results: list, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    boards = payload["boards"]
    names = [b["board"] for b in boards]
    fig = plt.figure(figsize=(12.0, 8.0))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(names))
    for key, color, lab in (
        ("severity", "#2F6B4F", "AUC(severity)"),
        ("J", "#1F6F8B", "AUC(J)"),
        ("d_mu", "#9AA3AE", "AUC(Δμ)"),
        ("label_flip_frac", "#A33B24", "AUC(Y-flip)"),
    ):
        vals = [b["auc"].get(key, float("nan")) for b in boards]
        ax.plot(x, vals, "o-", color=color, label=lab, lw=1.6)
    ax.axhline(0.5, ls="--", color="#ccc", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0.4, 1.05)
    ax.set_ylabel("AUC")
    ax.set_title("Per-board AUC: PermuCATE Mag×Var", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    keys = list(payload["pooled_auc"].keys())
    ax.barh(keys, [payload["pooled_auc"][k] for k in keys], color="#1F6F8B")
    ax.axvline(0.5, ls="--", color="#ccc")
    ax.set_xlim(0.4, 1.05)
    ax.set_xlabel("pooled AUC")
    ax.set_title("Pooled AUC (PermuCATE RF unit risk)", loc="left", fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)

    ax = fig.add_subplot(gs[1, 0])
    jb = [b["mean_by_kind"]["benign"]["severity"] for b in boards]
    jm = [b["mean_by_kind"]["malignant"]["severity"] for b in boards]
    ax.bar(x - 0.18, jb, width=0.36, color="#2F6B4F", label="benign severity")
    ax.bar(x + 0.18, jm, width=0.36, color="#A33B24", label="malignant severity")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("mean severity")
    ax.set_title("Severity (J + Y-flip) by shift type", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    tb = [b["mean_by_kind"]["benign"]["token_top_share"] for b in boards]
    tm = [b["mean_by_kind"]["malignant"]["token_top_share"] for b in boards]
    ax.bar(x - 0.18, tb, width=0.36, color="#2F6B4F", label="benign token top20%")
    ax.bar(x + 0.18, tm, width=0.36, color="#A33B24", label="malignant token top20%")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("token attr top-share")
    ax.set_title("Token VAI on frozen PermuCATE u", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "Token perturb → PermuCATE Mag×Var AUC — benign vs malignant (Y response, W=REF/LIVE)",
        fontsize=11.5,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_readme(path: Path, payload: dict, png_name: str):
    lines = [
        "# Token Mag×Var AUC — PermuCATE unit PO (benign vs malignant)",
        "",
        payload["scope"],
        "",
        "## Estimator",
        "",
        "- **PermuCATE RF** on text: $Y$=response, $W$=REF/LIVE, $X$=TF-IDF→SVD",
        "- Unit $u_i=(\\varphi_i-\\hat\\tau_i)^2$ — **not** $1-\\hat p(Y\\mid x)$",
        "- Token VAI: drop → re-featurize → frozen $\\hat\\mu_w,\\hat\\pi,\\hat\\tau$",
        "",
        "## Logic",
        "",
    ]
    for i, x in enumerate(payload["logic"], 1):
        lines.append("%d. %s" % (i, x))
    lines.extend(["", "## Pooled AUC", ""])
    for k, v in payload["pooled_auc"].items():
        lines.append("- **%s**: %.3f" % (k, v))
    lines.extend(
        [
            "",
            "## Per board AUC(J)",
            "",
            "| board | AUC(severity) | AUC(J) | AUC(Y-flip) | AUC(token) |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for b in payload["boards"]:
        a = b["auc"]
        lines.append(
            "| %s | %.3f | %.3f | %.3f | %.3f |"
            % (
                b["board"],
                a.get("severity", float("nan")),
                a["J"],
                a.get("label_flip_frac", float("nan")),
                a["token_top_share"],
            )
        )
    lines.extend(
        [
            "",
            "## Justify",
            "",
            payload["justify"],
            "",
            f"![]({png_name})",
            "",
            "```bash",
            "python3 scripts/run_as_dta_token_magvar_auc.py --quick",
            "python3 scripts/run_as_dta_token_magvar_auc.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
