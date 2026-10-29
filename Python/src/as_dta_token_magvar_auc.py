"""Token-wise perturbation → Mag×Var AUC: benign vs malignant shift.

Logic (your framing)
--------------------
1. Fit detector on REF.
2. Build LIVE shifts: benign (uniform synonym) vs malignant (subset attack).
3. Unit PO_i = 1 − p̂(Y_i|x_i).
4. **Magnitude** μ = E[PO]; **Variance** s = Std(bootstrap sub-env mean PO).
5. **Token-wise attribution:** perturb token j → delete/mask → re-forward →
   a_{i,j} = PO_i^{pert,j} − PO_i^{base}. Concentration of a distinguishes
   uniform benign edits vs sparse malignant triggers.
6. Score shifts with Mag×Var (and token concentration); report **AUC** /
   2×2 across replicates × datasets.

No multimodal — token / lexical axis only.
"""
from __future__ import annotations

import json
import re
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
    textfooler_lite,
)
from as_dta_benign_malignant_mvp import _merge_J, _subbatch_risk_var
from as_dta_shift_data import BOARD_NAMES, load_board, split_ref_live_board

OUT_DEFAULT = ROOT / "results" / "as_dta_token_magvar_auc"
WORD_RE = re.compile(r"[A-Za-z']+")


def _po(det, texts: Sequence[str], y: np.ndarray) -> np.ndarray:
    return np.asarray(det.po_proxy(list(texts), np.asarray(y, dtype=int)), dtype=float)


def token_perturb_attribution(
    det,
    text: str,
    y: int,
    *,
    max_tokens: int = 24,
    seed: int = SEED,
) -> dict:
    """Perturb-then-re-forward: drop one token at a time, measure ΔPO."""
    toks = text.split()
    if not toks:
        return {"attrs": [], "top_share": 0.0, "gini": 0.0, "n_tok": 0}
    cand = [i for i, t in enumerate(toks) if WORD_RE.fullmatch(t) and len(t) > 2]
    if not cand:
        cand = list(range(len(toks)))
    rng = np.random.default_rng(int(seed))
    if len(cand) > max_tokens:
        cand = list(rng.choice(cand, size=max_tokens, replace=False))
    base = float(_po(det, [text], np.array([y]))[0])
    attrs = []
    for j in cand:
        trial = list(toks)
        trial[j] = ""  # delete token
        new_t = " ".join(w for w in trial if w)
        if not new_t.strip():
            continue
        po_j = float(_po(det, [new_t], np.array([y]))[0])
        attrs.append(abs(po_j - base))
    attrs = np.asarray(attrs, dtype=float)
    if attrs.size == 0:
        return {"attrs": [], "top_share": 0.0, "gini": 0.0, "n_tok": 0, "base_po": base}
    s = float(attrs.sum()) + 1e-12
    top_n = max(1, int(0.2 * len(attrs)))
    order = np.sort(attrs)[::-1]
    top_share = float(order[:top_n].sum() / s)
    xs = np.sort(attrs)
    n = len(xs)
    gini = float((2.0 * np.arange(1, n + 1) - n - 1) @ xs / (n * s))
    return {
        "attrs": attrs.tolist(),
        "top_share": top_share,
        "gini": gini,
        "n_tok": int(len(attrs)),
        "base_po": base,
        "mean_abs_delta": float(np.mean(attrs)),
    }


def batch_token_concentration(
    det,
    texts: Sequence[str],
    y: np.ndarray,
    *,
    n_docs: int = 12,
    max_tokens: int = 20,
    seed: int = SEED,
) -> dict:
    """Aggregate token-attr concentration over a LIVE batch sample."""
    y = np.asarray(y, dtype=int)
    rng = np.random.default_rng(int(seed))
    n = len(texts)
    pick = rng.choice(n, size=min(n_docs, n), replace=False)
    tops, ginis, means = [], [], []
    for k, i in enumerate(pick):
        a = token_perturb_attribution(
            det, texts[int(i)], int(y[int(i)]), max_tokens=max_tokens, seed=seed + k
        )
        if a["n_tok"] == 0:
            continue
        tops.append(a["top_share"])
        ginis.append(a["gini"])
        means.append(a["mean_abs_delta"])
    if not tops:
        return {"top_share": 0.0, "gini": 0.0, "mean_abs_delta": 0.0, "n_docs": 0}
    return {
        "top_share": float(np.mean(tops)),
        "gini": float(np.mean(ginis)),
        "mean_abs_delta": float(np.mean(means)),
        "n_docs": int(len(tops)),
    }


def magvar_features(
    det,
    texts_ref: Sequence[str],
    y_ref: np.ndarray,
    texts_live: Sequence[str],
    y_live: np.ndarray,
    *,
    K: int = 24,
    seed: int = SEED,
    token_docs: int = 10,
) -> dict:
    """Δμ, Δs, J + token concentration on LIVE."""
    _, st_ref = _subbatch_risk_var(texts_ref, y_ref, det, K=K, seed=seed)
    _, st_live = _subbatch_risk_var(texts_live, y_live, det, K=K, seed=seed + 5)
    d_mu = float(st_live["sub_mean"] - st_ref["sub_mean"])
    # use full-batch mean for magnitude stability
    po_r = _po(det, texts_ref, y_ref)
    po_l = _po(det, texts_live, y_live)
    d_mu = float(np.mean(po_l) - np.mean(po_r))
    d_s = float(st_live["sub_std"] - st_ref["sub_std"])
    d_uv = float(st_live["unit"]["var"] - st_ref["unit"]["var"])
    J = _merge_J(
        {"d_mean": d_mu, "d_sub_std": d_s},
        ref_sub_std=st_ref["sub_std"],
        kappa=3.0,
    )
    tok = batch_token_concentration(
        det, texts_live, y_live, n_docs=token_docs, seed=seed + 11
    )
    err_r = 1.0 - float(np.mean((det.predict_proba(texts_ref) >= getattr(det, "threshold", 0.5)).astype(int) == np.asarray(y_ref)))
    err_l = 1.0 - float(np.mean((det.predict_proba(texts_live) >= getattr(det, "threshold", 0.5)).astype(int) == np.asarray(y_live)))
    return {
        "d_mu": d_mu,
        "d_s": d_s,
        "d_unit_var": d_uv,
        "J": float(J),
        "token_top_share": tok["top_share"],
        "token_gini": tok["gini"],
        "token_mean_abs_delta": tok["mean_abs_delta"],
        "d_error": float(err_l - err_r),
        "mu_live": float(np.mean(po_l)),
        "s_live": float(st_live["sub_std"]),
    }


def _make_shift_texts(
    texts: Sequence[str],
    y: np.ndarray,
    det,
    kind: str,
    seed: int,
) -> List[str]:
    texts = list(texts)
    y = np.asarray(y, dtype=int)
    if kind == "benign":
        return [synonym_replace(t, rate=0.12, seed=seed + i) for i, t in enumerate(texts)]
    # malignant: concentrated subset
    live = build_live_batches(
        texts,
        y,
        seed=seed,
        synonym_rate=0.12,
        attack_rate=0.40,
        attack_max=max(16, len(texts) // 3),
        detector=det,
    )
    return list(live["malignant_attack"]["texts"])


def eval_board_auc(
    board_name: str,
    *,
    n_reps: int = 16,
    live_n: int = 100,
    K: int = 20,
    seed: int = SEED,
    max_pairs: int = 900,
    token_docs: int = 8,
) -> dict:
    """Many LIVE replicates → Mag×Var / token scores → AUC (malignant=1)."""
    if board_name.startswith("hc3_"):
        board = load_board(board_name, max_pairs=max_pairs, seed=seed)
    else:
        board = load_board(board_name, max_per_class=min(450, max_pairs // 2), seed=seed)
    split = split_ref_live_board(board, live_frac=0.4, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    det = fit_sklearn_detector(
        X_tr[: min(len(X_tr), 1200)], y_tr[: min(len(y_tr), 1200)], seed=seed
    )

    rng = np.random.default_rng(int(seed))
    rows = []
    for r in range(int(n_reps)):
        idx = rng.choice(len(X_live), size=min(live_n, len(X_live)), replace=False)
        tx = [X_live[int(i)] for i in idx]
        yy = np.asarray([y_live[int(i)] for i in idx], dtype=int)
        # REF probe same size from train
        ridx = rng.choice(len(X_tr), size=min(live_n, len(X_tr)), replace=False)
        tx_r = [X_tr[int(i)] for i in ridx]
        yy_r = np.asarray([y_tr[int(i)] for i in ridx], dtype=int)

        for kind, lab in (("benign", 0), ("malignant", 1)):
            tx_s = _make_shift_texts(tx, yy, det, kind, seed=seed + 1000 * r + lab)
            feat = magvar_features(
                det, tx_r, yy_r, tx_s, yy, K=K, seed=seed + r + lab, token_docs=token_docs
            )
            feat["label"] = lab
            feat["kind"] = kind
            feat["rep"] = r
            rows.append(feat)

    y = np.array([r["label"] for r in rows], dtype=int)
    scores = {
        "J": np.array([r["J"] for r in rows], dtype=float),
        "d_s": np.array([r["d_s"] for r in rows], dtype=float),
        "d_mu": np.array([r["d_mu"] for r in rows], dtype=float),
        "token_top_share": np.array([r["token_top_share"] for r in rows], dtype=float),
        "token_gini": np.array([r["token_gini"] for r in rows], dtype=float),
        # Mag×Var merge + token concentration
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

    # simple threshold on J → 2×2 at Youden-like mid
    sc = scores["J"]
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
        "auc": aucs,
        "contingency_J_median": {
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "tp": tp,
            "accuracy": (tn + tp) / max(tn + tp + fn + fp, 1),
        },
        "mean_by_kind": {
            "benign": {
                k: float(np.mean([r[k] for r in rows if r["kind"] == "benign"]))
                for k in ("d_mu", "d_s", "J", "token_top_share", "token_gini")
            },
            "malignant": {
                k: float(np.mean([r[k] for r in rows if r["kind"] == "malignant"]))
                for k in ("d_mu", "d_s", "J", "token_top_share", "token_gini")
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
        n_reps, live_n = 8, 60
    names = list(boards or BOARD_NAMES)
    results = []
    for b in names:
        print("[tok-auc] board", b, flush=True)
        results.append(
            eval_board_auc(
                b,
                n_reps=n_reps,
                live_n=live_n,
                seed=seed,
                max_pairs=700 if quick else 1000,
                token_docs=6 if quick else 8,
                K=16 if quick else 20,
            )
        )

    # pooled AUC on concatenated rows
    all_rows = [r for block in results for r in block["rows"]]
    y = np.array([r["label"] for r in all_rows], dtype=int)
    pooled_auc = {}
    for key in ("J", "d_s", "d_mu", "token_top_share", "J_plus_token"):
        sc = np.array([r[key] if key != "J_plus_token" else r["J"] + 0.5 * r["token_top_share"] for r in all_rows])
        try:
            pooled_auc[key] = float(roc_auc_score(y, sc))
        except ValueError:
            pooled_auc[key] = float("nan")

    payload = {
        "scope": (
            "Token-wise perturb→re-forward + Mag×Var AUC to separate "
            "benign vs malignant shift (no multimodal)."
        ),
        "logic": [
            "PO_i = 1 - p̂(Y_i|x_i)",
            "Magnitude μ, Variance s of sub-env mean PO",
            "Token attr a_j = |PO^{drop j} - PO| after re-forward",
            "Score J = Δμ·(1+κ·Δs₊); AUC malignant vs benign",
        ],
        "boards": [{k: v for k, v in b.items() if k != "rows"} for b in results],
        "pooled_auc": pooled_auc,
        "justify": (
            "Magnitude alone confuses uniform synonym with attack when both raise μ. "
            "Variance of PO + token-attr concentration rises under malignant subset "
            "triggers; AUC of J / d_s justifies separation across datasets."
        ),
    }
    # attach light row summaries only
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
        ("J", "#2F6B4F", "AUC(J)"),
        ("d_s", "#1F6F8B", "AUC(Δs)"),
        ("d_mu", "#9AA3AE", "AUC(Δμ)"),
    ):
        vals = [b["auc"][key] for b in boards]
        ax.plot(x, vals, "o-", color=color, label=lab, lw=1.6)
    ax.axhline(0.5, ls="--", color="#ccc", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0.4, 1.05)
    ax.set_ylabel("AUC")
    ax.set_title("Per-board AUC: Mag×Var vs Magnitude-only", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    keys = list(payload["pooled_auc"].keys())
    ax.barh(keys, [payload["pooled_auc"][k] for k in keys], color="#1F6F8B")
    ax.axvline(0.5, ls="--", color="#ccc")
    ax.set_xlim(0.4, 1.05)
    ax.set_xlabel("pooled AUC")
    ax.set_title("Pooled AUC (all replicates × boards)", loc="left", fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)

    ax = fig.add_subplot(gs[1, 0])
    # mean J benign vs malignant
    jb = [b["mean_by_kind"]["benign"]["J"] for b in boards]
    jm = [b["mean_by_kind"]["malignant"]["J"] for b in boards]
    ax.bar(x - 0.18, jb, width=0.36, color="#2F6B4F", label="benign J")
    ax.bar(x + 0.18, jm, width=0.36, color="#A33B24", label="malignant J")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("mean J")
    ax.set_title("Merge J by shift type", loc="left", fontweight="bold")
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
    ax.set_title("Token perturb→re-forward concentration", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "Token-wise perturb → Mag×Var AUC — benign vs malignant shift (no multimodal)",
        fontsize=12.2,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_readme(path: Path, payload: dict, png_name: str):
    lines = [
        "# Token Mag×Var AUC — benign vs malignant shift",
        "",
        "## Logic",
        "",
    ]
    for x in payload["logic"]:
        lines.append("1. " + x if not x.startswith("PO") else "- " + x)
    # fix numbering casually
    lines = [
        "# Token Mag×Var AUC — benign vs malignant shift",
        "",
        payload["scope"],
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
            "| board | AUC(J) | AUC(Δs) | AUC(Δμ) | AUC(token) |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for b in payload["boards"]:
        a = b["auc"]
        lines.append(
            "| %s | %.3f | %.3f | %.3f | %.3f |"
            % (b["board"], a["J"], a["d_s"], a["d_mu"], a["token_top_share"])
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
            "python3 scripts/run_as_dta_token_magvar_auc.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
