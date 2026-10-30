"""Real-world applications to justify Var(u) product logic.

Boards (binary text, same Y-flip口径 as Mag×Var delivery):
  hc3_finance, hc3_medicine   — AIGC integrity
  news20_sci, news20_comp     — topic monitor
  ag_news                     — news topic (World vs Sci/Tech)
  imdb                        — sentiment

Per board: benign synonym (Y intact) vs malignant attack + Y-flip on
attacked rows. Product must call conformal on benign and VAI_freeze on
malignant under CF-debiased Var(u) + Mag×Var + Y-flip.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from po_risk_var._path import ensure_cfperm_src
from po_risk_var.estimators import monitor_batch

ensure_cfperm_src()

from as_dta_aigc import build_live_batches, fit_sklearn_detector, synonym_replace  # noqa: E402
from as_dta_permucate_magvar import featurize_texts  # noqa: E402
from as_dta_shift_data import load_board, split_ref_live_board  # noqa: E402

REAL_BOARDS = (
    "hc3_finance",
    "hc3_medicine",
    "news20_sci",
    "news20_comp",
    "ag_news",
    "imdb",
)


def _malignant_yflip(texts, y, det, seed: int):
    y = np.asarray(y, dtype=int).copy()
    live = build_live_batches(
        texts,
        y,
        seed=seed,
        synonym_rate=0.10,
        attack_rate=0.45,
        attack_max=max(10, len(texts) // 3),
        detector=det,
    )
    tx = list(live["malignant_attack"]["texts"])
    atk = [i for i, (a, b) in enumerate(zip(texts, tx)) if a != b]
    if not atk:
        rng = np.random.default_rng(int(seed))
        atk = list(rng.choice(len(texts), size=max(6, len(texts) // 4), replace=False))
    y2 = y.copy()
    y2[np.asarray(atk, dtype=int)] = 1 - y2[np.asarray(atk, dtype=int)]
    flip = float(np.mean(y2 != y))
    return tx, y2, flip


def eval_real_board(
    board_name: str,
    *,
    live_n: int = 56,
    seed: int = 2026,
    quick: bool = True,
    method: str = "cf",
) -> dict:
    board = load_board(
        board_name,
        max_pairs=320 if quick else 700,
        max_per_class=160 if quick else 350,
        seed=seed,
    )
    split = split_ref_live_board(board, live_frac=0.4, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    det = fit_sklearn_detector(
        X_tr[: min(800, len(X_tr))], y_tr[: min(800, len(y_tr))], seed=seed
    )
    rng = np.random.default_rng(int(seed))
    idx = rng.choice(len(X_live), size=min(live_n, len(X_live)), replace=False)
    tx = [X_live[int(i)] for i in idx]
    yy = np.asarray([y_live[int(i)] for i in idx], dtype=int)
    ridx = rng.choice(len(X_tr), size=min(live_n, len(X_tr)), replace=False)
    tx_r = [X_tr[int(i)] for i in ridx]
    yy_r = np.asarray([y_tr[int(i)] for i in ridx], dtype=int)

    rows = {}
    for kind in ("benign", "malignant"):
        if kind == "benign":
            tx_s = [synonym_replace(t, rate=0.14, seed=seed + i) for i, t in enumerate(tx)]
            yy_s = yy.copy()
            flip = 0.0
            oracle = "benign"
        else:
            tx_s, yy_s, flip = _malignant_yflip(tx, yy, det, seed=seed + 17)
            oracle = "malignant"

        feat = featurize_texts(
            tx_r,
            tx_s,
            max_features=350 if quick else 500,
            n_components=48 if quick else 64,
            seed=seed,
        )
        mon = monitor_batch(
            feat["X_ref"],
            yy_r,
            feat["X_live"],
            yy_s,
            method=method,
            y_flip=flip,
            quick=quick,
            seed=seed + (0 if kind == "benign" else 5),
        )
        pred = mon["decision"]["call"]
        ok = (oracle == "malignant" and pred == "malignant") or (
            oracle == "benign" and pred in ("quiet", "benign")
        )
        rows[kind] = {
            "oracle": oracle,
            "pred": pred,
            "actuator": mon["decision"]["actuator"],
            "correct": bool(ok),
            "y_flip": flip,
            "d_mu": mon["magvar"]["d_mu"],
            "d_s": mon["magvar"]["d_s"],
            "J": mon["magvar"]["J"],
            "V_debiased": mon["var_u"]["mean_V_debiased_live"],
            "V_naive": mon["var_u"]["mean_V_naive_live"],
            "bias_ratio": mon["var_u"]["bias_ratio_live"],
            "reason": mon["decision"]["reason"],
        }

    return {
        "board": board_name,
        "benign": rows["benign"],
        "malignant": rows["malignant"],
        "pair_ok": bool(rows["benign"]["correct"] and rows["malignant"]["correct"]),
        "V_gap": float(
            rows["malignant"]["V_debiased"] - rows["benign"]["V_debiased"]
        ),
    }


def run_realworld_apps(
    boards: Optional[Sequence[str]] = None,
    *,
    seed: int = 2026,
    quick: bool = True,
    method: str = "cf",
) -> dict:
    names = list(boards or REAL_BOARDS)
    if quick:
        # keep morning-check runtime tight but cover 4 families
        names = list(
            boards
            or ("hc3_finance", "news20_sci", "ag_news", "imdb")
        )
    rows: List[dict] = []
    for b in names:
        print("[real] ", b, flush=True)
        rows.append(eval_real_board(b, seed=seed, quick=quick, method=method))
    n_ok = sum(1 for r in rows if r["pair_ok"])
    return {
        "scope": (
            "Real-world boards — same Mag×Var + Y-flip口径 + CF-debiased Var(u) "
            "product triage (AIGC / topic / news / sentiment)."
        ),
        "method": method,
        "boards": rows,
        "n_boards": len(rows),
        "n_pair_ok": n_ok,
        "hit_rate": float(n_ok) / max(len(rows), 1),
    }
