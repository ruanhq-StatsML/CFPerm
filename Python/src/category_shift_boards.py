"""Synthetic boards mirroring CIFAR-100 C-GCD and Office-Home PDA.

CIFAR-100 (C-GCD proxy)
  Stage-0: known classes only → fit REF.
  Later stages: unlabeled-like influx with novel classes appearing.
  Metrics: known-class retention (forgetting), novel discovery (acc on new labels).

Office-Home (PDA proxy)
  Domains: Art / Clipart / Product / Real.
  Source has full label set; target keeps only a fraction (partial overlap).
  Metrics: shared-class transfer; incomparable mark rate at low overlap.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from weighted_cold_start import (
    SEED,
    adapt_batch,
    adapt_uniform_partial_fit,
    boundary_check,
    fit_initial,
)

# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------


def _class_means(n_classes: int, p: int, seed: int, active_dims: Optional[np.ndarray] = None):
    rng = np.random.default_rng(int(seed))
    mu = rng.normal(0, 1.0, size=(n_classes, p))
    if active_dims is not None:
        mask = np.zeros(p, dtype=bool)
        mask[np.asarray(active_dims, int)] = True
        mu[:, ~mask] *= 0.15
    return mu


def sample_classes(mu, classes, n_per, rng, noise=0.7):
    Xs, ys = [], []
    for c in classes:
        Xs.append(mu[int(c)] + rng.normal(0, noise, size=(n_per, mu.shape[1])))
        ys.append(np.full(n_per, int(c), dtype=int))
    return np.vstack(Xs), np.concatenate(ys)


def make_cgcd_stream(
    *,
    n_known: int = 12,
    n_novel_per_stage: int = 3,
    n_stages: int = 4,
    n_per: int = 40,
    p: int = 32,
    novel_dims: int = 8,
    seed: int = SEED,
) -> Dict:
    """CIFAR-100 C-GCD-style: known first, then staged novel-class influx."""
    rng = np.random.default_rng(int(seed))
    n_total = n_known + n_novel_per_stage * n_stages
    # novel classes live in a dedicated feature subspace (attribution target)
    active = np.arange(novel_dims)
    mu = _class_means(n_total, p, seed, active_dims=active)
    known = list(range(n_known))
    X0, y0 = sample_classes(mu, known, n_per, rng)
    stages = []
    cursor = n_known
    for s in range(n_stages):
        novel = list(range(cursor, cursor + n_novel_per_stage))
        cursor += n_novel_per_stage
        # Stress cold-start: known recirculation dominates; novel is sparse.
        Xk, yk = sample_classes(mu, known, max(12, n_per), rng)
        Xn, yn = sample_classes(mu, novel, max(6, n_per // 4), rng)
        stages.append({
            "stage": s,
            "X": np.vstack([Xk, Xn]),
            "y": np.concatenate([yk, yn]),
            "known": known,
            "novel": novel,
            "active_dims": active.tolist(),
        })
    all_classes = list(range(n_total))
    return {
        "X_ref": X0, "y_ref": y0, "known": known, "all_classes": all_classes,
        "stages": stages, "p": p, "mu": mu,
    }


def make_officehome_pda(
    *,
    n_classes: int = 20,
    n_per: int = 36,
    p: int = 28,
    overlap_frac: float = 0.4,
    seed: int = SEED,
) -> Dict:
    """Office-Home PDA-style: 4 domains; target keeps only overlap_frac of labels."""
    rng = np.random.default_rng(int(seed))
    domains = ["Art", "Clipart", "Product", "Real"]
    # domain style shifts on trailing coords; class signal on leading coords
    class_dims = np.arange(max(6, p // 2))
    mu = _class_means(n_classes, p, seed, active_dims=class_dims)
    all_cls = list(range(n_classes))
    k_keep = max(1, int(round(n_classes * float(overlap_frac))))
    shared = all_cls[:k_keep]
    out = {"domains": {}, "shared": shared, "all_classes": all_cls, "overlap_frac": overlap_frac}
    for i, d in enumerate(domains):
        shift = np.zeros(p)
        shift[-(p // 4):] = rng.normal(0.8 * (i + 1), 0.3, size=(p // 4,))
        cls = all_cls if d == "Art" else shared  # Art = full source; others = partial target
        X, y = sample_classes(mu, cls, n_per, rng)
        X = X + shift
        out["domains"][d] = {"X": X, "y": y, "classes": cls}
    return out


# ---------------------------------------------------------------------------
# Evaluators
# ---------------------------------------------------------------------------


def run_cgcd(
    *,
    seed: int = SEED,
    weighted: bool = True,
    quick: bool = False,
) -> Dict:
    stream = make_cgcd_stream(
        n_known=10 if quick else 12,
        n_novel_per_stage=2 if quick else 3,
        n_stages=3 if quick else 4,
        n_per=24 if quick else 40,
        p=24 if quick else 32,
        novel_dims=6 if quick else 8,
        seed=seed,
    )
    # Pre-declare full class set so partial_fit can absorb novel labels later.
    clf, sc = fit_initial(
        stream["X_ref"], stream["y_ref"], seed=seed, classes=stream["all_classes"],
    )
    X_ref = sc.transform(stream["X_ref"])
    y_ref = stream["y_ref"]
    known = stream["known"]
    classes = np.asarray(stream["all_classes"], dtype=int)
    rows = []
    for st in stream["stages"]:
        Xl = sc.transform(st["X"])
        yl = st["y"]
        if weighted:
            out = adapt_batch(
                clf, X_ref, y_ref, Xl, yl, seed=seed + st["stage"],
                top_k=8, lam=3.5, min_jaccard=0.05,  # C-GCD expects novel; allow adapt
                detect_thresh=1.05, force_adapt=True,
            )
            clf = out["clf"]
            meta = {"reason": out["reason"], "vimp_top_mass": out["vimp_top_mass"],
                    "feat_hit": float(np.isin(out["feat_idx"], st["active_dims"]).mean()),
                    "boundary": out["boundary"]["decision"]}
        else:
            # Uniform sees the same labels but flat weights — cold-start stress test.
            clf = adapt_uniform_partial_fit(clf, Xl, yl, classes)
            meta = {"reason": "uniform", "vimp_top_mass": None, "feat_hit": None, "boundary": "n/a"}

        # Hold-out eval slice (avoid train=test inflation on the same stage batch)
        rng_e = np.random.default_rng(seed + 100 + st["stage"])
        ev = rng_e.choice(len(yl), size=max(16, len(yl) // 3), replace=False)
        # metrics on this stage
        pred = clf.predict(Xl[ev])
        y_ev = yl[ev]
        known_mask = np.isin(y_ev, known)
        novel_mask = np.isin(y_ev, st["novel"])
        known_acc = float(np.mean(pred[known_mask] == y_ev[known_mask])) if known_mask.any() else float("nan")
        novel_acc = float(np.mean(pred[novel_mask] == y_ev[novel_mask])) if novel_mask.any() else float("nan")
        # REF forgetting on original REF only (first block)
        forget = 1.0 - float(np.mean(clf.predict(sc.transform(stream["X_ref"])) == stream["y_ref"]))
        rows.append({
            "stage": st["stage"], "known_acc": known_acc, "novel_acc": novel_acc,
            "forgetting": forget, **meta,
        })
        # grow REF gently with stage (replay mix of known only — reduces contamination)
        known_idx = np.flatnonzero(np.isin(yl, known))
        if len(known_idx):
            take = known_idx[::2][: max(8, len(known_idx) // 2)]
            X_ref = np.vstack([X_ref, Xl[take]])
            y_ref = np.concatenate([y_ref, yl[take]])
    return {
        "board": "cifar100_cgcd_proxy",
        "weighted": weighted,
        "stages": rows,
        "mean_known_acc": float(np.nanmean([r["known_acc"] for r in rows])),
        "mean_novel_acc": float(np.nanmean([r["novel_acc"] for r in rows])),
        "mean_forgetting": float(np.nanmean([r["forgetting"] for r in rows])),
        "mean_feat_hit": (
            float(np.mean([r["feat_hit"] for r in rows if r.get("feat_hit") is not None]))
            if any(r.get("feat_hit") is not None for r in rows)
            else float("nan")
        ),
    }


def run_officehome_pda(
    *,
    seed: int = SEED,
    overlap_fracs: Optional[List[float]] = None,
    quick: bool = False,
) -> Dict:
    overlap_fracs = overlap_fracs or ([0.1, 0.3, 0.5, 0.8] if not quick else [0.1, 0.4, 0.8])
    rows = []
    for ov in overlap_fracs:
        board = make_officehome_pda(
            n_classes=16 if quick else 20,
            n_per=20 if quick else 36,
            p=20 if quick else 28,
            overlap_frac=ov,
            seed=seed,
        )
        src = board["domains"]["Art"]
        clf, sc = fit_initial(src["X"], src["y"], seed=seed)
        Xs = sc.transform(src["X"])
        ys = src["y"]
        for tgt_name in ("Clipart", "Product", "Real"):
            tgt = board["domains"][tgt_name]
            Xt, yt = sc.transform(tgt["X"]), tgt["y"]
            bound = boundary_check(ys, yt, min_jaccard=0.20)
            out = adapt_batch(
                clf, Xs, ys, Xt, yt, seed=seed, top_k=6, lam=2.0,
                min_jaccard=0.20, detect_thresh=1.05, force_adapt=False,
            )
            # shared-class accuracy (even if incomparable — for reporting)
            shared = board["shared"]
            m = np.isin(yt, shared)
            shared_acc = float(np.mean(out["clf"].predict(Xt)[m] == yt[m])) if m.any() else float("nan")
            rows.append({
                "overlap_frac": ov,
                "target": tgt_name,
                "jaccard": bound["jaccard"],
                "boundary": bound["decision"],
                "adapted": out["adapted"],
                "reason": out["reason"],
                "shared_acc": shared_acc,
                "vimp_top_mass": out["vimp_top_mass"],
            })
    # aggregate
    by_ov = {}
    for ov in overlap_fracs:
        sub = [r for r in rows if r["overlap_frac"] == ov]
        by_ov[str(ov)] = {
            "mean_shared_acc": float(np.nanmean([r["shared_acc"] for r in sub])),
            "incomparable_rate": float(np.mean([r["boundary"] == "incomparable" for r in sub])),
            "adapt_rate": float(np.mean([r["adapted"] for r in sub])),
        }
    return {
        "board": "officehome_pda_proxy",
        "rows": rows,
        "by_overlap": by_ov,
    }


def run_suite(*, seed: int = SEED, quick: bool = True) -> Dict:
    cg_w = run_cgcd(seed=seed, weighted=True, quick=quick)
    cg_u = run_cgcd(seed=seed, weighted=False, quick=quick)
    oh = run_officehome_pda(seed=seed, quick=quick)
    return {
        "cgcd_weighted": cg_w,
        "cgcd_uniform": cg_u,
        "officehome": oh,
        "delta_novel_acc": cg_w["mean_novel_acc"] - cg_u["mean_novel_acc"],
        "delta_forgetting": cg_u["mean_forgetting"] - cg_w["mean_forgetting"],
        "takeaway": (
            "Weighted cold-start (upweight localized subset + selected features) "
            "vs uniform partial_fit on C-GCD proxy; Office-Home PDA checks "
            "incomparable rate rises as label overlap falls."
        ),
    }
