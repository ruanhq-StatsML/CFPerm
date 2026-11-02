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
# Evaluators — four axes:
#   1) novel_acc   2) known retention / forgetting
#   3) attribution (feat_hit, vimp_top_mass)
#   4) cold-start efficiency (n* novel samples to hit target)
# ---------------------------------------------------------------------------


def _summarize_cgcd(rows: List[dict], *, mode: str) -> Dict:
    def _mean(key):
        vals = [r[key] for r in rows if r.get(key) is not None and r[key] == r[key]]
        return float(np.mean(vals)) if vals else float("nan")

    return {
        "board": "cifar100_cgcd_proxy",
        "mode": mode,
        "stages": rows,
        "mean_known_acc": _mean("known_acc"),
        "mean_novel_acc": _mean("novel_acc"),
        "mean_forgetting": _mean("forgetting"),
        "mean_feat_hit": _mean("feat_hit"),
        "mean_vimp_top_mass": _mean("vimp_top_mass"),
    }


def run_cgcd(
    *,
    seed: int = SEED,
    mode: str = "weighted",
    quick: bool = False,
) -> Dict:
    """mode: weighted | membership | uniform

    weighted    — subspace membership ∪ private-live labels
    membership  — subspace membership only (no label fusion)
    uniform     — flat partial_fit
    """
    stream = make_cgcd_stream(
        n_known=10 if quick else 12,
        n_novel_per_stage=2 if quick else 3,
        n_stages=3 if quick else 4,
        n_per=24 if quick else 40,
        p=24 if quick else 32,
        novel_dims=6 if quick else 8,
        seed=seed,
    )
    clf, sc = fit_initial(
        stream["X_ref"], stream["y_ref"], seed=seed, classes=stream["all_classes"],
    )
    X_ref0 = sc.transform(stream["X_ref"])  # pinned REF for forgetting
    y_ref0 = np.asarray(stream["y_ref"]).ravel()
    X_ref, y_ref = X_ref0.copy(), y_ref0.copy()
    known = stream["known"]
    classes = np.asarray(stream["all_classes"], dtype=int)
    rows = []
    for st in stream["stages"]:
        Xl, yl = sc.transform(st["X"]), np.asarray(st["y"]).ravel()
        # train/eval split inside the stage (no train=test)
        rng_e = np.random.default_rng(seed + 100 + st["stage"])
        perm = rng_e.permutation(len(yl))
        n_ev = max(16, len(yl) // 3)
        ev, tr = perm[:n_ev], perm[n_ev:]
        Xtr, ytr, Xev, yev = Xl[tr], yl[tr], Xl[ev], yl[ev]

        if mode == "uniform":
            clf = adapt_uniform_partial_fit(clf, Xtr, ytr, classes)
            meta = {
                "reason": "uniform", "vimp_top_mass": None, "feat_hit": None,
                "boundary": "n/a", "mean_membership": None,
            }
        else:
            out = adapt_batch(
                clf, X_ref, y_ref, Xtr, ytr, seed=seed + st["stage"],
                top_k=8, lam=2.0, min_jaccard=0.05,
                detect_thresh=1.05, force_adapt=True,
                use_private_labels=(mode == "weighted"),
                weight_mode="fused" if mode == "weighted" else "sample",
                temper=True,
            )
            clf = out["clf"]
            meta = {
                "reason": out["reason"],
                "vimp_top_mass": out["vimp_top_mass"],
                "feat_hit": float(np.isin(out["feat_idx"], st["active_dims"]).mean()),
                "boundary": out["boundary"]["decision"],
                "mean_membership": out["mean_membership"],
                "subset_frac": out.get("subset_frac"),
                "temper_shrink": (out.get("weights") or {}).get("temper", {}).get("shrink"),
            }

        pred = clf.predict(Xev)
        known_mask, novel_mask = np.isin(yev, known), np.isin(yev, st["novel"])
        known_acc = float(np.mean(pred[known_mask] == yev[known_mask])) if known_mask.any() else float("nan")
        novel_acc = float(np.mean(pred[novel_mask] == yev[novel_mask])) if novel_mask.any() else float("nan")
        forget = 1.0 - float(np.mean(clf.predict(X_ref0) == y_ref0))
        rows.append({
            "stage": st["stage"], "known_acc": known_acc, "novel_acc": novel_acc,
            "forgetting": forget, "n_novel_train": int(np.isin(ytr, st["novel"]).sum()),
            **meta,
        })
        # replay known only into working REF (forgetting still on X_ref0)
        known_idx = np.flatnonzero(np.isin(ytr, known))
        if len(known_idx):
            take = known_idx[::2][: max(8, len(known_idx) // 2)]
            X_ref = np.vstack([X_ref, Xtr[take]])
            y_ref = np.concatenate([y_ref, ytr[take]])
    return _summarize_cgcd(rows, mode=mode)


def cold_start_efficiency(
    *,
    seed: int = SEED,
    budgets: Optional[List[int]] = None,
    target: float = 0.75,
    quick: bool = False,
) -> Dict:
    """n* = min novel-sample budget s.t. novel_acc ≥ target (holdout novel)."""
    budgets = budgets or ([4, 8, 12, 16, 24, 32] if not quick else [4, 8, 12, 16, 24])
    stream = make_cgcd_stream(
        n_known=10 if quick else 12, n_novel_per_stage=3, n_stages=1,
        n_per=48 if quick else 72, p=24 if quick else 32,
        novel_dims=6 if quick else 8, seed=seed,
    )
    st = stream["stages"][0]
    novel, known = st["novel"], stream["known"]
    classes = np.asarray(stream["all_classes"], dtype=int)
    _, sc = fit_initial(stream["X_ref"], stream["y_ref"], seed=seed, classes=classes)
    Xl, yl = sc.transform(st["X"]), np.asarray(st["y"]).ravel()
    X_ref = sc.transform(stream["X_ref"])
    y_ref = np.asarray(stream["y_ref"]).ravel()
    novel_idx = np.flatnonzero(np.isin(yl, novel))
    known_idx = np.flatnonzero(np.isin(yl, known))
    rng = np.random.default_rng(seed + 7)
    rng.shuffle(novel_idx)
    n_hold = max(8, len(novel_idx) // 3)
    hold, pool = novel_idx[:n_hold], novel_idx[n_hold:]

    curves = {}
    for mode in ("weighted", "membership", "uniform"):
        pts, n_star = [], None
        for b in budgets:
            clf, _ = fit_initial(stream["X_ref"], stream["y_ref"], seed=seed, classes=classes)
            take_n = pool[: min(int(b), len(pool))]
            take_k = known_idx[: max(12, 3 * len(take_n))]
            tr = np.concatenate([take_n, take_k])
            Xtr, ytr = Xl[tr], yl[tr]
            if mode == "uniform":
                clf = adapt_uniform_partial_fit(clf, Xtr, ytr, classes)
            else:
                out = adapt_batch(
                    clf, X_ref, y_ref, Xtr, ytr, seed=seed + int(b),
                    top_k=8, lam=2.0, min_jaccard=0.05, force_adapt=True,
                    use_private_labels=(mode == "weighted"),
                    weight_mode="fused" if mode == "weighted" else "sample",
                    temper=True,
                )
                clf = out["clf"]
            acc = float(np.mean(clf.predict(Xl[hold]) == yl[hold]))
            pts.append({"budget": int(b), "novel_acc": acc, "n_novel_train": int(len(take_n))})
            if n_star is None and acc >= float(target):
                n_star = int(b)
        curves[mode] = {"points": pts, "n_star": n_star, "target": float(target)}

    def _acc_at(mode, n):
        for p in curves[mode]["points"]:
            if p["budget"] == n:
                return float(p["novel_acc"])
        return -1.0

    reached = [
        (m, curves[m]["n_star"], _acc_at(m, curves[m]["n_star"]))
        for m in curves if curves[m]["n_star"] is not None
    ]
    # prefer smaller n*; break ties by higher novel_acc at n*
    winner = min(reached, key=lambda t: (t[1], -t[2]))[0] if reached else None
    return {
        "board": "cold_start_efficiency",
        "target": float(target),
        "budgets": list(budgets),
        "curves": curves,
        "winner": winner,
        "tie_break": "min n*, then max novel_acc at n*",
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
        _, sc = fit_initial(src["X"], src["y"], seed=seed)
        Xs, ys = sc.transform(src["X"]), src["y"]
        for tgt_name in ("Clipart", "Product", "Real"):
            clf, _ = fit_initial(src["X"], src["y"], seed=seed)
            tgt = board["domains"][tgt_name]
            Xt, yt = sc.transform(tgt["X"]), tgt["y"]
            bound = boundary_check(ys, yt, min_jaccard=0.20)
            out = adapt_batch(
                clf, Xs, ys, Xt, yt, seed=seed, top_k=6, lam=2.5,
                min_jaccard=0.20, detect_thresh=1.05, force_adapt=False,
                use_private_labels=True, weight_mode="fused",
            )
            shared = board["shared"]
            m = np.isin(yt, shared)
            shared_acc = float(np.mean(out["clf"].predict(Xt)[m] == yt[m])) if m.any() else float("nan")
            rows.append({
                "overlap_frac": ov, "target": tgt_name, "jaccard": bound["jaccard"],
                "boundary": bound["decision"], "adapted": out["adapted"],
                "reason": out["reason"], "shared_acc": shared_acc,
                "vimp_top_mass": out["vimp_top_mass"],
            })
    by_ov = {}
    for ov in overlap_fracs:
        sub = [r for r in rows if r["overlap_frac"] == ov]
        by_ov[str(ov)] = {
            "mean_shared_acc": float(np.nanmean([r["shared_acc"] for r in sub])),
            "incomparable_rate": float(np.mean([r["boundary"] == "incomparable" for r in sub])),
            "adapt_rate": float(np.mean([r["adapted"] for r in sub])),
        }
    return {"board": "officehome_pda_proxy", "rows": rows, "by_overlap": by_ov}


def run_suite(*, seed: int = SEED, quick: bool = True) -> Dict:
    cg_w = run_cgcd(seed=seed, mode="weighted", quick=quick)
    cg_m = run_cgcd(seed=seed, mode="membership", quick=quick)
    cg_u = run_cgcd(seed=seed, mode="uniform", quick=quick)
    eff = cold_start_efficiency(seed=seed, quick=quick, target=0.60)
    oh = run_officehome_pda(seed=seed, quick=quick)
    return {
        "protocol": {
            "axes": [
                "novel_acc (sparse discovery)",
                "known retention / forgetting on pinned REF",
                "attribution: feat_hit + vimp_top_mass",
                "cold-start efficiency n* (budget → target novel_acc)",
                "boundary: incomparable rate vs label overlap",
            ],
            "modes": {
                "weighted": "membership ∪ private-live labels",
                "membership": "subspace membership only",
                "uniform": "flat partial_fit",
            },
        },
        "cgcd_weighted": cg_w,
        "cgcd_membership": cg_m,
        "cgcd_uniform": cg_u,
        "efficiency": eff,
        "officehome": oh,
        "delta_novel_acc_vs_uniform": cg_w["mean_novel_acc"] - cg_u["mean_novel_acc"],
        "delta_novel_acc_membership": cg_m["mean_novel_acc"] - cg_u["mean_novel_acc"],
        "delta_forgetting_vs_uniform": cg_u["mean_forgetting"] - cg_w["mean_forgetting"],
        "takeaway": (
            "Directed adapt = gate + boundary + upweight localized/private samples. "
            "C-GCD reports novel/forget/feat_hit + n* efficiency; "
            "PDA reports incomparable rate rising as overlap falls."
        ),
    }
