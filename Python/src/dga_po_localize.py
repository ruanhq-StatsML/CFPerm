"""PO-risk post-hoc subset localization for DGA hop weights.

DGA answers a *training* question: which past domains should be mixed so that
one step most reduces L_spe. Subset localization answers a *diagnostic*
question: which pocket of X is driving the discrepancy. They share a slot
only after PO-risk localizes D_spe itself.

Pipeline (no new loss, ranking only — not a CFPerm p-value):

  1. Instance PO-risk on hop t-1 vs past (DR pseudo-outcome energy φ²).
  2. D_spe^loc = PO tail inside B_{t-1} (high-risk specialized pocket).
  3. DGA alignments vs ∇ℓ(θ, D_spe^loc), not the whole last batch.
  4. Feature mean-diff / Cohen's d (high-PO vs complement in B_{t-1})
     discretizes the subset — the README post-hoc subgroup step.

Locked signs stay with TSS (ĉ, δ̂). This module only changes how D_spe
is formed and how the pocket is described after DGA.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import Ridge

from amazon_continuous_batches import AmazonStream, SEED
from attribution_adapter import (
    _linear_mse_grad,
    _mse,
    _ridge_fit,
    _ridge_predict,
    _summary,
    dga_alignments,
    dga_mirror_step,
    dga_sample_weights,
)


def dr_pseudo_outcome(X, y, t, clip=0.05, ridge_alpha=3.0):
    """Doubly-robust pseudo-outcome φ for batch treatment T.

    φ = (μ1−μ0) + (T−π)(Y−μ_T) / (π(1−π)).

    π is the known batch share clip(mean(T)) — hop sizes are observed, so a
    logistic π is extra variance on small n. μ0 / μ1 are Ridge outcome models.
    Localization ranking, not a cross-fit inference target.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    t = np.asarray(t, dtype=int).ravel()
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    n0 = int((t == 0).sum())
    n1 = int((t == 1).sum())
    if n0 < 2 or n1 < 2:
        zeros = np.zeros(y.shape[0], dtype=float)
        pi = float(np.clip(n1 / max(n0 + n1, 1), clip, 1.0 - clip))
        return zeros, {"pi": pi, "mu0": zeros, "mu1": zeros, "ok": False}
    mu0 = Ridge(alpha=float(ridge_alpha)).fit(X[t == 0], y[t == 0]).predict(X)
    mu1 = Ridge(alpha=float(ridge_alpha)).fit(X[t == 1], y[t == 1]).predict(X)
    pi = float(np.clip(n1 / float(n0 + n1), clip, 1.0 - clip))
    mu_t = np.where(t == 1, mu1, mu0)
    phi = (mu1 - mu0) + (t - pi) * (y - mu_t) / (pi * (1.0 - pi))
    return np.asarray(phi, dtype=float), {
        "pi": pi,
        "mu0": np.asarray(mu0, dtype=float),
        "mu1": np.asarray(mu1, dtype=float),
        "ok": True,
    }


def instance_po_risk(phi, X, ridge_alpha=3.0, mode="phi2", mu0=None, y=None):
    """Instance score used to rank the specialized pocket.

    ``phi2`` (default): r_i = φ_i², the instance DR contrast. This is the
    localization score — large batch contrast at x.

    Do *not* use (φ − τ̂)² for subset localization: after a global τ̂ fits,
    leftover CATE residual often lands on the complement, anti-localizing
    the planted pocket. That quantity is estimator risk, not a shift map.

    ``control_resid``: (Y − μ0(X))², instance concept residual vs the past
    outcome — the row-level analogue of δ̂.
    """
    phi = np.asarray(phi, dtype=float).ravel()
    if mode == "control_resid":
        if mu0 is None or y is None:
            raise ValueError("control_resid needs mu0 and y")
        return np.square(np.asarray(y, dtype=float).ravel() - np.asarray(mu0, dtype=float).ravel())
    if mode == "cate_resid":
        X = np.asarray(X, dtype=float)
        if X.shape[0] < 4:
            return np.square(phi)
        tau = Ridge(alpha=float(ridge_alpha)).fit(X, phi).predict(X)
        return np.square(phi - tau)
    return np.square(phi)


def po_tail_mask(is_spe, risk, q=0.30, min_n=4):
    """Boolean mask: high-PO tail inside the specialized set."""
    is_spe = np.asarray(is_spe, dtype=bool).reshape(-1)
    risk = np.asarray(risk, dtype=float).ravel()
    out = np.zeros(is_spe.shape[0], dtype=bool)
    idx = np.flatnonzero(is_spe)
    if idx.size == 0:
        return out
    r = risk[idx]
    keep = max(int(min_n), int(np.ceil(float(q) * idx.size)))
    keep = min(keep, idx.size)
    if keep >= idx.size:
        out[idx] = True
        return out
    order = np.argpartition(r, -keep)[-keep:]
    out[idx[order]] = True
    return out


def cohens_d_coord(a, b):
    """Cohen's d for high-group minus low-group (signed)."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if a.size < 2 or b.size < 2:
        return 0.0
    sp2 = ((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1)) / max(
        a.size + b.size - 2, 1
    )
    return float((a.mean() - b.mean()) / np.sqrt(sp2 + 1e-12))


def subgroup_mean_diff(X, high_mask, names=None):
    """README post-hoc: mean difference / Cohen's d of high-PO vs complement."""
    X = np.asarray(X, dtype=float)
    high_mask = np.asarray(high_mask, dtype=bool).reshape(-1)
    if names is None:
        names = ["X%d" % j for j in range(X.shape[1])]
    rows = []
    for j, name in enumerate(names):
        hi = X[high_mask, j] if np.any(high_mask) else np.array([])
        lo = X[~high_mask, j] if np.any(~high_mask) else np.array([])
        mu_h = float(hi.mean()) if hi.size else float("nan")
        mu_l = float(lo.mean()) if lo.size else float("nan")
        rows.append(
            {
                "j": int(j),
                "name": str(name),
                "mean_high": mu_h,
                "mean_low": mu_l,
                "mean_diff": float(mu_h - mu_l) if hi.size and lo.size else 0.0,
                "cohens_d": cohens_d_coord(hi, lo),
                "n_high": int(hi.size),
                "n_low": int(lo.size),
            }
        )
    rows.sort(key=lambda r: abs(r["cohens_d"]), reverse=True)
    for rank, rec in enumerate(rows, start=1):
        rec["rank"] = rank
    return rows


def discretize_subset(ranked, k=3):
    """Turn top-k mean-diff features into a conjunction of midpoint splits.

    For feature j, split at 0.5 (mean_high + mean_low). Sign of Cohen's d
    picks {X_j ≥ c} vs {X_j ≤ c}. Intersection is the localized slice.
    """
    rules = []
    for rec in ranked[: int(k)]:
        if rec["n_high"] < 2 or rec["n_low"] < 2:
            continue
        c = 0.5 * (rec["mean_high"] + rec["mean_low"])
        ge = rec["cohens_d"] >= 0.0
        rules.append(
            {
                "j": rec["j"],
                "name": rec["name"],
                "threshold": float(c),
                "op": ">=" if ge else "<=",
                "cohens_d": rec["cohens_d"],
            }
        )
    return rules


def apply_rules(X, rules):
    """Boolean mask for the discretized conjunction (empty rules → all False)."""
    X = np.asarray(X, dtype=float)
    if not rules:
        return np.zeros(X.shape[0], dtype=bool)
    mask = np.ones(X.shape[0], dtype=bool)
    for rule in rules:
        col = X[:, int(rule["j"])]
        thr = float(rule["threshold"])
        if rule["op"] == ">=":
            mask &= col >= thr
        else:
            mask &= col <= thr
    return mask


def domain_po_table(batch, risk, alignments, alpha, domains):
    """Per-domain mean PO-risk joined with DGA a_i / α_i."""
    batch = np.asarray(batch, dtype=int)
    risk = np.asarray(risk, dtype=float).ravel()
    a = np.asarray(alignments, dtype=float).ravel()
    al = np.asarray(alpha, dtype=float).ravel()
    rows = []
    for j, s in enumerate(domains):
        idx = batch == int(s)
        rows.append(
            {
                "domain": int(s),
                "n": int(idx.sum()),
                "mean_po": float(risk[idx].mean()) if np.any(idx) else 0.0,
                "alignment": float(a[j]) if j < a.size else float("nan"),
                "alpha": float(al[j]) if j < al.size else float("nan"),
            }
        )
    return rows


def make_planted_subgroup_stream(
    n_batches=5,
    n_per=120,
    p=16,
    seed=SEED,
    subgroup_feat=0,
    subgroup_thr=0.35,
    shift_coords=(0, 1, 2),
    mean_bump=1.15,
    concept_flip=True,
    noise=0.40,
    plant_from=3,
    rank=4,
):
    """Past is a stable linear rating map. From ``plant_from`` onward, a
    pocket {X[subgroup_feat] > thr} gets a mean bump on ``shift_coords``
    and (optionally) a concept flip on those coordinates.

    The planted mask is stored on ``meta['planted_mask']`` so localization
    recovery is checkable.
    """
    rng = np.random.default_rng(seed)
    p = int(p)
    beta = np.zeros(p)
    beta[: int(rank)] = 0.85
    shift = np.array(list(shift_coords), dtype=int)
    planted = []
    rows, ys, batches = [], [], []
    for t in range(int(n_batches)):
        X = rng.normal(size=(int(n_per), p))
        y = 3.0 + X @ beta + rng.normal(scale=float(noise), size=int(n_per))
        mask = np.zeros(int(n_per), dtype=bool)
        if t >= int(plant_from):
            mask = X[:, int(subgroup_feat)] > float(subgroup_thr)
            if np.any(mask):
                X[np.ix_(mask, shift)] = X[np.ix_(mask, shift)] + float(mean_bump)
                if concept_flip:
                    y[mask] = (
                        3.0
                        - X[mask] @ beta
                        + rng.normal(scale=float(noise), size=int(mask.sum()))
                    )
        y = np.clip(y, 1.0, 5.0)
        rows.append(X)
        ys.append(y)
        batches.append(np.full(int(n_per), t, dtype=int))
        planted.append(mask)
    planted_mask = np.concatenate(planted)
    return AmazonStream(
        X=np.vstack(rows),
        y=np.concatenate(ys),
        batch=np.concatenate(batches),
        categories=tuple("plant_%d" % t for t in range(int(n_batches))),
        meta={
            "source": "planted_subgroup",
            "n_batches": int(n_batches),
            "n_per": int(n_per),
            "p": int(p),
            "seed": int(seed),
            "subgroup_feat": int(subgroup_feat),
            "subgroup_thr": float(subgroup_thr),
            "shift_coords": [int(j) for j in shift],
            "plant_from": int(plant_from),
            "planted_mask": planted_mask,
            "planted_frac": float(planted_mask.mean()),
        },
    )


def localize_one_hop(X, y, batch, spe, q=0.30, ridge_alpha=3.0, k_features=3):
    """PO-risk tail + mean-diff subgroup on one causal hop (spe = t-1)."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    batch = np.asarray(batch, dtype=int)
    t = (batch == int(spe)).astype(int)
    past = batch < int(spe)
    on_hop = past | (batch == int(spe))
    phi, nuis = dr_pseudo_outcome(X[on_hop], y[on_hop], t[on_hop], ridge_alpha=ridge_alpha)
    risk_hop = instance_po_risk(phi, X[on_hop], ridge_alpha=ridge_alpha, mode="phi2")
    risk = np.zeros(batch.shape[0], dtype=float)
    risk[on_hop] = risk_hop
    is_spe = batch == int(spe)
    tail = po_tail_mask(is_spe, risk, q=q)
    spe_high = tail[is_spe]
    ranked = subgroup_mean_diff(X[is_spe], spe_high)
    rules = discretize_subset(ranked, k=k_features)
    slice_mask = np.zeros(batch.shape[0], dtype=bool)
    slice_mask[is_spe] = apply_rules(X[is_spe], rules)
    resid0 = (
        np.square(y[on_hop] - nuis["mu0"])
        if nuis["ok"]
        else np.zeros(int(on_hop.sum()), dtype=float)
    )
    spe_po = float(risk[is_spe].mean()) if np.any(is_spe) else 0.0
    tail_po = float(risk[tail].mean()) if np.any(tail) else 0.0
    spe_on = t[on_hop] == 1
    return {
        "spe": int(spe),
        "risk": risk,
        "tail_mask": tail,
        "ranked_features": ranked,
        "rules": rules,
        "slice_mask": slice_mask,
        "n_spe": int(is_spe.sum()),
        "n_tail": int(tail.sum()),
        "mean_po_spe": spe_po,
        "mean_po_tail": tail_po,
        "pi": nuis["pi"],
        "ok": nuis["ok"],
        "mean_resid0_spe": float(resid0[spe_on].mean()) if nuis["ok"] and np.any(spe_on) else 0.0,
    }


def recovery_against_plant(stream, tail_mask, spe):
    """Precision / recall of the PO tail vs the planted pocket on B_spe."""
    plant = np.asarray(stream.meta.get("planted_mask"), dtype=bool)
    batch = np.asarray(stream.batch, dtype=int)
    on = batch == int(spe)
    if plant.size != batch.size or not np.any(on):
        return {"precision": float("nan"), "recall": float("nan"), "n_plant": 0}
    gold = plant[on]
    pred = np.asarray(tail_mask, dtype=bool)[on]
    tp = float(np.sum(gold & pred))
    prec = tp / max(float(pred.sum()), 1.0)
    rec = tp / max(float(gold.sum()), 1.0)
    return {
        "precision": prec,
        "recall": rec,
        "n_plant": int(gold.sum()),
        "n_tail": int(pred.sum()),
        "n_spe": int(on.sum()),
    }


@dataclass
class DGAPoConfig:
    eta: float = 1.0
    ema_beta: float = 0.35
    align: str = "cosine"
    ridge_alpha: float = 3.0
    po_q: float = 0.30
    k_features: int = 3


def run_dga_po_localize(stream, cfg=None, **kw):
    """DGA Ridge with PO-tail D_spe + post-hoc mean-diff localization.

    Same adapter slot as ``dga_ridge``: sample-weighted Ridge, no new loss.
    The specialized gradient is taken on the high-PO pocket of B_{t-1}.
    """
    cfg = cfg or DGAPoConfig()
    for key, val in kw.items():
        if hasattr(cfg, key):
            setattr(cfg, key, val)
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    ra = float(cfg.ridge_alpha)
    history = []
    alpha_ema = None
    for t in range(1, k):
        tr = batch < t
        ic = batch == t
        domains = list(range(t))
        spe = t - 1
        loc = localize_one_hop(
            X,
            y,
            batch,
            spe,
            q=cfg.po_q,
            ridge_alpha=ra,
            k_features=cfg.k_features,
        )
        if alpha_ema is None:
            probe_w = None
        else:
            prev = np.asarray(alpha_ema, dtype=float)
            if prev.size == t:
                mix = prev
            elif prev.size == t - 1:
                mean_m = float(prev.mean()) if prev.size else 1.0
                mix = np.concatenate([prev, [mean_m]])
                mix = mix / mix.sum()
            else:
                mix = np.ones(t, dtype=float) / float(t)
            probe_w = dga_sample_weights(mix, batch, domains)[tr]
        clf = _ridge_fit(X[tr], y[tr], sample_weight=probe_w, alpha=ra)
        spe_mask_tr = loc["tail_mask"][tr]
        if not np.any(spe_mask_tr):
            spe_mask_tr = (batch[tr] == spe)
        a = dga_alignments(
            clf,
            X[tr],
            y[tr],
            batch[tr],
            domains,
            spe,
            align=cfg.align,
            spe_mask=spe_mask_tr,
        )
        alpha_inst = dga_mirror_step(np.ones(t, dtype=float) / float(t), a, eta=cfg.eta)
        if alpha_ema is None:
            alpha_ema = alpha_inst.copy()
        else:
            b = float(np.clip(cfg.ema_beta, 0.0, 1.0))
            prev = np.asarray(alpha_ema, dtype=float)
            if prev.size == t - 1:
                mean_m = float(prev.mean()) if prev.size else 1.0
                prev = np.concatenate([prev, [mean_m]])
                prev = prev / prev.sum()
            elif prev.size != t:
                prev = np.ones(t, dtype=float) / float(t)
            alpha_ema = (1.0 - b) * prev + b * alpha_inst
            alpha_ema = alpha_ema / alpha_ema.sum()
        w = dga_sample_weights(alpha_ema, batch, domains)
        pred = _ridge_predict(X[tr], y[tr], X[ic], sample_weight=w[tr], alpha=ra)
        recov = recovery_against_plant(stream, loc["tail_mask"], spe)
        domains_tbl = domain_po_table(
            batch[tr], loc["risk"][tr], a, alpha_ema, domains
        )
        top = loc["ranked_features"][: cfg.k_features]
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "alpha_ema": alpha_ema.tolist(),
                "alpha_inst": alpha_inst.tolist(),
                "alignments": a.tolist(),
                "w_mean": float(w[tr].mean()),
                "w_max": float(w[tr].max()),
                "spe_domain": int(spe),
                "n": int(ic.sum()),
                "n_spe": loc["n_spe"],
                "n_tail": loc["n_tail"],
                "mean_po_spe": loc["mean_po_spe"],
                "mean_po_tail": loc["mean_po_tail"],
                "top_features": [
                    {"name": r["name"], "cohens_d": r["cohens_d"], "j": r["j"]}
                    for r in top
                ],
                "rules": loc["rules"],
                "domains": domains_tbl,
                "recovery": recov,
                "slice_frac": float(loc["slice_mask"][batch == spe].mean())
                if loc["n_spe"]
                else 0.0,
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    last = history[-1] if history else {}
    return _summary(
        "dga_po_ridge",
        online,
        history,
        stream,
        extras={
            "eta": float(cfg.eta),
            "ema_beta": float(cfg.ema_beta),
            "align": cfg.align,
            "alpha": ra,
            "po_q": float(cfg.po_q),
            "k_features": int(cfg.k_features),
            "last_top_features": last.get("top_features", []),
            "last_rules": last.get("rules", []),
            "last_recovery": last.get("recovery", {}),
            "method_ref": "DGA + PO-risk post-hoc subset localization",
        },
    )


def top_feature_hit(ranked, coords, k=3):
    """Whether planted shift coordinates appear in the top-k mean-diff list."""
    got = {int(r["j"]) for r in ranked[: int(k)]}
    want = {int(j) for j in coords}
    return {
        "hit": sorted(got & want),
        "miss": sorted(want - got),
        "topk": sorted(got),
        "recall": float(len(got & want) / max(len(want), 1)),
    }
