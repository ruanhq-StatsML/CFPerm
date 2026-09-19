"""Recsys order-stream sliced by grain: OnlineRFPerm + RFPerm/CFPerm + FSDS.

X is order / merchant / user / all. Y is conversion, never a feature.
T is the batch label. Localization, not a unique decomposition.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from graph_fsds_localize import (
    assert_not_outcome,
    cmean_x,
    cmean_y,
    fsds_rank_columns,
    slice_stream,
    three_metrics,
)
from online_fdr import addis, saffron
from online_rfperm import FrozenRFPerm, onset_from_rows, probe_mse
from stream_dgps import GRAIN_FEATS, planted_how_for_kind
from streaming_po_risk import TabularPORisk, pack_ref_new, rbf_bandwidth, mmd_vs_reference

DIMS = ("order", "merchant", "user", "all")
PO_MIN_N = 40

# DGP: amount walks by `shift`, channel by 0.7. GMV is not planted on the row —
# it is a past-amount EMA, so it follows amount on south merchants.
# Concept plants amount in Y|X only. Ranking recovery uses this order, not Shapley.
PLANTED_MAG = {
    "covariate_south": {"amount": 1.0, "channel": 0.7, "merchant_gmv": 0.8},
    "concept_south": {"amount": 1.0},
    "both": {"amount": 1.0, "channel": 0.7, "merchant_gmv": 0.8},
}


def planted_in_grain(kind: str, names: Sequence[str]) -> dict[str, float]:
    mag = dict(PLANTED_MAG.get(str(kind)) or {})
    if not mag:
        mag = {k: 1.0 for k in planted_how_for_kind(kind)}
    return {n: float(mag[n]) for n in names if n in mag}


def mark_delay(t_hat, onset_true: int) -> dict:
    """OnlineRFPerm paper: delay = first rejection − labeled onset. Negative = FAR."""
    if t_hat is None:
        return {"t": None, "delay": None, "status": "miss"}
    t_hat = int(t_hat)
    d = t_hat - int(onset_true)
    return {"t": t_hat, "delay": int(d), "status": "FAR" if d < 0 else "hit"}


def kendall_vs_planted(rank_rows: Sequence[Mapping], planted: Mapping[str, float], *, score_key: str) -> float | None:
    """Kendall-τ of estimated scores vs planted magnitudes (MetaLearner §ranking)."""
    if not planted:
        return None
    from scipy.stats import kendalltau

    scores = []
    truth = []
    for r in rank_rows:
        name = str(r["feature"])
        scores.append(float(r[score_key]))
        truth.append(float(planted.get(name, 0.0)))
    if len(set(truth)) < 2:
        return None
    tau, _ = kendalltau(scores, truth)
    if tau is None or not np.isfinite(tau):
        return None
    return float(tau)


def _grain_block(window: Mapping, dim: str) -> tuple[np.ndarray, tuple[str, ...]]:
    """One native grain. Columns are the catalog for that grain only."""
    feats = GRAIN_FEATS[dim]
    X = np.asarray(window[f"X_{dim}"], dtype=float)
    names = tuple(window[f"names_{dim}"])
    assert_not_outcome(names)
    missing = [f for f in feats if f not in names]
    extra = [n for n in names if n not in feats]
    if missing or extra:
        raise ValueError(f"{dim} grain columns mismatch catalog: extra={extra} missing={missing}")
    idx = [names.index(f) for f in feats]
    return X[:, idx], tuple(feats)


def slice_xy(window: Mapping, dim: str) -> tuple[np.ndarray, tuple[str, ...]]:
    """Native-grain X. Region / id / Y stay out. Grains do not share columns."""
    dim = str(dim)
    if dim in GRAIN_FEATS:
        return _grain_block(window, dim)
    if dim == "all":
        blocks = [_grain_block(window, g) for g in GRAIN_FEATS]
        names = tuple(n for _, ng in blocks for n in ng)
        if len(names) != len(set(names)):
            raise ValueError(f"grain columns overlap: {names}")
        assert_not_outcome(names)
        return np.hstack([X for X, _ in blocks]), names
    raise ValueError(f"dim must be one of {DIMS}, got {dim!r}")


def rfperm_po_vimp(
    X_ref,
    Y_ref,
    X_new,
    Y_new,
    names: Sequence[str],
    *,
    seed: int = 2026,
    n_perm: int = 12,
) -> dict:
    """PermuCATE-style VIMP on frozen φ. T = batch. Not CATE, not Shapley.

    Nuisances fit once. Feature j's v_j = extra MSE of predicting φ after
    permuting column j. CFPerm-style reject: max v_j vs T-permuted null.
    """
    names = tuple(str(n) for n in names)
    assert_not_outcome(names)
    X, Y, T = pack_ref_new(X_ref, Y_ref, X_new, Y_new)
    est = TabularPORisk(seed=seed)
    packed = est.risk(X, Y, T)
    mu = np.asarray(packed["mu"], dtype=float).ravel()
    e = np.asarray(packed["e"], dtype=float).ravel()
    phi = (Y - mu) * (T - e)
    v_obs = _phi_vimp(X, phi, seed=seed)
    rng = np.random.default_rng(int(seed) + 9)
    null_max = []
    for b in range(int(n_perm)):
        Tp = rng.permutation(T)
        phi_p = (Y - mu) * (Tp - e)
        vp = _phi_vimp(X, phi_p, seed=seed + 1 + b)
        null_max.append(max(vp) if vp else 0.0)
    mx = float(max(v_obs) if v_obs else 0.0)
    q95 = float(np.quantile(null_max, 0.95)) if null_max else 0.0
    ranked = sorted(
        [{"feature": names[j], "vimp": float(v_obs[j])} for j in range(len(names))],
        key=lambda r: -r["vimp"],
    )
    return {
        "rank": ranked,
        "po_risk": float(packed["po_risk"]),
        "max_vimp": mx,
        "null_q95": q95,
        "reject": bool(mx > q95 and mx > 1e-8),
        "top": [r["feature"] for r in ranked[:3]],
    }


def rfperm_mse_vimp(model, X_new, Y_new, names: Sequence[str], *, seed: int = 2026) -> dict:
    """PermOOB-style VIMP on the frozen OnlineRFPerm model. Not CATE.

    v_j = MSE(Y, f_ref(X with column j shuffled)) − MSE(Y, f_ref(X)).
    """
    names = tuple(str(n) for n in names)
    assert_not_outcome(names)
    X_new = np.asarray(X_new, dtype=float)
    Y_new = np.asarray(Y_new, dtype=float).ravel()
    base = probe_mse(model, X_new, Y_new)
    rng = np.random.default_rng(int(seed) + 5)
    scores = []
    for j in range(X_new.shape[1]):
        Xp = X_new.copy()
        rng.shuffle(Xp[:, j])
        scores.append(float(probe_mse(model, Xp, Y_new) - base))
    ranked = sorted(
        [{"feature": names[j], "vimp": float(scores[j])} for j in range(len(names))],
        key=lambda r: -r["vimp"],
    )
    return {
        "rank": ranked,
        "base_mse": float(base),
        "top": [r["feature"] for r in ranked[:3]],
    }


def _phi_vimp(X, phi, *, seed: int) -> list[float]:
    from sklearn.ensemble import RandomForestRegressor

    X = np.asarray(X, dtype=float)
    phi = np.asarray(phi, dtype=float).ravel()
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=5, random_state=int(seed), n_jobs=1
    )
    rf.fit(X, phi)
    base = float(np.mean((phi - rf.predict(X)) ** 2))
    rng = np.random.default_rng(int(seed) + 3)
    out = []
    for j in range(X.shape[1]):
        Xp = X.copy()
        rng.shuffle(Xp[:, j])
        out.append(float(np.mean((phi - rf.predict(Xp)) ** 2) - base))
    return out


def posthoc_region(ref: Mapping, new: Mapping, dim: str, *, seed: int = 2026) -> dict:
    """South vs north on this grain. Subset key is region, not Y."""
    X_ref, names = slice_xy(ref, dim)
    X_new, _ = slice_xy(new, dim)
    Y_ref = np.asarray(ref["Y"], dtype=float).ravel()
    Y_new = np.asarray(new["Y"], dtype=float).ravel()
    south = np.asarray(new["region"]) == "south"
    north = ~south
    sigma = rbf_bandwidth(X_ref, seed=seed)
    out = {
        "names": list(names),
        "south": three_metrics(
            X_ref[np.asarray(ref["region"]) == "south"] if (np.asarray(ref["region"]) == "south").any() else X_ref,
            Y_ref[np.asarray(ref["region"]) == "south"] if (np.asarray(ref["region"]) == "south").any() else Y_ref,
            X_new[south],
            Y_new[south],
            sigma=sigma,
            seed=seed,
            with_po=int(south.sum()) >= PO_MIN_N,
        )
        if int(south.sum()) >= 8
        else {},
        "north": three_metrics(
            X_ref[np.asarray(ref["region"]) == "north"] if (np.asarray(ref["region"]) == "north").any() else X_ref,
            Y_ref[np.asarray(ref["region"]) == "north"] if (np.asarray(ref["region"]) == "north").any() else Y_ref,
            X_new[north],
            Y_new[north],
            sigma=sigma,
            seed=seed,
            with_po=int(north.sum()) >= PO_MIN_N,
        )
        if int(north.sum()) >= 8
        else {},
        "pair_mmd": float(mmd_vs_reference(X_new[north], X_new[south], sigma=sigma, seed=seed))
        if int(south.sum()) >= 8 and int(north.sum()) >= 8
        else None,
        "pair_cmean_x": cmean_x(X_new[south], X_new[north])
        if int(south.sum()) >= 8 and int(north.sum()) >= 8
        else None,
        "pair_cmean_y": cmean_y(Y_new[south], Y_new[north])
        if int(south.sum()) >= 8 and int(north.sum()) >= 8
        else None,
        "n_south": int(south.sum()),
        "n_north": int(north.sum()),
    }
    return out


def run_dim_stream(
    tables: Mapping,
    dim: str,
    *,
    seed: int = 2026,
    with_po: bool = True,
) -> dict:
    """OnlineRFPerm on one grain, then last-batch RFPerm/FSDS/post-hoc."""
    meta = tables["meta"]
    n_batches = int(meta["n_batches"])
    onset_true = int(meta["onset_batch"])
    cut0 = slice_stream(tables, 0)
    X_ref, names = slice_xy(cut0["ref"], dim)
    Y_ref = np.asarray(cut0["ref"]["Y"], dtype=float).ravel()
    probe = FrozenRFPerm(X_ref, Y_ref, seed=seed)
    sigma = rbf_bandwidth(X_ref, seed=seed)
    rows = []
    last_new = None
    for t in range(n_batches):
        cut = slice_stream(tables, t)
        X_new, _ = slice_xy(cut["new"], dim)
        Y_new = np.asarray(cut["new"]["Y"], dtype=float).ravel()
        rec = probe.step(X_new, Y_new)
        rec["t"] = int(t)
        rec["onset"] = t >= onset_true
        rec["mse"] = rec["rfperm_mse"]
        rec["mmd"] = float(mmd_vs_reference(X_ref, X_new, sigma=sigma, seed=seed))
        rec["cmean_x"] = cmean_x(X_new, X_ref)
        rec["cmean_y"] = cmean_y(Y_new, Y_ref)
        rec["po"] = None
        if with_po and t == n_batches - 1 and len(X_new) >= PO_MIN_N:
            mets = three_metrics(X_ref, Y_ref, X_new, Y_new, sigma=sigma, seed=seed, with_po=True)
            rec["po"] = mets["po"]
        rows.append(rec)
        last_new = cut["new"]
    p_fdr = [float(r["rfperm_p"]) if float(r["rfperm_T"]) > 0 else 1.0 for r in rows]
    ad = addis(p_fdr)
    sf = saffron(p_fdr)
    onset = onset_from_rows(rows)
    addis_t = next((int(rows[i]["t"]) for i, r in enumerate(ad["reject"]) if r), None)
    saffron_t = next((int(rows[i]["t"]) for i, r in enumerate(sf["reject"]) if r), None)
    fsds = fsds_rank_columns(
        X_ref,
        slice_xy(last_new, dim)[0],
        Y_ref,
        np.asarray(last_new["Y"], dtype=float).ravel(),
        names,
        seed=seed,
    )
    X_last, _ = slice_xy(last_new, dim)
    Y_last = np.asarray(last_new["Y"], dtype=float).ravel()
    vimp = rfperm_po_vimp(X_ref, Y_ref, X_last, Y_last, names, seed=seed)
    mse_vimp = rfperm_mse_vimp(probe.probe, X_last, Y_last, names, seed=seed)
    loc = posthoc_region(cut0["ref"], last_new, dim, seed=seed)
    planted = planted_in_grain(meta["kind"], names)
    recovered = [r["feature"] for r in fsds if r.get("loud") and r["feature"] in planted]
    hop_mark = mark_delay(onset["onset_hat"], onset_true)
    rank_mark = mark_delay(onset["onset_rank"], onset_true)
    addis_mark = mark_delay(addis_t, onset_true)
    saffron_mark = mark_delay(saffron_t, onset_true)
    return {
        "dim": dim,
        "kind": meta["kind"],
        "names": list(names),
        "onset_true": onset_true,
        "onset_hat": onset["onset_hat"],
        "onset_rank": onset["onset_rank"],
        "addis_t": addis_t,
        "saffron_t": saffron_t,
        "n_hop": onset["n_rfperm_hop"],
        "hop_delay": hop_mark["delay"],
        "hop_status": hop_mark["status"],
        "rank_delay": rank_mark["delay"],
        "rank_status": rank_mark["status"],
        "addis_delay": addis_mark["delay"],
        "addis_status": addis_mark["status"],
        "saffron_delay": saffron_mark["delay"],
        "saffron_status": saffron_mark["status"],
        "rows": rows,
        "fsds": fsds,
        "vimp": vimp,
        "mse_vimp": mse_vimp,
        "localization": loc,
        "planted": planted,
        "fsds_recovered": recovered,
        "vimp_top": vimp["top"],
        "vimp_reject": vimp["reject"],
        "mse_vimp_top": mse_vimp["top"],
        "tau_fsds": kendall_vs_planted(fsds, planted, score_key="score"),
        "tau_cfperm": kendall_vs_planted(vimp["rank"], planted, score_key="vimp"),
        "tau_rfperm": kendall_vs_planted(mse_vimp["rank"], planted, score_key="vimp"),
        "y_in_X": any(n in ("y", "Y", "label") for n in names),
    }


def run_all_dims(tables: Mapping, *, seed: int = 2026, with_po: bool = True) -> dict:
    return {dim: run_dim_stream(tables, dim, seed=seed, with_po=with_po) for dim in DIMS}
