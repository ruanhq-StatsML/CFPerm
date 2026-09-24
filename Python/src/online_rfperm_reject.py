"""OnlineRFPerm with a selection readout, and a quantile probe.

This is not the style-transfer router.

Is it rejection inference? No. Credit-style reject inference *imputes* Y
on S=0 so a scorecard can train on the full applicant file. We do not
write a fake Y. Last column stays the outcome; non-finite means unshown.

OnlineRFPerm's frozen probe is MSE − E_ref on rows where Y is observed.
If the serving policy S changes who gets a label, complete-case MSE can
hop while P(Y|X) did not. That is a selection hop, P(S=1|X), not concept.

Three frozen probes vs D_ref, same T=batch clock:
  sel   Brier of a frozen e(S|X)          — policy / who is labeled
  x     MMD²(X_new, X_ref) on full X      — covariate, rejects included
  y     IPS-weighted MSE on S=1           — P(Y|X) on the labeled slice

IPS assumes selection on observables, S ⊥ Y | X. That is an assumption,
not a CATE we identify. Shares are not Shapley.

A second probe, different discipline: pinball at τ. MSE is L2 / mean.
A tail hop can move pinball while complete-case MSE stays quiet.

df contract: last column is Y. Non-finite Y ⇒ S=0. Y is never an X column.

    rec = onlinePermOOB_reject(df, ref_batch_size=..., batch_size=...)
    rec = onlinePermOOB_quantile(df, tau=0.9, ...)
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression

from online_rfperm import fit_frozen_rf, probe_mse
from streaming_po_risk import large_deviation, mmd_vs_reference, rbf_bandwidth

CLIP = 0.05
IPS_CAP = 10.0


def split_xy_s(df):
    """Last column Y. S = 1{Y finite}. X is every column except Y."""
    df = np.asarray(df, dtype=float)
    if df.ndim != 2 or df.shape[1] < 2:
        raise ValueError("df needs X columns and Y last")
    X = df[:, :-1]
    Y = df[:, -1]
    S = np.isfinite(Y).astype(float)
    return X, Y, S


def clip_e(e, lo: float = CLIP):
    e = np.asarray(e, dtype=float).ravel()
    return np.clip(e, float(lo), 1.0 - float(lo))


def fit_propensity(X, S, seed: int = 2026):
    X = np.asarray(X, dtype=float)
    S = np.asarray(S, dtype=float).ravel().astype(int)
    m = LogisticRegression(max_iter=400, random_state=int(seed))
    if len(np.unique(S)) < 2:
        class _Const:
            def predict_proba(self, X_new):
                p = float(S.mean()) if len(S) else 0.5
                p = float(np.clip(p, CLIP, 1.0 - CLIP))
                n = len(np.asarray(X_new))
                return np.column_stack([np.full(n, 1.0 - p), np.full(n, p)])

        return _Const()
    m.fit(X, S)
    return m


def e_hat(model, X) -> np.ndarray:
    proba = np.asarray(model.predict_proba(np.asarray(X, dtype=float)))
    if proba.ndim == 2 and proba.shape[1] == 2:
        return clip_e(proba[:, 1])
    return clip_e(proba.ravel())


def brier(S, p) -> float:
    S = np.asarray(S, dtype=float).ravel()
    p = np.asarray(p, dtype=float).ravel()
    return float(np.mean((S - p) ** 2))


def complete_case_mse(model, X, Y, S) -> float:
    obs = np.asarray(S, dtype=float).ravel() > 0.5
    if obs.sum() < 2:
        return float("nan")
    return probe_mse(model, np.asarray(X)[obs], np.asarray(Y)[obs])


def ips_mse(model, X, Y, S, e) -> float:
    """Hájek IPS of (Y−μ)² on S=1. Not an imputed reject label."""
    S = np.asarray(S, dtype=float).ravel()
    Y = np.asarray(Y, dtype=float).ravel()
    e = clip_e(e)
    obs = S > 0.5
    if obs.sum() < 2:
        return float("nan")
    pred = np.asarray(model.predict(np.asarray(X)[obs]), dtype=float).ravel()
    err = (Y[obs] - pred) ** 2
    w = np.clip(1.0 / e[obs], 1.0, float(IPS_CAP))
    den = float(w.sum())
    if den <= 1e-12:
        return float("nan")
    return float(np.sum(w * err) / den)


def pinball(y, mu, tau: float) -> float:
    y = np.asarray(y, dtype=float).ravel()
    mu = np.asarray(mu, dtype=float).ravel()
    d = y - mu
    t = float(tau)
    return float(np.mean(np.maximum(t * d, (t - 1.0) * d)))


def localize_call(hop_sel, hop_x, hop_y) -> str:
    """Heuristic on the three hops. Localization, not a unique decomp."""
    if hop_sel and not hop_y:
        return "selection"
    if hop_y and not hop_sel and not hop_x:
        return "concept_on_shown"
    if hop_x and not hop_y:
        return "covariate"
    if hop_y and hop_sel:
        return "both"
    return "keep"


def onlinePermOOB_reject(
    df,
    ref_batch_size=120,
    batch_size=40,
    seed=2026,
    gate=2.0,
):
    """Frozen RF on shown D_ref, frozen e(S|X), then the trail.

    Last column is Y. Non-finite Y is S=0, never an X feature.
    """
    X, Y, S = split_xy_s(df)
    n_ref = int(ref_batch_size)
    bs = int(batch_size)
    Xr, Yr, Sr = X[:n_ref], Y[:n_ref], S[:n_ref]
    cut = max(int(0.7 * n_ref), 16)
    shown_fit = (Sr > 0.5) & (np.arange(n_ref) < cut)
    shown_hold = (Sr > 0.5) & (np.arange(n_ref) >= cut)
    if shown_fit.sum() < 8 or shown_hold.sum() < 4:
        shown_fit = Sr > 0.5
        shown_hold = shown_fit
    y_model = fit_frozen_rf(Xr[shown_fit], Yr[shown_fit], seed=seed)
    s_model = fit_propensity(Xr[:cut], Sr[:cut], seed=seed)
    e_hold = e_hat(s_model, Xr[cut:] if cut < n_ref else Xr)
    S_hold = Sr[cut:] if cut < n_ref else Sr
    mse_cc_ref = complete_case_mse(y_model, Xr[shown_hold], Yr[shown_hold], np.ones(shown_hold.sum()))
    mse_ips_ref = ips_mse(
        y_model,
        Xr[cut:] if cut < n_ref else Xr,
        Yr[cut:] if cut < n_ref else Yr,
        S_hold,
        e_hold,
    )
    rate_ref = float(Sr.mean())
    brier_ref = brier(S_hold, e_hold)
    sigma = rbf_bandwidth(Xr, seed=seed)
    Xr_s = Xr[Sr > 0.5]
    sigma_s = rbf_bandwidth(Xr_s, seed=seed)
    rng = np.random.default_rng(int(seed) + 3)
    i1 = rng.choice(len(Xr), size=min(40, len(Xr)), replace=True)
    i2 = rng.choice(len(Xr), size=min(40, len(Xr)), replace=True)
    mmd_base = max(float(mmd_vs_reference(Xr[i1], Xr[i2], sigma=sigma, seed=seed)), 0.015)
    j1 = rng.choice(len(Xr_s), size=min(40, len(Xr_s)), replace=True)
    j2 = rng.choice(len(Xr_s), size=min(40, len(Xr_s)), replace=True)
    mmd_s_base = max(
        float(mmd_vs_reference(Xr_s[j1], Xr_s[j2], sigma=sigma_s, seed=seed)), 0.015
    )

    cc, ips, sel, mmd, calls = [], [], [], [], []
    hop_cc, hop_ips, hop_sel, hop_x = [], [], [], []
    i = n_ref
    while i < len(X):
        Xb, Yb, Sb = X[i : i + bs], Y[i : i + bs], S[i : i + bs]
        if len(Xb) < 8:
            break
        eb = e_hat(s_model, Xb)
        cc_t = complete_case_mse(y_model, Xb, Yb, Sb)
        ips_t = ips_mse(y_model, Xb, Yb, Sb, eb)
        sel_t = brier(Sb, eb)
        mmd_t = float(mmd_vs_reference(Xr, Xb, sigma=sigma, seed=seed))
        obs_b = Sb > 0.5
        if obs_b.sum() >= 8:
            mmd_s = float(mmd_vs_reference(Xr_s, Xb[obs_b], sigma=sigma_s, seed=seed))
        else:
            mmd_s = 0.0
        h_cc = bool(np.isfinite(cc_t) and large_deviation(cc_t, mse_cc_ref, ratio=float(gate)))
        h_ips = bool(np.isfinite(ips_t) and large_deviation(ips_t, mse_ips_ref, ratio=float(gate)))
        h_sel = bool(
            abs(float(Sb.mean()) - rate_ref) >= 0.15
            and large_deviation(mmd_s, mmd_s_base, ratio=float(gate))
        )
        h_x = bool(large_deviation(mmd_t, mmd_base, ratio=float(gate)))
        cc.append(cc_t)
        ips.append(ips_t)
        sel.append(sel_t)
        mmd.append(mmd_t)
        hop_cc.append(h_cc)
        hop_ips.append(h_ips)
        hop_sel.append(h_sel)
        hop_x.append(h_x)
        calls.append(localize_call(h_sel, h_x, h_ips))
        i += bs

    def _first(flags):
        for t, v in enumerate(flags):
            if v:
                return int(t)
        return -1

    return {
        "mse_cc_ref": float(mse_cc_ref),
        "mse_ips_ref": float(mse_ips_ref),
        "brier_sel_ref": float(brier_ref),
        "MSE_cc_list": np.asarray(cc, dtype=float),
        "MSE_ips_list": np.asarray(ips, dtype=float),
        "Brier_sel_list": np.asarray(sel, dtype=float),
        "MMD_list": np.asarray(mmd, dtype=float),
        "hop_cc": np.asarray(hop_cc, dtype=bool),
        "hop_ips": np.asarray(hop_ips, dtype=bool),
        "hop_sel": np.asarray(hop_sel, dtype=bool),
        "hop_x": np.asarray(hop_x, dtype=bool),
        "call": calls,
        "cc_1": _first(hop_cc),
        "ips_1": _first(hop_ips),
        "sel_1": _first(hop_sel),
        "x_1": _first(hop_x),
        "n_shown_ref": int((Sr > 0.5).sum()),
        "y_in_x": False,
    }


def onlinePermOOB_quantile(
    df,
    tau=0.9,
    ref_batch_size=120,
    batch_size=40,
    seed=2026,
    gate=2.0,
):
    """Frozen mean-MSE and frozen pinball(τ) on a complete Y.

    Last column is Y, must be finite. Tail hops can move pinball first.
    """
    df = np.asarray(df, dtype=float)
    X, Y = df[:, :-1], df[:, -1]
    if not np.all(np.isfinite(Y)):
        raise ValueError("quantile probe wants a complete Y; use reject() for NaNs")
    n_ref = int(ref_batch_size)
    bs = int(batch_size)
    cut = max(int(0.7 * n_ref), 16)
    Xr, Yr = X[:cut], Y[:cut]
    Xh, Yh = X[cut:n_ref], Y[cut:n_ref]
    if len(Xh) < 8:
        Xr, Yr = X[:n_ref], Y[:n_ref]
        Xh, Yh = Xr, Yr
    mean_m = fit_frozen_rf(Xr, Yr, seed=seed)
    q_m = GradientBoostingRegressor(
        loss="quantile",
        alpha=float(tau),
        n_estimators=40,
        max_depth=2,
        learning_rate=0.1,
        random_state=int(seed),
    )
    q_m.fit(Xr, Yr)
    mse_ref = probe_mse(mean_m, Xh, Yh)
    pin_ref = pinball(Yh, q_m.predict(Xh), tau)

    mse, pin = [], []
    hop_m, hop_q = [], []
    i = n_ref
    while i < len(X):
        Xb, Yb = X[i : i + bs], Y[i : i + bs]
        if len(Xb) < 8:
            break
        m_t = probe_mse(mean_m, Xb, Yb)
        q_t = pinball(Yb, q_m.predict(Xb), tau)
        hop_m.append(bool(large_deviation(m_t, mse_ref, ratio=float(gate))))
        hop_q.append(bool(large_deviation(q_t, pin_ref, ratio=float(gate))))
        mse.append(m_t)
        pin.append(q_t)
        i += bs

    def _first(flags):
        for t, v in enumerate(flags):
            if v:
                return int(t)
        return -1

    return {
        "tau": float(tau),
        "mse_ref": float(mse_ref),
        "pinball_ref": float(pin_ref),
        "MSE_list": np.asarray(mse, dtype=float),
        "pinball_list": np.asarray(pin, dtype=float),
        "hop_mse": np.asarray(hop_m, dtype=bool),
        "hop_pinball": np.asarray(hop_q, dtype=bool),
        "mse_1": _first(hop_m),
        "pinball_1": _first(hop_q),
        "y_in_x": False,
    }


def make_reject_df(n=480, p=6, kind="quiet", onset=240, seed=0):
    """Y* always exists; last column is NaN when S=0.

    quiet     : frozen policy
    select    : after onset, S tighter on X0 (policy). Y*|X frozen
    concept   : after onset, Y* += c X1. S frozen
    covariate : after onset, X0 mean walks. S and Y*|X frozen
    """
    rng = np.random.default_rng(int(seed))
    n = int(n)
    X = rng.normal(size=(n, int(p)))
    w = np.zeros(int(p))
    w[: min(3, int(p))] = np.array([1.2, -0.8, 0.6][: min(3, int(p))])
    Ystar = X @ w + rng.normal(0.0, 0.35, size=n)
    t = np.arange(n)
    onset = int(onset)
    if kind == "covariate":
        X[t >= onset, 0] += 2.2
        Ystar = X @ w + rng.normal(0.0, 0.35, size=n)
    if kind == "concept":
        Ystar[t >= onset] = Ystar[t >= onset] + 1.6 * X[t >= onset, 1]
    b0, b1 = -0.2, 1.1
    logit = b0 + b1 * X[:, 0]
    S = (rng.uniform(size=n) < 1.0 / (1.0 + np.exp(-logit))).astype(float)
    if kind == "select":
        # cream-skimming: only the right tail of X0 is labeled
        S = S.copy()
        S[t >= onset] = (X[t >= onset, 0] > 1.05).astype(float)
    Y = np.where(S > 0.5, Ystar, np.nan)
    return np.column_stack([X, Y])


def make_tail_df(n=480, p=6, onset=240, seed=0):
    """Complete Y. After onset, upper tail on X0>0 moves. Mean can stay quieter."""
    rng = np.random.default_rng(int(seed))
    n = int(n)
    X = rng.normal(size=(n, int(p)))
    w = np.zeros(int(p))
    w[0] = 0.8
    Y = X @ w + rng.normal(0.0, 0.4, size=n)
    tail = (np.arange(n) >= int(onset)) & (X[:, 0] > 0.5)
    Y = Y + 3.2 * tail.astype(float) * np.abs(rng.normal(size=n))
    return np.column_stack([X, Y])
