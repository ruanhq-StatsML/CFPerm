"""Attribution as next-batch training adapter (Amazon rating MSE).

FSDS / TSS already estimates the hop (ĉ from the batch-relationship heatmap,
δ̂ from excess risk after mean-align). This module turns that pair into the
*training weights* for the next batch — not a new loss.

Existing pieces only:
  - heatmap hop cosine → sample weights for closed-form Ridge (domain proximity IW)
  - modality-π hop: w ∝ exp(γ Σ_m π_m (cos(μ_s^m, μ_{t-1}^m)−1)) for multimodal
  - δ̂ / ĉ → gate between instance bank (Wu NW) and hop-weighted Ridge (TSS signs)
  - rolling EWMA of per-predictor MSE → online stacking weight
  - river.PARegressor as an off-the-shelf online baseline

Locked recipe:
  Amazon (scalar text)  → hop_ridge (uniform heatmap IW) — best MSE so far
  Multimodal (MSR-VTT)  → modality-π hop; π_m via ``benchmark_feature_selection``
  PO-risk VIMP          → mask drill-down attribution, not feature reweighting

Locked signs (unchanged):
  ĉ large  → prefer bank / shrink SGD step (do not chase P(X|W))
  δ̂ large  → prefer Ridge refit / raise step (relearn P(Y|X))
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge

from amazon_continuous_batches import (
    CACHE_DIR,
    CATEGORIES,
    ETA0,
    batch_mean_cosine,
    concept_intensity_mse,
    featurize_reviews,
    heatmap_hop_c,
    load_amazon_reviews,
    make_amazon_like_stream,
    run_amazon_method,
)
from amazon_mse_prototype import _mse, _nw_predict, describe_shapes, run_bank_mse
from msrvtt_multimodal_attribution import (
    SEED,
    benchmark_feature_selection,
    modality_mass,
)

METHODS = (
    "plateau",
    "bank",
    "ridge_past",
    "hop_ridge",
    "dga_ridge",
    "attr_adapter",
    "river_pa",
)


def load_amazon_stream(n_per=240, n_batches=9, seed=SEED, max_features=128):
    raw = load_amazon_reviews(
        categories=CATEGORIES[: int(n_batches)],
        n_per=n_per,
        seed=seed,
        cache_dir=CACHE_DIR,
    )
    return featurize_reviews(raw, max_features=max_features)


def rolling_hop_stats(stream):
    """Online-causal (ĉ_t, δ̂_t) on each hop t-1 → t, plus batch means."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    R = batch_mean_cosine(X, batch, n_batches=k)
    mus = np.stack([X[batch == t].mean(axis=0) for t in range(k)])
    c = np.zeros(k)
    d = np.zeros(k)
    for t in range(1, k):
        ip, ic = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        c[t] = heatmap_hop_c(R, t)
        d[t], _ = concept_intensity_mse(X[ip], y[ip], X[ic], y[ic])
    return {"R": R, "mus": mus, "c": c, "delta": d, "n_batches": k}


def hop_sample_weights(mus, batch, t, gamma=4.0):
    """Importance weights for past rows: exp(γ (cos(μ_s, μ_{t-1}) − 1)).

    Causal: the arriving batch mean is unknown, so the reference is μ_{t-1}
    (last observed category). Same cosine board as the FSDS heatmap.
    """
    batch = np.asarray(batch, dtype=int)
    w = np.zeros(batch.shape[0], dtype=float)
    ref = mus[t - 1]
    ref = ref / (np.linalg.norm(ref) + 1e-12)
    for s in range(t):
        mu = mus[s] / (np.linalg.norm(mus[s]) + 1e-12)
        cos = float(np.dot(mu, ref))
        w[batch == s] = float(np.exp(float(gamma) * (cos - 1.0)))
    return np.maximum(w, 1e-3)


def block_means(X, batch, groups):
    """Per-batch, per-modality mean vectors for modality-π hop."""
    batch = np.asarray(batch, dtype=int)
    k = int(batch.max()) + 1
    out = {g: np.zeros((k, sl.stop - sl.start)) for g, sl in groups.items()}
    for t in range(k):
        idx = batch == t
        if not np.any(idx):
            continue
        Xt = X[idx]
        for g, sl in groups.items():
            out[g][t] = Xt[:, sl].mean(axis=0)
    return out


def modality_hop_weights(block_mus, batch, t, shares, gamma=4.0):
    """Modality-specific hop IW: exp(γ Σ_m π_m (cos(μ_s^m, μ_{t-1}^m) − 1)).

    π_m is the modality contribution (from ``modality_pi_shares`` /
    ``benchmark_feature_selection``). Causal reference is still μ_{t-1}.
    With one modality this equals ``hop_sample_weights``.
    """
    batch = np.asarray(batch, dtype=int)
    w = np.zeros(batch.shape[0], dtype=float)
    names = [g for g in shares if g in block_mus]
    if not names:
        return np.ones(batch.shape[0], dtype=float)
    pi = np.asarray([max(float(shares[g]), 1e-6) for g in names], dtype=float)
    pi = pi / pi.sum()
    for s in range(t):
        score = 0.0
        for j, g in enumerate(names):
            a = block_mus[g][s]
            b = block_mus[g][t - 1]
            na = np.linalg.norm(a) + 1e-12
            nb = np.linalg.norm(b) + 1e-12
            score += float(pi[j]) * float(np.dot(a, b) / (na * nb))
        w[batch == s] = float(np.exp(float(gamma) * (score - 1.0)))
    return np.maximum(w, 1e-3)


def modality_pi_shares(
    X0,
    X1,
    seed=SEED,
    n_estimators=40,
    selector=None,
    prev_shares=None,
    ewma=0.0,
):
    """π_m for modality-π hop via ``benchmark_feature_selection``.

    Drop-in: pass another ``selector(X0, X1, seed=, n_estimators=) -> (vimp, meta)``
    (default is RF-Domain VIMP). Shares come from ``modality_mass(vimp)``.
    Optional EWMA with ``prev_shares`` (``ewma`` in [0,1], weight on previous).
    """
    sel = selector or benchmark_feature_selection
    vimp, meta = sel(X0, X1, seed=seed, n_estimators=n_estimators)
    _, shares = modality_mass(vimp)
    if prev_shares is not None and float(ewma) > 0:
        a = float(np.clip(ewma, 0.0, 1.0))
        shares = {
            g: (1.0 - a) * float(shares.get(g, 0.0)) + a * float(prev_shares.get(g, 0.0))
            for g in shares
        }
        tot = sum(shares.values()) + 1e-12
        shares = {g: shares[g] / tot for g in shares}
    return shares, vimp, meta


def _ridge_predict(Xtr, ytr, Xte, sample_weight=None, alpha=3.0):
    clf = Ridge(alpha=float(alpha))
    if sample_weight is None:
        clf.fit(Xtr, ytr)
    else:
        clf.fit(Xtr, ytr, sample_weight=sample_weight)
    return np.clip(clf.predict(Xte), 1.0, 5.0)


def _ridge_fit(Xtr, ytr, sample_weight=None, alpha=3.0):
    clf = Ridge(alpha=float(alpha))
    if sample_weight is None:
        clf.fit(Xtr, ytr)
    else:
        clf.fit(Xtr, ytr, sample_weight=sample_weight)
    return clf


def _linear_mse_grad(clf, X, y):
    """∇_θ of mean squared error for a fitted linear model (coef ‖ intercept).

    Used by DGA (Fan, Grangier, Ablin 2024): domain weight ∝ gradient alignment
    with the specialized set, not a density ratio.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    n = max(int(X.shape[0]), 1)
    pred = clf.predict(X)
    resid = pred - y
    g_coef = (X.T @ resid) / float(n)
    g_int = float(np.mean(resid))
    return np.concatenate([np.asarray(g_coef, dtype=float).ravel(), [g_int]])


def dga_alignments(clf, X, y, batch, domains, spe_domain, align="dot"):
    """Per-domain gradient alignment a_i = ⟨∇ℓ(θ, D_i), ∇ℓ(θ, D_spe)⟩.

    ``align='cosine'`` uses cosine of the two gradients (scale-stable surrogate
    still in the DGA/DoGE family). ``align='dot'`` is the paper inner product.
    """
    batch = np.asarray(batch, dtype=int)
    is_spe = batch == int(spe_domain)
    if not np.any(is_spe):
        return np.zeros(len(domains), dtype=float)
    g_spe = _linear_mse_grad(clf, X[is_spe], y[is_spe])
    n_spe = float(np.linalg.norm(g_spe) + 1e-12)
    out = np.zeros(len(domains), dtype=float)
    for j, s in enumerate(domains):
        idx = batch == int(s)
        if not np.any(idx):
            continue
        g = _linear_mse_grad(clf, X[idx], y[idx])
        if align == "cosine":
            out[j] = float(np.dot(g, g_spe) / ((np.linalg.norm(g) + 1e-12) * n_spe))
        else:
            out[j] = float(np.dot(g, g_spe))
    return out


def dga_mirror_step(alpha, alignments, eta=1.0):
    """Simplex mirror descent: α ← normalize(α ⊙ exp(η a)).

    Upweights domains whose gradient *aligns* with D_spe (Fan et al. 2024 §2.4
    Taylor argument). Note: their Algorithm 1 prints ``exp(-η a)``; that sign
    contradicts the surrounding derivation that maximizes alignment — we follow
    the derivation / Eq. for increasing α on large ⟨∇L_i, ∇L_spe⟩.
    """
    a = np.asarray(alignments, dtype=float).copy()
    alpha = np.asarray(alpha, dtype=float)
    if alpha.shape != a.shape:
        raise ValueError("alpha/alignments shape mismatch")
    if a.size == 0:
        return alpha
    # center + scale so η stays O(1); single-domain / flat a → no-op
    a = a - float(np.mean(a))
    scale = float(np.max(np.abs(a)))
    if scale < 1e-12:
        return alpha / alpha.sum()
    a = a / scale
    # clip exponent to avoid overflow when η is large
    log_hat = np.log(np.maximum(alpha, 1e-12)) + float(eta) * a
    log_hat = log_hat - float(np.max(log_hat))
    hat = np.exp(np.clip(log_hat, -60.0, 0.0))
    hat = np.maximum(hat, 1e-12)
    return hat / hat.sum()


def dga_sample_weights(alpha, batch, domains):
    """Broadcast domain simplex weights onto rows."""
    batch = np.asarray(batch, dtype=int)
    w = np.ones(batch.shape[0], dtype=float)
    for j, s in enumerate(domains):
        w[batch == int(s)] = float(alpha[j])
    return np.maximum(w, 1e-6)


def run_dga_ridge(
    stream,
    alpha=3.0,
    eta=1.0,
    ema_beta=0.35,
    align="cosine",
    ridge_alpha=None,
):
    """DGA domain reweighting + closed-form Ridge (Fan, Grangier, Ablin 2024).

    Mapping onto the Amazon hop board (existing method, not a new loss):
      generic domains D_i  = past rating categories / batches
      specialized set D_spe = last observed batch (causal proxy for the next hop)
      α_i via gradient alignment + mirror descent + EMA
      train = sample-weighted Ridge with w_row = α_{batch(row)}

    Why this instead of density-ratio IW: DGA never estimates p_te/p_tr. It
    upweights domains whose *training gradient* currently aligns with the
    specialized set — the same reason IS/embedding hop weights fail when the
    specialized pocket is tiny or a domain is overfit (Fan et al. §1, §3.1).

    Streaming note: the paper uses a *fixed* domain simplex. Here the number of
    past batches grows, so each hop takes one mirror step from the uniform
    prior on the current simplex (≡ softmax of tempered alignments). EMA
    carries temporal memory of α without freezing mass on early batches.
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    ra = float(alpha if ridge_alpha is None else ridge_alpha)
    history = []
    alpha_ema = None
    for t in range(1, k):
        tr = batch < t
        ic = batch == t
        domains = list(range(t))
        spe = t - 1
        if alpha_ema is None:
            probe_w = None
        else:
            prev = np.asarray(alpha_ema, dtype=float)
            if prev.size == t:
                mix = prev
            elif prev.size == t - 1:
                # new domain: give it the mean mass of existing domains
                mean_m = float(prev.mean()) if prev.size else 1.0
                mix = np.concatenate([prev, [mean_m]])
                mix = mix / mix.sum()
            else:
                mix = np.ones(t, dtype=float) / float(t)
            probe_w = dga_sample_weights(mix, batch, domains)[tr]
        clf = _ridge_fit(X[tr], y[tr], sample_weight=probe_w, alpha=ra)
        a = dga_alignments(clf, X[tr], y[tr], batch[tr], domains, spe, align=align)
        # one MD step from uniform on the *current* simplex (avoids sticky early mass)
        alpha_inst = dga_mirror_step(np.ones(t, dtype=float) / float(t), a, eta=eta)
        if alpha_ema is None:
            alpha_ema = alpha_inst.copy()
        else:
            b = float(np.clip(ema_beta, 0.0, 1.0))
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
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary(
        "dga_ridge",
        online,
        history,
        stream,
        extras={
            "eta": float(eta),
            "ema_beta": float(ema_beta),
            "align": align,
            "alpha": ra,
            "method_ref": "Fan, Grangier, Ablin 2024 — Dynamic Gradient Alignment (arXiv:2410.02498)",
        },
    )


def run_ridge_past(stream, alpha=3.0):
    """Closed-form Ridge on all past batches (uniform weights)."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k):
        tr, ic = batch < t, batch == t
        pred = _ridge_predict(X[tr], y[tr], X[ic], alpha=alpha)
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "n": int(ic.sum()),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary("ridge_past", online, history, stream)


def run_hop_ridge(stream, alpha=3.0, gamma=4.0):
    """Heatmap-hop importance-weighted Ridge — attribution as sample weights."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    stats = rolling_hop_stats(stream)
    history = []
    for t in range(1, stats["n_batches"]):
        tr, ic = batch < t, batch == t
        w = hop_sample_weights(stats["mus"], batch, t, gamma=gamma)
        pred = _ridge_predict(X[tr], y[tr], X[ic], sample_weight=w[tr], alpha=alpha)
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "c": float(stats["c"][t]),
                "delta": float(stats["delta"][t]),
                "w_mean": float(w[tr].mean()),
                "w_max": float(w[tr].max()),
                "n": int(ic.sum()),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary("hop_ridge", online, history, stream, extras={"gamma": gamma, "alpha": alpha})


def typed_bank_weight(c, delta, ewma_bank=None, ewma_ridge=None):
    """TSS signs → mixture weight on the instance bank.

    ĉ large / δ̂ quiet → bank up. δ̂ large → bank down (Ridge refit).
    Optional EWMA stack blends with inverse recent error.
    """
    c = float(c)
    delta = float(delta)
    w_cd = (1.0 + 2.0 * c) / (1.0 + 2.0 * c + 4.0 * delta)
    if ewma_bank is None or ewma_ridge is None:
        return float(np.clip(w_cd, 0.05, 0.70))
    inv_b = 1.0 / (float(ewma_bank) + 1e-3)
    inv_r = 1.0 / (float(ewma_ridge) + 1e-3)
    w_err = inv_b / (inv_b + inv_r)
    return float(np.clip(0.5 * w_cd + 0.5 * w_err, 0.05, 0.70))


def attribution_ewma_rate(c, delta, pi_change=0.0):
    """How much *new* hop MSE enters the ensemble EWMA (λ_new ∈ [0.15, 0.75]).

    Attribution guides next-batch memory:
      ĉ large / δ̂ quiet (covariate hop) → small λ_new → keep past ensemble
      δ̂ large (concept hop)             → large λ_new → forget fast, track Ridge
      π_m jumps                          → bump λ_new (modality mix moved)

    Locked TSS signs only — no new loss. ``pi_change`` = L1(|π_t − π_{t-1}|)/2 ∈ [0,1].
    """
    c = float(c)
    d = float(delta)
    react = (1.0 + 4.0 * d) / (1.0 + 2.0 * c + 4.0 * d)
    react = react + 0.5 * float(np.clip(pi_change, 0.0, 1.0))
    return float(np.clip(0.15 + 0.60 * react, 0.15, 0.75))


def pi_l1_change(prev, curr):
    """Half L1 distance between modality share dicts (∈ [0,1])."""
    if not prev or not curr:
        return 0.0
    keys = set(prev) | set(curr)
    return 0.5 * float(sum(abs(float(curr.get(g, 0.0)) - float(prev.get(g, 0.0))) for g in keys))


def run_attr_adapter(stream, alpha=3.0, gamma=4.0, tau=0.20, cap=800):
    """Next-batch adapter: typed mix of Wu bank and hop-weighted Ridge.

    Predict the arriving batch with
        ŷ = w · bank_NW + (1−w) · Ridge_hop(past),
    w from (ĉ, δ̂) and rolling predictor EWMA. The EWMA *rate* λ_new is set by
    ``attribution_ewma_rate(ĉ, δ̂)`` so attribution guides next-hop ensemble
    memory (concept shift → react; covariate shift → remember).
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    stats = rolling_hop_stats(stream)
    i0 = np.flatnonzero(batch == 0)
    Z, yz = X[i0].copy(), y[i0].copy()
    ewma_bank = ewma_ridge = None
    history = []
    for t in range(1, stats["n_batches"]):
        ic = np.flatnonzero(batch == t)
        tr = batch < t
        Xt, yt = X[ic], y[ic]
        c, d = float(stats["c"][t]), float(stats["delta"][t])
        w = typed_bank_weight(c, d, ewma_bank, ewma_ridge)
        pb = _nw_predict(Xt, Z, yz, tau=tau)
        sw = hop_sample_weights(stats["mus"], batch, t, gamma=gamma)
        pr = _ridge_predict(X[tr], y[tr], Xt, sample_weight=sw[tr], alpha=alpha)
        pred = w * pb + (1.0 - w) * pr
        mb, mr = _mse(pb, yt), _mse(pr, yt)
        # attribution → next-hop ensemble memory
        lam = attribution_ewma_rate(c, d)
        ewma_bank = mb if ewma_bank is None else (1.0 - lam) * ewma_bank + lam * mb
        ewma_ridge = mr if ewma_ridge is None else (1.0 - lam) * ewma_ridge + lam * mr
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, yt),
                "mse_bank": mb,
                "mse_ridge": mr,
                "w_bank": w,
                "ewma_rate": lam,
                "c": c,
                "delta": d,
                "n": int(ic.size),
            }
        )
        Z = np.vstack([Z, Xt])
        yz = np.concatenate([yz, yt])
        if len(yz) > int(cap):
            Z, yz = Z[-int(cap) :], yz[-int(cap) :]
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary(
        "attr_adapter",
        online,
        history,
        stream,
        extras={
            "alpha": alpha,
            "gamma": gamma,
            "mean_w_bank": float(np.mean([h["w_bank"] for h in history])),
            "mean_ewma_rate": float(np.mean([h["ewma_rate"] for h in history])),
        },
    )


def run_river_pa(stream, C=0.01):
    """river PassiveAwareRegressor — package baseline on the same stream."""
    from river import compose, linear_model, preprocessing

    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    p = X.shape[1]
    model = compose.Pipeline(
        preprocessing.StandardScaler(),
        linear_model.PARegressor(C=float(C), mode=2),
    )

    def row(i):
        return {j: float(X[i, j]) for j in range(p)}

    for i in np.flatnonzero(batch == 0):
        model.learn_one(row(i), float(y[i]))
    history = []
    for t in range(1, k):
        ic = np.flatnonzero(batch == t)
        preds = []
        for i in ic:
            yhat = model.predict_one(row(i))
            preds.append(float(np.clip(3.0 if yhat is None else yhat, 1.0, 5.0)))
        history.append({"round": int(t), "online_mse": _mse(preds, y[ic]), "n": int(ic.size)})
        for i in ic:
            model.learn_one(row(i), float(y[i]))
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary("river_pa", online, history, stream, extras={"C": C, "package": "river"})


def _summary(method, online, history, stream, extras=None):
    out = {
        "method": method,
        "n_batches": int(np.asarray(stream.batch).max()) + 1,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "cum_mse": float(online.sum()) if online.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "online_path": online.tolist(),
        "history": history,
        "meta": dict(getattr(stream, "meta", {}) or {}),
    }
    if extras:
        out.update(extras)
    return out


def run_adapter_method(stream, method, seed=SEED, **kw):
    if method == "plateau":
        rec, _ = run_amazon_method(stream, method="plateau", seed=seed, **kw)
        return {
            "method": "plateau",
            "n_batches": rec["n_batches"],
            "online_mse": rec["online_mse"],
            "cum_mse": rec["cum_mse"],
            "last_online_mse": rec["last_online_mse"],
            "online_path": rec.get("online_path", []),
            "history": rec.get("history", []),
            "meta": dict(getattr(stream, "meta", {}) or {}),
        }
    if method == "bank":
        return run_bank_mse(stream)
    if method == "ridge_past":
        return run_ridge_past(stream)
    if method == "hop_ridge":
        return run_hop_ridge(stream)
    if method == "dga_ridge":
        return run_dga_ridge(stream)
    if method == "attr_adapter":
        return run_attr_adapter(stream)
    if method == "river_pa":
        return run_river_pa(stream)
    raise ValueError("unknown method %s" % method)


def run_adapter_suite(
    seeds=None,
    n_per=240,
    n_batches=9,
    source="amazon",
    methods=None,
    steps_per_batch=8,
):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    methods = list(methods or METHODS)
    rows = []
    traces = {m: [] for m in methods}
    shapes = None
    hop_logs = []
    for s in seeds:
        if source == "synthetic":
            stream = make_amazon_like_stream(
                n_batches=n_batches, n_per=n_per, seed=int(s), cov=0.35
            )
        else:
            stream = load_amazon_stream(n_per=n_per, n_batches=n_batches, seed=int(s))
        if shapes is None:
            shapes = describe_shapes(stream)
        hops = rolling_hop_stats(stream)
        hop_logs.append(
            {
                "seed": int(s),
                "c": hops["c"].tolist(),
                "delta": hops["delta"].tolist(),
            }
        )
        for method in methods:
            kw = {}
            if method == "plateau":
                kw["steps_per_batch"] = steps_per_batch
            rec = run_adapter_method(stream, method, seed=int(s), **kw)
            traces[method].append(rec)
            rows.append(
                {
                    "seed": int(s),
                    "method": method,
                    "online_mse": rec["online_mse"],
                    "cum_mse": rec["cum_mse"],
                    "last_online_mse": rec["last_online_mse"],
                }
            )
    table = {}
    for method in methods:
        v = np.array([r["online_mse"] for r in rows if r["method"] == method], dtype=float)
        table[method] = {
            "online_mse": {
                "mean": float(v.mean()),
                "sd": float(v.std(ddof=1) if len(v) > 1 else 0.0),
            },
            "n": int(len(v)),
        }
    # paired lift vs bank / plateau
    lifts = {}
    for baseline in ("bank", "plateau"):
        if baseline not in methods:
            continue
        lifts[baseline] = {}
        for method in methods:
            if method == baseline:
                continue
            diffs = []
            for s in seeds:
                a = next(r["online_mse"] for r in rows if r["seed"] == s and r["method"] == method)
                b = next(r["online_mse"] for r in rows if r["seed"] == s and r["method"] == baseline)
                diffs.append(b - a)  # positive => method better (lower MSE)
            d = np.asarray(diffs, dtype=float)
            lifts[baseline][method] = {
                "mean_mse_drop": float(d.mean()),
                "sd": float(d.std(ddof=1) if len(d) > 1 else 0.0),
            }
    return {
        "table": table,
        "rows": rows,
        "traces": traces,
        "shapes": shapes,
        "lifts": lifts,
        "hop_logs": hop_logs,
        "seeds": [int(s) for s in seeds],
        "source": source,
        "methods": methods,
    }


def plot_adapter_suite(suite, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import GRID, INK, MUTED, _save, _style

    _style()
    colors = {
        "plateau": "#C4892A",
        "bank": "#2F6B4F",
        "ridge_past": "#6B7C8A",
        "hop_ridge": "#C45C26",
        "dga_ridge": "#5B3A8C",
        "attr_adapter": "#2C4A6E",
        "river_pa": "#9AA3AE",
    }
    labels = {
        "plateau": "plateau RidgeSGD",
        "bank": "instance bank",
        "ridge_past": "Ridge on past",
        "hop_ridge": "hop-weighted Ridge",
        "dga_ridge": "DGA-weighted Ridge",
        "attr_adapter": "attr adapter (bank⊕hop)",
        "river_pa": "river PA",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), gridspec_kw={"width_ratios": [1.15, 1.0]})
    ax = axes[0]
    methods = [m for m in suite["methods"] if m in suite["traces"]]
    for method in methods:
        paths = np.array([r["online_path"] for r in suite["traces"][method]], dtype=float)
        if paths.size == 0 or paths.shape[1] == 0:
            continue
        t = np.arange(1, paths.shape[1] + 1)
        mu, sd = paths.mean(0), paths.std(0)
        ax.plot(t, mu, color=colors.get(method, INK), lw=2.1, label=labels.get(method, method))
        ax.fill_between(t, mu - sd, mu + sd, color=colors.get(method, INK), alpha=0.12, lw=0)
    ax.set_xlabel("batch")
    ax.set_ylabel("online MSE")
    ax.set_title("Attribution adapter on Amazon ratings", loc="left", fontsize=12.5, fontweight="bold")
    ax.grid(True, color=GRID)
    ax.legend(frameon=False, fontsize=7.8)

    ax = axes[1]
    means = [suite["table"][m]["online_mse"]["mean"] for m in methods]
    sds = [suite["table"][m]["online_mse"]["sd"] for m in methods]
    ypos = np.arange(len(methods))[::-1]
    ax.barh(
        ypos,
        means,
        xerr=sds,
        color=[colors.get(m, INK) for m in methods],
        height=0.62,
        error_kw={"ecolor": "#555", "lw": 1.0, "capsize": 2},
    )
    ax.set_yticks(ypos)
    ax.set_yticklabels([labels.get(m, m) for m in methods], fontsize=8.2)
    ax.set_xlabel("mean online MSE")
    ax.set_title("Lower is better", loc="left", fontsize=12.5, fontweight="bold")
    ax.grid(True, color=GRID, axis="x")
    fig.text(
        0.04,
        0.01,
        "hop-weighted Ridge uses heatmap cos(μ_s, μ_{t−1}) as sample weights. "
        "DGA-weighted Ridge uses Fan et al. 2024 gradient alignment (+EMA) instead of density-ratio IW. "
        "attr adapter mixes Wu bank and hop Ridge with TSS signs on (ĉ, δ̂).",
        fontsize=8.0,
        color=MUTED,
    )
    return _save(fig, Path(path))


def write_adapter_tex(suite, path):
    methods = suite["methods"]
    labels = {
        "plateau": "plateau RidgeSGD",
        "bank": "instance bank (Wu NW)",
        "ridge_past": "Ridge on all past",
        "hop_ridge": "hop-weighted Ridge",
        "dga_ridge": "DGA-weighted Ridge",
        "attr_adapter": "attr.\\ adapter (bank$\\oplus$hop)",
        "river_pa": "river PARegressor",
    }
    lines = [
        r"% Attribution adapter on Amazon. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Multimodal / batch attribution as a next-batch training adapter on Amazon Reviews 2023.",
        r"Online rating MSE (predict the arriving category, then update). Hop-weighted Ridge reuses the",
        r"batch-relationship heatmap as sample weights $w_s=\exp(\gamma(\cos(\bar x_s,\bar x_{t-1})-1))$.",
        r"DGA-weighted Ridge (Fan--Grangier--Ablin 2024) replaces that density-ratio surrogate by",
        r"mirror-descent domain weights from gradient alignment with $D_{\mathrm{spe}}=B_{t-1}$.",
        r"The typed adapter mixes the Wu instance bank with hop Ridge using the locked TSS signs on $(\hat c,\hat\delta)$.}",
        r"\label{tab:attr-adapter-amazon}",
        r"\small",
        r"\begin{tabular}{@{}lcc@{}}\toprule",
        r"Method & online MSE & vs bank \\",
        r"\midrule",
    ]
    bank_mu = suite["table"].get("bank", {}).get("online_mse", {}).get("mean", float("nan"))
    for m in methods:
        cell = suite["table"][m]["online_mse"]
        drop = bank_mu - cell["mean"] if np.isfinite(bank_mu) else float("nan")
        lines.append(
            r"%s & $%.3f$ ($%.3f$) & $%+.3f$ \\"
            % (labels.get(m, m), cell["mean"], cell["sd"], drop)
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path
