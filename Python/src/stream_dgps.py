"""Gradual concept-drift vs gradual covariate-shift DGPs.

Concept: P(X) fixed, P(Y|X) rotates after a labeled onset batch.
Covariate: P(Y|X) fixed, mean of X walks after the same onset.

Labeled onset is a stream-batch index (0-based, after D_ref).
OnlineRFPerm marks WHEN. PO × MSE × MMD vs X_ref marks WHAT.
"""
from __future__ import annotations

import numpy as np

ONSET_BATCH = 2


def _sigmoid(z):
    z = np.clip(np.asarray(z, dtype=float), -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(-z))


def _alpha(n_ref: int, n_new: int, n_batches: int, onset_batch: int) -> np.ndarray:
    n = int(n_ref) + int(n_new) * int(n_batches)
    idx = np.arange(n)
    start = int(n_ref) + int(onset_batch) * int(n_new)
    span = max(n - start - 1, 1)
    a = np.clip((idx - start) / span, 0.0, 1.0)
    a[idx < start] = 0.0
    return a


def make_gradual_concept(
    n_ref: int = 10_000,
    n_new: int = 2_000,
    n_batches: int = 12,
    onset_batch: int = ONSET_BATCH,
    p: int = 8,
    seed: int = 2026,
):
    """P(X)~N(0,I). β rotates from β0 to −β0 after onset_batch."""
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    X = rng.normal(size=(n, p))
    beta0 = np.zeros(p)
    beta0[:4] = 1.2
    beta1 = -beta0
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    logit = (X @ beta0) * (1.0 - a) + (X @ beta1) * a
    Y = rng.binomial(1, _sigmoid(logit)).astype(float)
    meta = {
        "kind": "concept",
        "onset_batch": int(onset_batch),
        "title": f"DGP gradual concept (β flips after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, meta


def make_gradual_covariate(
    n_ref: int = 10_000,
    n_new: int = 2_000,
    n_batches: int = 12,
    onset_batch: int = ONSET_BATCH,
    p: int = 8,
    shift: float = 2.5,
    seed: int = 2026,
):
    """P(Y|X) is a radial bump on the first 4 coords. Mean of those coords walks.

    Same f before and after. After onset the mass leaves the bump, so a model
    fit on D_ref pays MSE while PO-risk (same f) stays quieter than concept.
    """
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    X0 = rng.normal(size=(n, p))
    delta = np.zeros(p)
    delta[:4] = float(shift)
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    X = X0 + a[:, None] * delta
    r = np.sqrt(np.sum(X[:, :4] ** 2, axis=1))
    Y = rng.binomial(1, _sigmoid(2.2 - 1.6 * r)).astype(float)
    meta = {
        "kind": "covariate",
        "onset_batch": int(onset_batch),
        "title": f"DGP gradual covariate (bump; μ walks after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, meta


TRIMODAL_GROUPS = {
    "video": slice(0, 4),
    "audio": slice(4, 8),
    "text": slice(8, 12),
}


def make_trimodal_stream(
    n_ref: int = 800,
    n_new: int = 200,
    n_batches: int = 6,
    onset_batch: int = ONSET_BATCH,
    kind: str = "concept_video",
    n_slices: int = 4,
    shift: float = 2.2,
    seed: int = 2026,
):
    """Three feature blocks standing in for video / audio / text.

    Y depends on all three blocks. After onset:
      concept_video — flip the video coefficients; P(X) fixed
      covariate_audio — walk the audio mean; same f
      both — both of the above
    Only slice == (n_slices-1) is shifted, so subset tables have a planted peak.
    """
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    p = 12
    X = rng.normal(size=(n, p))
    n_slices = max(2, int(n_slices))
    slice_id = rng.integers(0, n_slices, size=n)
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    hot = (slice_id == (n_slices - 1)).astype(float)
    w = a * hot

    kind = str(kind)
    if kind in ("covariate_audio", "both"):
        X[:, 4:8] = X[:, 4:8] + w[:, None] * float(shift)

    beta_v = np.array([1.2, 0.4, 0.0, 0.0])
    beta_a = np.array([0.45, 0.0, 0.0, 0.0])
    beta_t = np.array([0.45, 0.0, 0.0, 0.0])
    logit0 = X[:, 0:4] @ beta_v + X[:, 4:8] @ beta_a + X[:, 8:12] @ beta_t
    if kind in ("concept_video", "both"):
        logit = logit0 * (1.0 - w) + (
            X[:, 0:4] @ (-beta_v) + X[:, 4:8] @ beta_a + X[:, 8:12] @ beta_t
        ) * w
    else:
        logit = logit0
    Y = rng.binomial(1, _sigmoid(logit)).astype(float)
    meta = {
        "kind": kind,
        "onset_batch": int(onset_batch),
        "groups": {k: list(range(s.start, s.stop)) for k, s in TRIMODAL_GROUPS.items()},
        "shifted_slice": int(n_slices - 1),
        "n_slices": int(n_slices),
        "title": f"trimodal {kind} (slice {n_slices - 1} after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, slice_id, meta


def make_trimodal_gradual_concept(
    n_ref: int = 800,
    n_new: int = 40,
    n_batches: int = 16,
    onset_batch: int = 3,
    seed: int = 2026,
):
    """Continuous-time gradual concept: P(X) fixed, video β rotates for every row.

    Clock is the observation index τ. α(τ) walks 0→1 after labeled onset.
    There is no jump. Last-two hops on a slow walk should stay quiet.
    """
    rng = np.random.default_rng(seed)
    n_ref = int(n_ref)
    n_new = int(n_new)
    n_batches = int(n_batches)
    n = n_ref + n_new * n_batches
    X = rng.normal(size=(n, 12))
    alpha = _alpha(n_ref, n_new, n_batches, onset_batch)
    tau = np.arange(n, dtype=float)
    beta_v = np.array([1.2, 0.4, 0.0, 0.0])
    beta_a = np.array([0.45, 0.0, 0.0, 0.0])
    beta_t = np.array([0.45, 0.0, 0.0, 0.0])
    logit0 = X[:, 0:4] @ beta_v + X[:, 4:8] @ beta_a + X[:, 8:12] @ beta_t
    logit1 = X[:, 0:4] @ (-beta_v) + X[:, 4:8] @ beta_a + X[:, 8:12] @ beta_t
    logit = logit0 * (1.0 - alpha) + logit1 * alpha
    Y = rng.binomial(1, _sigmoid(logit)).astype(float)
    onset_tau = n_ref + int(onset_batch) * n_new
    meta = {
        "kind": "gradual_concept_video",
        "onset_batch": int(onset_batch),
        "onset_tau": int(onset_tau),
        "clock": "observation_index",
        "groups": {k: list(range(s.start, s.stop)) for k, s in TRIMODAL_GROUPS.items()},
        "title": f"gradual concept on video (α walks after τ={onset_tau})",
        "n_ref": n_ref,
        "n_new": n_new,
        "n_batches": n_batches,
    }
    return X, Y, alpha, tau, meta


ORDER_FEATS = ("amount", "hour", "n_items", "channel")
MERCHANT_FEATS = ("merchant_cat", "merchant_gmv", "n_skus")
USER_FEATS = ("user_tenure", "user_hist_freq")
GRAIN_FEATS = {
    "order": ORDER_FEATS,
    "merchant": MERCHANT_FEATS,
    "user": USER_FEATS,
}
SOUTH_MERCHANTS = (4, 5, 6, 7)
GMV_EMA_LAM = 0.25

# Native-grain catalog. Y is not a feature. planted_how is filled per kind in meta.
FEATURE_LIBRARY = (
    {"feature": "amount", "grain": "order", "role": "本单金额（进 logit）"},
    {"feature": "hour", "grain": "order", "role": "下单时刻（噪声列）"},
    {"feature": "n_items", "grain": "order", "role": "件数（噪声列）"},
    {"feature": "channel", "grain": "order", "role": "本单渠道"},
    {"feature": "merchant_cat", "grain": "merchant", "role": "商户类目（进 logit，不种 shift）"},
    {
        "feature": "merchant_gmv",
        "grain": "merchant",
        "role": "商户 GMV：上一窗订单金额 EMA，不含本窗、不含 Y、不含未来",
    },
    {"feature": "n_skus", "grain": "merchant", "role": "SKU 数（噪声列）"},
    {"feature": "user_tenure", "grain": "user", "role": "进流前 tenure（进 logit，不种 shift）"},
    {
        "feature": "user_hist_freq",
        "grain": "user",
        "role": "用户历史频次：进流前快照，本窗不更新（不含本单、不含 Y、不含未来）",
    },
)


def causal_ema_lookup(ids, values, n_entities, prior, lam: float = GMV_EMA_LAM):
    """Serve entity state, then update. No current row, no future, no Y."""
    ids = np.asarray(ids, dtype=int)
    values = np.asarray(values, dtype=float).ravel()
    state = np.asarray(prior, dtype=float).ravel().copy()
    if state.shape[0] != int(n_entities):
        raise ValueError("prior length != n_entities")
    if ids.shape[0] != values.shape[0]:
        raise ValueError("ids/values length mismatch")
    lam = float(lam)
    out = np.empty(ids.shape[0], dtype=float)
    for i, e in enumerate(ids):
        out[i] = state[e]
        state[e] = (1.0 - lam) * state[e] + lam * values[i]
    return out


def causal_batch_ema(
    ids, values, n_entities, prior, n_ref, n_new, n_batches, lam: float = GMV_EMA_LAM
):
    """GMV at batch t is the close of t-1. Current batch, Y, and future stay out."""
    ids = np.asarray(ids, dtype=int)
    values = np.asarray(values, dtype=float).ravel()
    n_ref = int(n_ref)
    n_new = int(n_new)
    n_batches = int(n_batches)
    state = np.asarray(prior, dtype=float).ravel().copy()
    out = np.empty(ids.shape[0], dtype=float)
    lam = float(lam)
    for i in range(n_ref):
        e = ids[i]
        out[i] = state[e]
        state[e] = (1.0 - lam) * state[e] + lam * values[i]
    for b in range(n_batches):
        lo = n_ref + b * n_new
        hi = lo + n_new
        snap = state.copy()
        for i in range(lo, hi):
            out[i] = snap[ids[i]]
        for i in range(lo, hi):
            e = ids[i]
            state[e] = (1.0 - lam) * state[e] + lam * values[i]
    return out


def causal_batch_share(ids, n_entities, prior, n_ref, n_new, n_batches, win: int):
    """User share at batch t uses only ids before t. No current batch, no Y."""
    from collections import deque

    ids = np.asarray(ids, dtype=int)
    n_ref = int(n_ref)
    n_new = int(n_new)
    n_batches = int(n_batches)
    win = max(int(win), 1)
    prior = np.asarray(prior, dtype=float).ravel()
    ps = np.exp(prior - float(np.max(prior)))
    ps = ps / ps.sum()
    buf = deque()
    counts = np.zeros(int(n_entities), dtype=float)
    out = np.empty(ids.shape[0], dtype=float)

    def _share(e):
        return counts[e] / len(buf) if buf else float(ps[e])

    def _push(e):
        buf.append(int(e))
        counts[e] += 1.0
        if len(buf) > win:
            old = buf.popleft()
            counts[old] -= 1.0

    for i in range(n_ref):
        e = ids[i]
        out[i] = _share(e)
        _push(e)
    for b in range(n_batches):
        lo = n_ref + b * n_new
        hi = lo + n_new
        snap = np.array([_share(e) for e in range(int(n_entities))], dtype=float)
        for i in range(lo, hi):
            out[i] = snap[ids[i]]
        for i in range(lo, hi):
            _push(ids[i])
    return out


def causal_rolling_share(ids, n_entities, prior, win: int):
    """Serve this entity's share of the last `win` ids, then push. No Y."""
    from collections import deque

    ids = np.asarray(ids, dtype=int)
    prior = np.asarray(prior, dtype=float).ravel()
    if prior.shape[0] != int(n_entities):
        raise ValueError("prior length != n_entities")
    win = max(int(win), 1)
    ps = np.exp(prior - float(np.max(prior)))
    ps = ps / ps.sum()
    buf = deque()
    counts = np.zeros(int(n_entities), dtype=float)
    out = np.empty(ids.shape[0], dtype=float)
    for i, e in enumerate(ids):
        out[i] = counts[e] / len(buf) if buf else float(ps[e])
        buf.append(int(e))
        counts[e] += 1.0
        if len(buf) > win:
            old = buf.popleft()
            counts[old] -= 1.0
    return out


def planted_how_for_kind(kind: str) -> dict[str, str]:
    """How a catalog column is planted. Empty → not planted. Never includes Y."""
    kind = str(kind)
    if kind == "covariate_south":
        return {"amount": "x", "channel": "x", "merchant_gmv": "x"}
    if kind == "concept_south":
        return {"amount": "y|x"}
    if kind == "both":
        return {"amount": "both", "channel": "x", "merchant_gmv": "x"}
    return {}


def south_merchant_ids(n_merchants: int) -> np.ndarray:
    """Second half of merchant ids. n=8 → {4,5,6,7}, same as SOUTH_MERCHANTS."""
    n_merchants = int(n_merchants)
    return np.arange(n_merchants // 2, n_merchants)


def make_order_graph_stream(
    n_ref: int = 480,
    n_new: int = 120,
    n_batches: int = 4,
    onset_batch: int = 1,
    n_merchants: int = 8,
    n_users: int = 40,
    kind: str = "covariate_south",
    shift: float = 2.4,
    seed: int = 2026,
):
    """Order–merchant–user graph. Y is never a feature.

    Edges: order→merchant, order→user. Region is a merchant attribute
    (north = ids 0..3, south = 4..7). After onset, only **south** orders
    are shifted:

      covariate_south — amount and channel walk; f fixed. GMV follows
                        from past amounts (not planted on the row).
      concept_south   — amount coefficient flips; P(X) fixed
      both            — both of the above

    Grain features are computed as they would be at serving:
      order    — this ticket (amount, hour, items, channel)
      merchant — lookup + GMV EMA of past amounts for that store
      user     — pre-period tenure and pre-period frequency (frozen snapshots)

    GMV / hist_freq never see the current row, Y, or future rows.
    Subset key is region / merchant_group, not Y.
    """
    rng = np.random.default_rng(seed)
    n_ref = int(n_ref)
    n_new = int(n_new)
    n_batches = int(n_batches)
    n_merchants = int(n_merchants)
    n_users = int(n_users)
    n = n_ref + n_new * n_batches
    kind = str(kind)

    merchant_id = rng.integers(0, n_merchants, size=n)
    user_id = rng.integers(0, n_users, size=n)
    south_ids = south_merchant_ids(n_merchants)
    south = np.isin(merchant_id, south_ids)
    region = np.where(south, "south", "north")

    merchant_cat = rng.normal(size=n_merchants)
    merchant_gmv0 = rng.normal(size=n_merchants)
    n_skus = rng.normal(size=n_merchants)
    user_tenure = rng.normal(size=n_users)
    user_hist_freq = rng.normal(size=n_users)

    amount0 = rng.normal(size=n) + 0.35 * merchant_cat[merchant_id]
    hour = rng.normal(size=n)
    n_items = rng.normal(size=n)
    channel0 = rng.normal(size=n)

    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    w = a * south.astype(float)
    amount = amount0.copy()
    channel = channel0.copy()
    if kind in ("covariate_south", "both"):
        amount = amount + w * float(shift)
        channel = channel + w * (0.7 * float(shift))

    # Merchant GMV and user frequency are serving lookups: update after serve.
    merchant_gmv = causal_batch_ema(
        merchant_id,
        amount,
        n_merchants,
        merchant_gmv0,
        n_ref,
        n_new,
        n_batches,
        lam=GMV_EMA_LAM,
    )
    # Frequency is a closed-period snapshot. Updating it online would write
    # the clock into the user table (D_ref cold-start vs later batches).
    user_freq = user_hist_freq[user_id]

    X_order = np.column_stack([amount, hour, n_items, channel])
    X_merchant = np.column_stack(
        [merchant_cat[merchant_id], merchant_gmv, n_skus[merchant_id]]
    )
    X_user = np.column_stack([user_tenure[user_id], user_freq])

    logit0 = 0.95 * amount + 0.55 * merchant_cat[merchant_id] + 0.25 * user_tenure[user_id]
    if kind in ("concept_south", "both"):
        # Flip amount and drop the intercept so E[Y|south] actually moves.
        logit = logit0 * (1.0 - w) + (
            -0.95 * amount
            + 0.55 * merchant_cat[merchant_id]
            + 0.25 * user_tenure[user_id]
            - 1.15
        ) * w
    else:
        logit = logit0
    Y = rng.binomial(1, _sigmoid(logit)).astype(float)

    merchant_table = {
        "merchant_id": np.arange(n_merchants, dtype=int),
        "region": np.array(
            ["south" if i in set(south_ids.tolist()) else "north" for i in range(n_merchants)]
        ),
        "X": np.column_stack([merchant_cat, merchant_gmv0, n_skus]),
        "names": MERCHANT_FEATS,
    }
    user_table = {
        "user_id": np.arange(n_users, dtype=int),
        "X": np.column_stack([user_tenure, user_hist_freq]),
        "names": USER_FEATS,
    }
    meta = {
        "kind": kind,
        "onset_batch": int(onset_batch),
        "planted_subset": "south",
        "planted_order_feats": (
            ("amount", "channel") if kind in ("covariate_south", "both") else ("amount",)
        ),
        "planted_merchant_feats": (
            ("merchant_gmv",) if kind in ("covariate_south", "both") else ()
        ),
        "n_ref": n_ref,
        "n_new": n_new,
        "n_batches": n_batches,
        "n_merchants": n_merchants,
        "n_users": n_users,
        "title": f"order stream {kind} (south after batch {onset_batch})",
        "leakage": (
            "Y is outcome only; subset key is region, not Y; "
            "GMV is previous-batch amount EMA (no current batch/Y/future); "
            "hist_freq is a pre-period snapshot (no current/Y/future/clock)"
        ),
        "gmv_ema_lam": float(GMV_EMA_LAM),
        "feature_library": list(FEATURE_LIBRARY),
        "planted_how": planted_how_for_kind(kind),
    }
    tables = {
        "order_id": np.arange(n, dtype=int),
        "merchant_id": merchant_id,
        "user_id": user_id,
        "region": region,
        "X_order": X_order,
        "X_merchant": X_merchant,
        "X_user": X_user,
        "names_order": ORDER_FEATS,
        "names_merchant": MERCHANT_FEATS,
        "names_user": USER_FEATS,
        "Y": Y,
        "merchant_table": merchant_table,
        "user_table": user_table,
        "meta": meta,
    }
    return tables
