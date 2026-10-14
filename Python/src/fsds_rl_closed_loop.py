"""FSDS-empowered RL closed loop: three executable control modules.

Locked FSDS/TSS signs (do not invert):
  covariate P(X|W) large  → shrink η / prefer IW / refresh memory
  concept   P(Y|X) large  → raise η / relearn map / refresh reward components
  both quiet              → hold

Modules (existing controls, not a new loss):
  1. Meta-adapt   — few-shot inner SGD on a linear policy, conditioned on
                    FSDS task descriptor (top-K drift features + typed scores)
  2. Reward refresh — multi-component reward weights from typed drift
  3. Region LR boost — TSS map η_k = η0 (1+β δ_k)/(1+λ c_k) on feature groups

Runnable boards in this checkout (Meta-World / Procgen / D4RL not installed):
  bandit_cov     — contextual bandit, covariate shift on context mean
  bandit_concept — reward weights / arm map change (concept)
  bandit_mixed   — both channels
  maze_proxy     — D4RL-like offline support mismatch (state density shift)
  causal_proxy   — CausalWorld-like mass/friction parameter jump
  procgen_proxy  — Procgen-like layout/style OOD + sparse reward flip
  meta_vision    — Meta-Dataset-like cross-domain few-shot arm map

These proxies encode the *same drift taxonomy* the full simulators would
expose; swap the stream loader when gymnasium + those packages land.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import Ridge

SEED = 2026
ETA0 = 0.15
BETA = 2.5
LAMBDA = 1.5
ETA_MAX_MULT = 2.5
TAU_C = 0.15
TAU_D = 0.08
TOP_K = 8

BOARDS = (
    "bandit_cov",
    "bandit_concept",
    "bandit_mixed",
    "maze_proxy",
    "causal_proxy",
    "procgen_proxy",
    "meta_vision",
)
METHODS = (
    "static",           # no FSDS control
    "lr_boost_only",    # module 3
    "reward_only",      # module 2
    "meta_only",        # module 1
    "fsds_full",        # all three
    "naive_boost",      # WRONG sign: raise η on covariate (ablation)
)


# ---------------------------------------------------------------------------
# Streams: synthetic RL proxies with known typed drift
# ---------------------------------------------------------------------------


@dataclass
class RLStream:
    """Batch-indexed contextual stream for online policy / reward adaptation."""

    X: np.ndarray          # contexts / states
    A: np.ndarray          # actions taken (logged)
    R: np.ndarray          # observed scalar reward
    R_comp: np.ndarray     # (n, K) reward components before mixing
    batch: np.ndarray
    true_w: np.ndarray     # true reward mix weights per batch (n_batches, K)
    feature_groups: dict   # name -> slice into X columns
    meta: dict = field(default_factory=dict)


def _softmax(z, t=1.0):
    z = np.asarray(z, dtype=float) / max(float(t), 1e-8)
    z = z - z.max()
    e = np.exp(z)
    return e / (e.sum() + 1e-12)


def make_bandit_stream(
    kind="bandit_mixed",
    n_batches=8,
    n_per=200,
    p=24,
    n_arms=4,
    n_comp=3,
    seed=SEED,
):
    """Logged contextual bandit with typed drift between batches."""
    rng = np.random.default_rng(int(seed))
    groups = {
        "physics": slice(0, p // 3),
        "goal": slice(p // 3, 2 * p // 3),
        "agent": slice(2 * p // 3, p),
    }
    # base linear reward maps per arm / component
    Theta = rng.normal(scale=0.6, size=(n_arms, n_comp, p))
    w0 = np.array([0.5, 0.3, 0.2])
    Xs, As, Rs, Rcs, batches, true_ws = [], [], [], [], [], []
    for t in range(int(n_batches)):
        # --- covariate: move context mean on physics / goal ---
        mu = np.zeros(p)
        if kind in ("bandit_cov", "bandit_mixed"):
            mu[groups["physics"]] = 0.55 * np.sin(0.7 * t)
            mu[groups["goal"]] = 0.45 * (t / max(n_batches - 1, 1))
        X = rng.normal(size=(n_per, p)) + mu
        # --- concept: rotate reward component maps / mix weights ---
        w = w0.copy()
        Theta_t = Theta.copy()
        if kind in ("bandit_concept", "bandit_mixed"):
            w = _softmax(np.log(w0 + 1e-8) + 0.35 * np.array([t, -0.5 * t, 0.2 * t]))
            # flip goal-related columns on component 0
            rot = 1.0 + 0.25 * np.sin(0.9 * t)
            Theta_t[:, 0, groups["goal"]] *= rot
            if t >= n_batches // 2:
                Theta_t[:, 1, :] *= -0.3  # stale secondary reward head
        # score arm a by sum_k w_k * X·θ_{a,k}; behavior = noisier softmax
        arm_logits = np.zeros((n_per, n_arms))
        for a in range(n_arms):
            comps = X @ Theta_t[a].T  # (n, K)
            arm_logits[:, a] = comps @ w
        probs = np.exp(arm_logits - arm_logits.max(1, keepdims=True))
        probs = probs / probs.sum(1, keepdims=True)
        A = np.array([rng.choice(n_arms, p=probs[i]) for i in range(n_per)])
        R_comp = np.zeros((n_per, n_comp))
        for i in range(n_per):
            R_comp[i] = X[i] @ Theta_t[A[i]].T + 0.15 * rng.normal(size=n_comp)
        R = R_comp @ w + 0.05 * rng.normal(size=n_per)
        Xs.append(X)
        As.append(A)
        Rs.append(R)
        Rcs.append(R_comp)
        batches.append(np.full(n_per, t, dtype=int))
        true_ws.append(w)
    return RLStream(
        X=np.vstack(Xs),
        A=np.concatenate(As),
        R=np.concatenate(Rs),
        R_comp=np.vstack(Rcs),
        batch=np.concatenate(batches),
        true_w=np.stack(true_ws),
        feature_groups=groups,
        meta={"board": kind, "n_arms": n_arms, "n_comp": n_comp, "p": p},
    )


def make_maze_proxy(n_batches=8, n_per=200, p=20, seed=SEED):
    """D4RL-like support mismatch: early batches narrow maze, later wider.

    Action = discrete corridor choice; reward = progress − collision penalty.
    Covariate: state density expands; concept quiet unless flagged.
    """
    rng = np.random.default_rng(int(seed) + 3)
    groups = {"position": slice(0, 8), "velocity": slice(8, 14), "lidar": slice(14, p)}
    Xs, As, Rs, Rcs, batches, true_ws = [], [], [], [], [], []
    w = np.array([0.7, 0.2, 0.1])  # progress, energy, collision
    for t in range(int(n_batches)):
        scale = 0.4 + 0.9 * (t / max(n_batches - 1, 1))  # support expands
        X = rng.normal(scale=scale, size=(n_per, p))
        X[:, groups["position"]] += 0.3 * t
        A = rng.integers(0, 3, size=n_per)
        progress = X[:, 0] * (A == 0) + X[:, 1] * (A == 1) + 0.1 * rng.normal(size=n_per)
        energy = -0.2 * np.abs(X[:, groups["velocity"]]).mean(1)
        collision = -0.5 * (np.abs(X[:, groups["lidar"]]).max(1) > (1.2 + 0.1 * t))
        R_comp = np.c_[progress, energy, collision]
        R = R_comp @ w + 0.05 * rng.normal(size=n_per)
        Xs.append(X); As.append(A); Rs.append(R); Rcs.append(R_comp)
        batches.append(np.full(n_per, t, dtype=int)); true_ws.append(w.copy())
    return RLStream(
        X=np.vstack(Xs), A=np.concatenate(As), R=np.concatenate(Rs),
        R_comp=np.vstack(Rcs), batch=np.concatenate(batches),
        true_w=np.stack(true_ws), feature_groups=groups,
        meta={"board": "maze_proxy", "surrogate": "D4RL-support-mismatch"},
    )


def make_causal_proxy(n_batches=8, n_per=200, p=18, seed=SEED):
    """CausalWorld-like: mass/friction jumps (cov) + reward formula change (concept)."""
    rng = np.random.default_rng(int(seed) + 7)
    groups = {"mass": slice(0, 4), "friction": slice(4, 8), "joints": slice(8, p)}
    Xs, As, Rs, Rcs, batches, true_ws = [], [], [], [], [], []
    for t in range(int(n_batches)):
        mass = 1.0 + 0.4 * (t >= n_batches // 2)  # mid-stream mass jump
        fric = 0.3 + 0.25 * np.sin(0.8 * t)
        X = rng.normal(size=(n_per, p))
        X[:, groups["mass"]] *= mass
        X[:, groups["friction"]] *= fric
        A = rng.integers(0, 4, size=n_per)
        # components: reach, energy, safety
        reach = -(X[:, 8] - 0.2 * A) ** 2
        energy = -0.1 * (X[:, groups["joints"]] ** 2).mean(1)
        safety = -0.3 * (np.abs(X[:, groups["friction"]].mean(1)) > 0.5)
        R_comp = np.c_[reach, energy, safety]
        if t < n_batches // 2:
            w = np.array([0.6, 0.25, 0.15])
        else:
            # concept: safety becomes primary objective
            w = np.array([0.25, 0.15, 0.60])
        R = R_comp @ w + 0.05 * rng.normal(size=n_per)
        Xs.append(X); As.append(A); Rs.append(R); Rcs.append(R_comp)
        batches.append(np.full(n_per, t, dtype=int)); true_ws.append(w.copy())
    return RLStream(
        X=np.vstack(Xs), A=np.concatenate(As), R=np.concatenate(Rs),
        R_comp=np.vstack(Rcs), batch=np.concatenate(batches),
        true_w=np.stack(true_ws), feature_groups=groups,
        meta={"board": "causal_proxy", "surrogate": "CausalWorld-mass-friction"},
    )


def make_procgen_proxy(n_batches=8, n_per=200, p=22, seed=SEED):
    """Procgen-like: early easy layouts → later OOD layouts + reward sparsity flip.

    Groups mimic visual style / layout / dynamics; mid-stream concept changes
    which reward head dominates (dense progress → sparse gem).
    """
    rng = np.random.default_rng(int(seed) + 11)
    groups = {"style": slice(0, 7), "layout": slice(7, 14), "dynamics": slice(14, p)}
    Xs, As, Rs, Rcs, batches, true_ws = [], [], [], [], [], []
    for t in range(int(n_batches)):
        # covariate: layout support expands + style palette drifts
        style_mu = 0.5 * np.sin(1.1 * t)
        layout_scale = 0.5 + 1.1 * (t / max(n_batches - 1, 1))
        X = rng.normal(size=(n_per, p))
        X[:, groups["style"]] += style_mu
        X[:, groups["layout"]] *= layout_scale
        X[:, groups["dynamics"]] += 0.15 * t
        A = rng.integers(0, 4, size=n_per)
        progress = X[:, 7] * (A % 2 == 0) - 0.05 * np.abs(X[:, groups["dynamics"]]).mean(1)
        gem = (X[:, 10] > (0.4 + 0.05 * t)).astype(float) * (A == 2) - 0.1
        hazard = -0.4 * (np.abs(X[:, groups["style"]].mean(1)) > 0.8)
        R_comp = np.c_[progress, gem, hazard]
        if t < n_batches // 2:
            w = np.array([0.65, 0.20, 0.15])  # dense progress game
        else:
            w = np.array([0.20, 0.65, 0.15])  # sparse gem game (concept)
        R = R_comp @ w + 0.05 * rng.normal(size=n_per)
        Xs.append(X); As.append(A); Rs.append(R); Rcs.append(R_comp)
        batches.append(np.full(n_per, t, dtype=int)); true_ws.append(w.copy())
    return RLStream(
        X=np.vstack(Xs), A=np.concatenate(As), R=np.concatenate(Rs),
        R_comp=np.vstack(Rcs), batch=np.concatenate(batches),
        true_w=np.stack(true_ws), feature_groups=groups,
        meta={"board": "procgen_proxy", "surrogate": "Procgen-layout-OOD", "n_arms": 4},
    )


def make_meta_vision(n_batches=8, n_per=200, p=24, n_arms=5, seed=SEED):
    """Meta-Dataset-like: domain hops (cov) + class/arm map rotation (concept).

    Few-shot meta-adapt should win when concept rotates mid-stream; region LR
    should shrink on pure style covariate hops.
    """
    rng = np.random.default_rng(int(seed) + 13)
    groups = {"domain": slice(0, 8), "style": slice(8, 16), "classproto": slice(16, p)}
    n_comp = 3
    Theta = rng.normal(scale=0.55, size=(n_arms, n_comp, p))
    Xs, As, Rs, Rcs, batches, true_ws = [], [], [], [], [], []
    w0 = np.array([0.55, 0.30, 0.15])
    for t in range(int(n_batches)):
        # domain mean jumps every 2 batches (ImageNet→Aircraft→Birds…)
        domain_id = t // 2
        mu = np.zeros(p)
        mu[groups["domain"]] = 0.7 * domain_id
        mu[groups["style"]] = 0.35 * np.sin(1.3 * domain_id)
        X = rng.normal(size=(n_per, p)) + mu
        Theta_t = Theta.copy()
        w = w0.copy()
        if t >= n_batches // 2:
            # concept: permute arm preference via classproto columns
            Theta_t[:, :, groups["classproto"]] = -Theta_t[:, :, groups["classproto"]]
            w = _softmax(np.log(w0 + 1e-8) + np.array([0.4 * t, -0.3 * t, 0.1 * t]))
        arm_logits = np.zeros((n_per, n_arms))
        for a in range(n_arms):
            arm_logits[:, a] = (X @ Theta_t[a].T) @ w
        probs = np.exp(arm_logits - arm_logits.max(1, keepdims=True))
        probs = probs / probs.sum(1, keepdims=True)
        A = np.array([rng.choice(n_arms, p=probs[i]) for i in range(n_per)])
        R_comp = np.zeros((n_per, n_comp))
        for i in range(n_per):
            R_comp[i] = X[i] @ Theta_t[A[i]].T + 0.12 * rng.normal(size=n_comp)
        R = R_comp @ w + 0.04 * rng.normal(size=n_per)
        Xs.append(X); As.append(A); Rs.append(R); Rcs.append(R_comp)
        batches.append(np.full(n_per, t, dtype=int)); true_ws.append(w.copy())
    return RLStream(
        X=np.vstack(Xs), A=np.concatenate(As), R=np.concatenate(Rs),
        R_comp=np.vstack(Rcs), batch=np.concatenate(batches),
        true_w=np.stack(true_ws), feature_groups=groups,
        meta={"board": "meta_vision", "surrogate": "Meta-Dataset-domain-hop",
              "n_arms": n_arms, "n_comp": n_comp, "p": p},
    )


def load_board(board, **kw):
    board = str(board).lower()
    if board.startswith("bandit"):
        return make_bandit_stream(kind=board, **kw)
    if board == "maze_proxy":
        return make_maze_proxy(**kw)
    if board == "causal_proxy":
        return make_causal_proxy(**kw)
    if board == "procgen_proxy":
        return make_procgen_proxy(**kw)
    if board == "meta_vision":
        return make_meta_vision(**kw)
    raise ValueError("unknown board %s" % board)


# ---------------------------------------------------------------------------
# FSDS-style detectors on feature groups (reuse taxonomy, lightweight)
# ---------------------------------------------------------------------------


def group_covariate(X0, X1, groups):
    """|Cohen's d| on group means + RF domain VIMP mass as share."""
    c, share = {}, {}
    X0 = np.asarray(X0, dtype=float)
    X1 = np.asarray(X1, dtype=float)
    for name, sl in groups.items():
        a = X0[:, sl].mean(1)
        b = X1[:, sl].mean(1)
        pooled = float(np.sqrt(0.5 * (a.var() + b.var()) + 1e-8))
        c[name] = abs(float(a.mean() - b.mean()) / pooled)
    # RF-Domain VIMP over concatenated coords with group aggregation
    n0, n1 = len(X0), len(X1)
    X = np.vstack([X0, X1])
    y = np.array([0] * n0 + [1] * n1)
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=4, min_samples_leaf=3, random_state=SEED, n_jobs=1
    )
    clf.fit(X, y)
    vimp = clf.feature_importances_
    mass = {}
    for name, sl in groups.items():
        mass[name] = float(vimp[sl].sum())
    tot = sum(mass.values()) + 1e-12
    share = {k: mass[k] / tot for k in mass}
    # mix absolute intensity with share for descriptor stability
    for k in c:
        c[k] = float(0.5 * c[k] + 0.5 * share[k])
    return c, share


def group_concept(X0, y0, X1, y1, groups, A0=None, A1=None):
    """Excess risk after mean-align: concept intensity per feature group.

    Fit Ridge on old group features → residual on new (mean-aligned) features.
    """
    delta = {}
    X0 = np.asarray(X0, dtype=float)
    X1 = np.asarray(X1, dtype=float)
    y0 = np.asarray(y0, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    for name, sl in groups.items():
        Z0 = X0[:, sl].copy()
        Z1 = X1[:, sl].copy()
        # mean-align new to old support (covariate accounted for)
        Z1 = Z1 - Z1.mean(0) + Z0.mean(0)
        ridge = Ridge(alpha=1.0)
        ridge.fit(Z0, y0)
        e0 = float(np.mean((ridge.predict(Z0) - y0) ** 2))
        e1 = float(np.mean((ridge.predict(Z1) - y1) ** 2))
        delta[name] = float(max(e1 - e0, 0.0))
    # normalize to comparable scale
    mx = max(delta.values()) + 1e-12
    delta = {k: v / mx for k, v in delta.items()}
    return delta


def topk_drift_features(X0, X1, y0, y1, groups, k=TOP_K):
    """Coordinate-level drift scores for meta task descriptor."""
    c, _ = group_covariate(X0, X1, groups)
    d = group_concept(X0, y0, X1, y1, groups)
    # expand to coords: score = group typed intensity
    p = X0.shape[1]
    scores = np.zeros(p)
    types = np.zeros(p, dtype=int)  # 0 cov-dominant, 1 concept-dominant
    for name, sl in groups.items():
        s = float(c[name] + d[name])
        scores[sl] = s
        types[sl] = 1 if d[name] > c[name] else 0
    idx = np.argsort(-scores)[: int(k)]
    return {
        "indices": idx.tolist(),
        "scores": scores[idx].tolist(),
        "types": types[idx].tolist(),  # 0 cov, 1 concept
        "c": c,
        "delta": d,
    }


def task_descriptor(drift, dim=16):
    """Encode top-K drift into a fixed vector for meta conditioning."""
    v = np.zeros(dim, dtype=float)
    scores = np.asarray(drift["scores"], dtype=float)
    types = np.asarray(drift["types"], dtype=float)
    n = min(len(scores), dim // 2)
    v[:n] = scores[:n]
    v[dim // 2 : dim // 2 + n] = types[:n]
    # append aggregate typed intensities
    cvals = list(drift["c"].values())
    dvals = list(drift["delta"].values())
    if dim >= 4:
        v[-4] = float(np.mean(cvals))
        v[-3] = float(np.max(cvals))
        v[-2] = float(np.mean(dvals))
        v[-1] = float(np.max(dvals))
    return v


# ---------------------------------------------------------------------------
# Module 3: typed region LR (LOCKED signs)
# ---------------------------------------------------------------------------


def tss_region_lr(c, delta, eta0=ETA0, beta=BETA, lam=LAMBDA):
    """η_g = clip( η0 (1+βδ)/(1+λc) ). Cov↓ η, concept↑ η."""
    out = {}
    for g in c:
        cg = float(c.get(g, 0.0))
        dg = float(delta.get(g, 0.0))
        if cg < TAU_C and dg < TAU_D:
            out[g] = float(eta0)  # hold-ish baseline
        else:
            eta = float(eta0) * (1.0 + float(beta) * dg) / (1.0 + float(lam) * cg)
            out[g] = float(np.clip(eta, 0.05 * eta0, ETA_MAX_MULT * eta0))
    return out


def naive_boost_lr(c, delta, eta0=ETA0):
    """Ablation: raise η on any drift (WRONG for covariate)."""
    out = {}
    for g in c:
        s = float(c.get(g, 0.0) + delta.get(g, 0.0))
        out[g] = float(np.clip(eta0 * (1.0 + 1.5 * s), 0.05 * eta0, ETA_MAX_MULT * eta0))
    return out


# ---------------------------------------------------------------------------
# Module 2: reward component weight refresh
# ---------------------------------------------------------------------------


def _component_typed_scores(Rc0, Rc1, R0, R1, seed=SEED):
    """Per-component (c_k, δ_k): value-shift vs mix-role change vs observed R."""
    Rc0 = np.asarray(Rc0, dtype=float)
    Rc1 = np.asarray(Rc1, dtype=float)
    R0 = np.asarray(R0, dtype=float)
    R1 = np.asarray(R1, dtype=float)
    K = Rc0.shape[1]
    c_scores, d_scores = [], []
    for k in range(K):
        a, b = Rc0[:, k], Rc1[:, k]
        pooled = float(np.sqrt(0.5 * (a.var() + b.var()) + 1e-8))
        c_k = abs(float(a.mean() - b.mean()) / pooled)
        # concept: correlation / ridge role of component k in explaining R jumps
        # fit R ~ Rc on old; residual on new (after mean-aligning Rc)
        ridge = Ridge(alpha=0.5)
        ridge.fit(Rc0, R0)
        e0 = float(np.mean((ridge.predict(Rc0) - R0) ** 2))
        Rc1_aligned = Rc1 - Rc1.mean(0) + Rc0.mean(0)
        e1 = float(np.mean((ridge.predict(Rc1_aligned) - R1) ** 2))
        # leave-one-component ablation excess on new batch
        mask = np.ones(K, dtype=bool)
        mask[k] = False
        if mask.any():
            r_ab = Ridge(alpha=0.5).fit(Rc0[:, mask], R0)
            e_ab = float(np.mean((r_ab.predict(Rc1_aligned[:, mask]) - R1) ** 2))
            importance_jump = max(e_ab - e1, 0.0)  # k still useful on new
        else:
            importance_jump = 0.0
        d_k = max(e1 - e0, 0.0) / (e0 + 1e-6) + 0.25 * importance_jump
        c_scores.append(float(c_k))
        d_scores.append(float(d_k))
    return np.asarray(c_scores), np.asarray(d_scores)


def refresh_reward_weights(Rc0, Rc1, R0, R1, old_w, alpha_concept=0.75, beta_cov=0.25, seed=SEED):
    """Typed refresh of mix weights over reward components.

    Locked polarity:
      concept δ large → move mix toward empirical Ridge(R ~ Rc) on the *new*
                        batch (relearn the mix / stale head replaced).
      covariate c large, δ quiet → keep old mix, mild upweight of shifted
                        heads (reweight, do not relearn).
    """
    Rc0 = np.asarray(Rc0, dtype=float)
    Rc1 = np.asarray(Rc1, dtype=float)
    R0 = np.asarray(R0, dtype=float)
    R1 = np.asarray(R1, dtype=float)
    K = Rc0.shape[1]
    old_w = np.asarray(old_w, dtype=float)
    old_w = old_w / (old_w.sum() + 1e-12)

    # empirical non-negative mix on new batch via NNLS-style clipped ridge coefs
    # fit on first half, validate on second half — reject if worse than old_w
    n = len(R1)
    split = max(n // 2, 8)
    emp_raw = Ridge(alpha=1.0, fit_intercept=False).fit(Rc1[:split], R1[:split]).coef_.astype(float)
    emp_raw = np.clip(emp_raw, 1e-4, None)
    emp = emp_raw / emp_raw.sum()
    mse_emp = float(np.mean((Rc1[split:] @ emp - R1[split:]) ** 2))
    mse_old = float(np.mean((Rc1[split:] @ old_w - R1[split:]) ** 2))
    accept_emp = mse_emp <= mse_old * 1.05  # allow tiny slack

    c_k, d_k = _component_typed_scores(Rc0, Rc1, R0, R1, seed=seed)
    d_bar = float(np.mean(d_k))
    c_bar = float(np.mean(c_k))
    # gate: how hard to trust empirical re-fit (0 if holdout rejects)
    gate = float(np.clip(alpha_concept * np.tanh(d_bar) / (1.0 + 0.5 * c_bar), 0.0, 0.95))
    if not accept_emp:
        gate = 0.0
    new_w = (1.0 - gate) * old_w + gate * emp
    # covariate-only polish: mild raise of high-c components when concept quiet
    if d_bar < TAU_D and c_bar >= TAU_C:
        bump = 1.0 + beta_cov * np.tanh(c_k)
        new_w = new_w * bump
    new_w = np.clip(new_w, 1e-3, None)
    return new_w / new_w.sum()


# ---------------------------------------------------------------------------
# Module 1: meta-adapt — few-shot linear policy conditioned on descriptor
# ---------------------------------------------------------------------------


class LinearPolicy:
    """Softmax policy over arms: logits = X @ W[:, a] + b[a]."""

    def __init__(self, p, n_arms, seed=SEED):
        rng = np.random.default_rng(int(seed))
        self.W = rng.normal(scale=0.05, size=(p, n_arms))
        self.b = np.zeros(n_arms)

    def logits(self, X):
        return np.asarray(X, dtype=float) @ self.W + self.b

    def probs(self, X):
        z = self.logits(X)
        z = z - z.max(1, keepdims=True)
        e = np.exp(z)
        return e / (e.sum(1, keepdims=True) + 1e-12)

    def act(self, X, rng):
        P = self.probs(X)
        return np.array([rng.choice(P.shape[1], p=P[i]) for i in range(len(X))])

    def copy(self):
        q = LinearPolicy.__new__(LinearPolicy)
        q.W = self.W.copy()
        q.b = self.b.copy()
        return q


def policy_grad_step(policy, X, A, adv, lr_vec, groups):
    """REINFORCE-style step with per-group LR from TSS."""
    X = np.asarray(X, dtype=float)
    A = np.asarray(A, dtype=int)
    adv = np.asarray(adv, dtype=float)
    P = policy.probs(X)
    n, p = X.shape
    n_arms = P.shape[1]
    # d log π / d W[:,a] = 1[A=a] X - P[:,a] X
    dW = np.zeros_like(policy.W)
    db = np.zeros_like(policy.b)
    for i in range(n):
        for a in range(n_arms):
            ind = 1.0 if A[i] == a else 0.0
            coef = (ind - P[i, a]) * adv[i]
            dW[:, a] += coef * X[i]
            db[a] += coef
    dW /= max(n, 1)
    db /= max(n, 1)
    # apply group LRs
    for g, sl in groups.items():
        lr = float(lr_vec.get(g, ETA0))
        policy.W[sl, :] += lr * dW[sl, :]
    policy.b += float(np.mean(list(lr_vec.values()))) * db
    return policy


def meta_adapt(policy, X, A, R, descriptor, lr_vec, groups, n_inner=5):
    """Few-shot inner updates; descriptor scales η further (meta conditioning)."""
    pol = policy.copy()
    # descriptor aggregates: last two entries are mean/max concept
    scale = 1.0 + 0.5 * float(descriptor[-2]) + 0.25 * float(descriptor[-1])
    lr_scaled = {g: float(np.clip(lr_vec[g] * scale, 1e-4, ETA_MAX_MULT * ETA0)) for g in lr_vec}
    adv = R - float(R.mean())
    for _ in range(int(n_inner)):
        pol = policy_grad_step(pol, X, A, adv, lr_scaled, groups)
    return pol


# ---------------------------------------------------------------------------
# Online closed loop
# ---------------------------------------------------------------------------


def evaluate_policy(policy, X, A_star_scores):
    """Expected reward under policy if arm scores are known (oracle eval)."""
    P = policy.probs(X)
    return float((P * A_star_scores).sum(1).mean())


def offline_value(policy, X, A, R):
    """Importance-weighted offline value of logged rewards (self-norm IW)."""
    P = policy.probs(X)
    # behavior approx uniform-ish fallback
    n_arms = P.shape[1]
    pi_b = np.full(len(A), 1.0 / n_arms)
    pi_a = P[np.arange(len(A)), A]
    w = pi_a / (pi_b + 1e-8)
    w = w / (w.mean() + 1e-12)
    return float(np.mean(w * R))


def run_closed_loop(stream, method="fsds_full", seed=SEED):
    """Hop over batches: detect → refresh reward → LR → meta-adapt → score."""
    X = stream.X
    A = stream.A
    R = stream.R
    Rc = stream.R_comp
    batch = stream.batch
    groups = stream.feature_groups
    k = int(batch.max()) + 1
    n_arms = int(stream.meta.get("n_arms", int(A.max()) + 1))
    p = X.shape[1]
    Kc = Rc.shape[1]
    policy = LinearPolicy(p, n_arms, seed=seed)
    w_reward = np.ones(Kc) / Kc
    rng = np.random.default_rng(int(seed))
    history = []
    # warm start on batch 0
    i0 = batch == 0
    adv0 = R[i0] - R[i0].mean()
    lr0 = {g: ETA0 for g in groups}
    for _ in range(8):
        policy = policy_grad_step(policy, X[i0], A[i0], adv0, lr0, groups)

    for t in range(1, k):
        ip, ic = batch == (t - 1), batch == t
        X0, X1 = X[ip], X[ic]
        R0, R1 = R[ip], R[ic]
        A0, A1 = A[ip], A[ic]
        Rc0, Rc1 = Rc[ip], Rc[ic]

        drift = topk_drift_features(X0, X1, R0, R1, groups)
        c, delta = drift["c"], drift["delta"]
        desc = task_descriptor(drift)

        # --- module 2 ---
        if method in ("reward_only", "fsds_full"):
            w_reward = refresh_reward_weights(
                Rc0, Rc1, R0, R1, w_reward, seed=seed + t
            )
        R1_hat = Rc1 @ w_reward
        # reward-mix fidelity (lower better)
        mix_mse = float(np.mean((R1_hat - R1) ** 2))

        # --- module 3 ---
        if method == "naive_boost":
            lr_vec = naive_boost_lr(c, delta)
        elif method in ("lr_boost_only", "fsds_full", "meta_only"):
            lr_vec = tss_region_lr(c, delta)
        else:
            lr_vec = {g: ETA0 for g in groups}

        # --- module 1 ---
        # reward_only / fsds_full train on refreshed mix; others on observed R
        train_R = R1_hat if method in ("reward_only", "fsds_full") else R1
        if method in ("meta_only", "fsds_full"):
            policy = meta_adapt(policy, X1, A1, train_R, desc, lr_vec, groups, n_inner=5)
        elif method in ("lr_boost_only", "reward_only", "naive_boost"):
            adv = train_R - train_R.mean()
            for _ in range(5):
                policy = policy_grad_step(policy, X1, A1, adv, lr_vec, groups)
        else:  # static: tiny update with base lr only
            adv = R1 - R1.mean()
            lr_base = {g: ETA0 for g in groups}
            policy = policy_grad_step(policy, X1, A1, adv, lr_base, groups)

        v_off = offline_value(policy, X1, A1, R1)
        # oracle arm scores: batch mean reward per arm as cheap online proxy
        arm_scores = np.zeros((len(A1), n_arms))
        for a in range(n_arms):
            mask = A1 == a
            mu_a = float(R1[mask].mean()) if mask.any() else float(R1.mean())
            arm_scores[:, a] = mu_a
        v_on = evaluate_policy(policy, X1, arm_scores)

        w_err = float(np.linalg.norm(w_reward - stream.true_w[t]))
        history.append(
            {
                "round": int(t),
                "offline_value": v_off,
                "online_proxy": v_on,
                "reward_w_err": w_err,
                "mix_mse": mix_mse,
                "c": {kk: float(vv) for kk, vv in c.items()},
                "delta": {kk: float(vv) for kk, vv in delta.items()},
                "lr": {kk: float(vv) for kk, vv in lr_vec.items()},
                "w_reward": w_reward.tolist(),
                "mean_c": float(np.mean(list(c.values()))),
                "mean_delta": float(np.mean(list(delta.values()))),
            }
        )
    offline = np.array([h["offline_value"] for h in history], dtype=float)
    mix_mses = np.array([h["mix_mse"] for h in history], dtype=float)
    return {
        "method": method,
        "board": stream.meta.get("board"),
        "online_value": float(offline.mean()) if offline.size else float("nan"),
        "cum_value": float(offline.sum()) if offline.size else float("nan"),
        "last_value": float(offline[-1]) if offline.size else float("nan"),
        "mean_w_err": float(np.mean([h["reward_w_err"] for h in history])) if history else float("nan"),
        "mean_mix_mse": float(np.mean(mix_mses)) if mix_mses.size else float("nan"),
        "path": offline.tolist(),
        "history": history,
        "meta": dict(stream.meta),
    }


def run_suite(boards=None, methods=None, seeds=None, n_per=200, n_batches=8):
    boards = list(boards or BOARDS)
    methods = list(methods or METHODS)
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    rows = []
    traces = {b: {m: [] for m in methods} for b in boards}
    for board in boards:
        for s in seeds:
            stream = load_board(board, n_per=n_per, n_batches=n_batches, seed=int(s))
            for method in methods:
                rec = run_closed_loop(stream, method=method, seed=int(s))
                traces[board][method].append(rec)
                rows.append(
                    {
                        "board": board,
                        "seed": int(s),
                        "method": method,
                        "online_value": rec["online_value"],
                        "cum_value": rec["cum_value"],
                        "mean_w_err": rec["mean_w_err"],
                        "mean_mix_mse": rec["mean_mix_mse"],
                    }
                )
    table = {}
    for board in boards:
        table[board] = {}
        for method in methods:
            v = np.array(
                [r["online_value"] for r in rows if r["board"] == board and r["method"] == method],
                dtype=float,
            )
            we = np.array(
                [r["mean_w_err"] for r in rows if r["board"] == board and r["method"] == method],
                dtype=float,
            )
            mm = np.array(
                [r["mean_mix_mse"] for r in rows if r["board"] == board and r["method"] == method],
                dtype=float,
            )
            table[board][method] = {
                "online_value": {
                    "mean": float(np.nanmean(v)),
                    "sd": float(np.nanstd(v, ddof=1) if len(v) > 1 else 0.0),
                },
                "mean_w_err": {
                    "mean": float(np.nanmean(we)),
                    "sd": float(np.nanstd(we, ddof=1) if len(we) > 1 else 0.0),
                },
                "mean_mix_mse": {
                    "mean": float(np.nanmean(mm)),
                    "sd": float(np.nanstd(mm, ddof=1) if len(mm) > 1 else 0.0),
                },
                "n": int(len(v)),
            }
    return {
        "boards": boards,
        "methods": methods,
        "seeds": seeds,
        "table": table,
        "rows": rows,
        "traces": traces,
    }


def plot_suite(suite, path):
    import matplotlib.pyplot as plt

    boards = suite["boards"]
    methods = suite["methods"]
    colors = {
        "static": "#9AA3AE",
        "lr_boost_only": "#C4892A",
        "reward_only": "#2F6B4F",
        "meta_only": "#5B3A8C",
        "fsds_full": "#C45C26",
        "naive_boost": "#A33B24",
    }
    fig, axes = plt.subplots(1, len(boards), figsize=(3.2 * len(boards), 4.2), squeeze=False)
    for j, board in enumerate(boards):
        ax = axes[0, j]
        means = [suite["table"][board][m]["online_value"]["mean"] for m in methods]
        sds = [suite["table"][board][m]["online_value"]["sd"] for m in methods]
        ypos = np.arange(len(methods))[::-1]
        ax.barh(
            ypos, means, xerr=sds,
            color=[colors.get(m, "#333") for m in methods], height=0.62,
            error_kw={"ecolor": "#555", "lw": 0.9, "capsize": 2},
        )
        ax.set_yticks(ypos)
        ax.set_yticklabels(methods, fontsize=7.5)
        ax.set_xlabel("offline IW value ↑")
        ax.set_title(board, loc="left", fontsize=10, fontweight="bold")
        ax.grid(True, axis="x", alpha=0.35)
    fig.suptitle(
        "FSDS→RL closed loop (typed TSS signs) — offline value",
        fontsize=12, fontweight="bold", y=0.98,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def write_readme(suite, path):
    lines = [
        "# FSDS-empowered RL closed loop",
        "",
        "Three modules driven by FSDS typed drift `(c, δ)` — **locked signs**:",
        "covariate large → shrink η; concept large → raise η / refresh reward heads.",
        "",
        "## Modules",
        "",
        "| id | role |",
        "| --- | --- |",
        "| meta_adapt | few-shot policy inner-loop conditioned on drift descriptor |",
        "| reward_refresh | mix weights over reward components from typed drift |",
        "| region_lr | TSS η per feature group |",
        "",
        "## Boards (proxies — full sims not installed here)",
        "",
        "| board | stands in for | drift |",
        "| --- | --- | --- |",
        "| bandit_cov | Meta-World task params | covariate |",
        "| bandit_concept | reward / map change | concept |",
        "| bandit_mixed | both | both |",
        "| maze_proxy | D4RL support mismatch | covariate |",
        "| causal_proxy | CausalWorld mass/friction + objective flip | both |",
        "| procgen_proxy | Procgen layout/style OOD + gem reward flip | both |",
        "| meta_vision | Meta-Dataset domain hop + class map rotation | both |",
        "",
        "## Results (offline IW value ↑ better)",
        "",
    ]
    for board in suite["boards"]:
        lines.append("### `%s`" % board)
        lines.append("")
        lines.append("| method | value ↑ | w_err ↓ | mix_mse ↓ |")
        lines.append("| --- | --- | --- | --- |")
        for m in suite["methods"]:
            t = suite["table"][board][m]
            lines.append(
                "| %s | %.4f (%.4f) | %.4f | %.4f |"
                % (
                    m,
                    t["online_value"]["mean"],
                    t["online_value"]["sd"],
                    t["mean_w_err"]["mean"],
                    t["mean_mix_mse"]["mean"],
                )
            )
        lines.append("")
    lines.extend(
        [
            "## Run",
            "",
            "```bash",
            "python3 scripts/run_fsds_rl_closed_loop.py",
            "python3 scripts/run_fsds_rl_closed_loop.py --quick",
            "```",
            "",
            "Code: `Python/src/fsds_rl_closed_loop.py`",
            "",
            "## Opportunity notes",
            "",
            "- Install `gymnasium` + Meta-World / Procgen / D4RL / CausalWorld to replace proxies.",
            "- Existing FSDS boards (MSR-VTT modality π, Amazon hop IW) already implement",
            "  module-3 style control; this package adds reward refresh + meta inner-loop.",
            "- Do **not** use `naive_boost` in production — it raises η on covariate shift.",
            "",
        ]
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _c(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))

    slim = {
        "boards": obj["boards"],
        "methods": obj["methods"],
        "seeds": obj["seeds"],
        "table": obj["table"],
        "rows": obj["rows"],
    }
    path.write_text(json.dumps(slim, indent=2, default=_c))
    return path
