"""Instrumental-variable selection + 2SLS on synthetic confounded uplift.

This is *instrumental variable* (Z ⟂ ε | …, Z → T), not “important variable”
VIMP. Under unobserved confounding U → (T, Y), naive meta-learners / OLS for
the treatment effect are biased; a valid instrument recovers the LATE / ATE
linear target.

Synthetic argument (this module):
  1. DGP with confounder U, valid Z, invalid Z (exclusion broken), weak Z
  2. Naive OLS / T-side residualization vs 2SLS
  3. IV selection by first-stage partial-F (+ exclusion residual screen)
  4. Suite: oracle-valid / selected / all-Z / invalid-only / naive
  5. E-value + Rosenbaum Γ sensitivity on the naive association
  6. Subset localization: IV strength heterogeneity across segments

Locked takeaway for the FSDS product case:
  when treatment is confounded, uplift monitors need an IV channel;
  selection must prefer strong + exclusion-plausible Z, not raw correlation;
  report how fragile the naive effect is (E-value / Γ) and where Z is strong.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression

SEED = 2026
F_WEAK = 10.0  # Staiger–Stock rule-of-thumb for first-stage strength
ROSENBAUM_GAMMAS = (
    1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0,
)


# ---------------------------------------------------------------------------
# Synthetic DGP
# ---------------------------------------------------------------------------


@dataclass
class IVBundle:
    X: np.ndarray          # observed covariates (no U)
    Z: np.ndarray          # candidate instruments (valid | invalid | weak)
    T: np.ndarray          # endogenous treatment
    Y: np.ndarray          # outcome
    U: np.ndarray          # hidden confounder (eval only)
    true_ate: float
    z_meta: dict           # column roles
    meta: dict


def make_iv_dgp(
    n=2000,
    p_x=6,
    n_valid=3,
    n_invalid=3,
    n_weak=2,
    tau=1.5,
    conf_strength=1.0,
    valid_strength=1.0,
    invalid_y=0.8,
    weak_strength=0.05,
    het_strength=1.2,
    seed=SEED,
):
    """Linear IV world with one hidden confounder + segment IV-strength het.

    T = X β_t + s(X0)·Z_valid γ_v + Z_weak γ_w + Z_invalid γ_i + α U + ε_t
    Y = τ T + X β_y + Z_invalid δ + α U + ε_y

    ``s(X0) = 1 + het_strength·σ(X0)`` scales the valid first stage so high-X0
    segments have stronger instruments (subset localization target).

    Valid Z: affect T, excluded from Y.
    Invalid Z: appear in Y directly (exclusion violated) and may weakly enter T.
    Weak Z: tiny first stage — fail F-screen.
    """
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(n, p_x))
    U = rng.normal(size=n)
    Zv = rng.normal(size=(n, n_valid))
    Zi = rng.normal(size=(n, n_invalid))
    Zw = rng.normal(size=(n, n_weak))

    beta_t = rng.normal(scale=0.3, size=p_x)
    beta_y = rng.normal(scale=0.3, size=p_x)
    gamma_v = valid_strength * (0.6 + 0.4 * rng.random(n_valid))
    gamma_w = weak_strength * (0.5 + rng.random(n_weak))
    # invalid instruments: small leak into T + direct path into Y
    gamma_i = 0.15 * rng.normal(size=n_invalid)
    delta_i = invalid_y * (0.5 + 0.5 * rng.random(n_invalid))

    # segment scaler on valid first stage (logistic in X0)
    seg_scale = 1.0 + float(het_strength) / (1.0 + np.exp(-X[:, 0]))
    eps_t = rng.normal(scale=0.5, size=n)
    eps_y = rng.normal(scale=0.5, size=n)
    T = (
        X @ beta_t
        + (Zv * seg_scale[:, None]) @ gamma_v
        + Zw @ gamma_w
        + Zi @ gamma_i
        + float(conf_strength) * U
        + eps_t
    )
    Y = (
        float(tau) * T
        + X @ beta_y
        + Zi @ delta_i
        + float(conf_strength) * U
        + eps_y
    )
    Z = np.hstack([Zv, Zi, Zw])
    roles = (
        ["valid"] * n_valid + ["invalid"] * n_invalid + ["weak"] * n_weak
    )
    # quartile segment labels on X0 for subset localization
    edges = np.quantile(X[:, 0], [0.0, 0.25, 0.5, 0.75, 1.0])
    seg = np.digitize(X[:, 0], edges[1:-1], right=True)  # 0..3
    return IVBundle(
        X=X,
        Z=Z,
        T=T,
        Y=Y,
        U=U,
        true_ate=float(tau),
        z_meta={
            "roles": roles,
            "valid_idx": list(range(0, n_valid)),
            "invalid_idx": list(range(n_valid, n_valid + n_invalid)),
            "weak_idx": list(
                range(n_valid + n_invalid, n_valid + n_invalid + n_weak)
            ),
            "names": [f"Z_{r}{i}" for i, r in enumerate(roles)],
            "segment": seg.astype(int),
            "segment_edges": edges.tolist(),
            "segment_names": ["Q1_low", "Q2", "Q3", "Q4_high"],
            "seg_scale": seg_scale,
        },
        meta={
            "n": n,
            "p_x": p_x,
            "conf_strength": float(conf_strength),
            "valid_strength": float(valid_strength),
            "invalid_y": float(invalid_y),
            "weak_strength": float(weak_strength),
            "het_strength": float(het_strength),
            "seed": int(seed),
        },
    )


# ---------------------------------------------------------------------------
# First stage / IV diagnostics
# ---------------------------------------------------------------------------


def _resid(y, X):
    y = np.asarray(y, dtype=float).reshape(-1)
    if X is None or np.asarray(X).size == 0:
        return y - y.mean()
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    reg = LinearRegression().fit(X, y)
    return y - reg.predict(X)


def partial_f_stat(T, Z_j, controls):
    """Partial F for one instrument after residualizing on controls."""
    T = np.asarray(T, dtype=float).reshape(-1)
    Z_j = np.asarray(Z_j, dtype=float).reshape(-1)
    n = len(T)
    # residualize on controls
    Tr = _resid(T, controls)
    Zr = _resid(Z_j, controls)
    # restricted: Tr ~ 1  (already demeaned by resid)
    ssr_r = float(np.sum(Tr ** 2))
    # unrestricted: Tr ~ Zr
    if float(np.var(Zr)) < 1e-12:
        return 0.0, 0.0
    b = float(np.dot(Zr, Tr) / (np.dot(Zr, Zr) + 1e-12))
    fitted = b * Zr
    ssr_u = float(np.sum((Tr - fitted) ** 2))
    df_u = max(n - 2 - (0 if controls is None else np.asarray(controls).shape[1]), 1)
    f = ((ssr_r - ssr_u) / 1.0) / (ssr_u / df_u + 1e-12)
    r2 = 1.0 - ssr_u / (ssr_r + 1e-12)
    return float(max(f, 0.0)), float(max(r2, 0.0))


def _slope(y, x):
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    v = float(np.dot(x, x))
    if v < 1e-12:
        return 0.0
    return float(np.dot(x, y) / v)


def exclusion_excess_score(Y, T, Z_j, X, tau_pilot):
    """Exclusion screen via reduced-form excess.

    After residualizing on X:
      FS = slope(T̃, Z̃),  RF = slope(Ỹ, Z̃)
    Valid IV ⇒ RF ≈ τ · FS.  Invalid IV adds a direct δ ⇒ |RF − τ·FS| large.
    ``tau_pilot`` is a rough scale (e.g. OLS); only ranks candidates.
    """
    Tr = _resid(T, X)
    Yr = _resid(Y, X)
    Zr = _resid(Z_j, X)
    fs = _slope(Tr, Zr)
    rf = _slope(Yr, Zr)
    return float(abs(rf - float(tau_pilot) * fs))


def score_instruments(bundle, f_min=F_WEAK, excl_max=None):
    """Score every Z column: first-stage F, partial R², exclusion excess."""
    X, Z, T, Y = bundle.X, bundle.Z, bundle.T, bundle.Y
    k = Z.shape[1]
    tau_pilot = ols_ate(Y, T, X)
    raw = []
    for j in range(k):
        f, r2 = partial_f_stat(T, Z[:, j], X)
        excl = exclusion_excess_score(Y, T, Z[:, j], X, tau_pilot=tau_pilot)
        raw.append((j, f, r2, excl))
    # adaptive exclusion cutoff: keep the low-excess cluster among strong Z
    strong_ex = [excl for _, f, _, excl in raw if f >= float(f_min)]
    if excl_max is None:
        if strong_ex:
            excl_max = float(np.median(strong_ex) + 0.75 * np.std(strong_ex))
            excl_max = max(excl_max, 0.05)
        else:
            excl_max = 0.15
    rows = []
    for j, f, r2, excl in raw:
        strong = f >= float(f_min)
        clean = excl <= float(excl_max)
        rows.append(
            {
                "j": int(j),
                "name": bundle.z_meta["names"][j],
                "role": bundle.z_meta["roles"][j],
                "partial_F": float(f),
                "partial_R2": float(r2),
                "exclusion_score": float(excl),
                "excl_max": float(excl_max),
                "pass_strength": bool(strong),
                "pass_exclusion": bool(clean),
                "selected": bool(strong and clean),
            }
        )
    return rows


def select_instruments(scores, mode="auto"):
    """Return selected column indices.

    modes:
      auto     — pass F-screen AND exclusion screen
      strong   — F-screen only
      oracle   — roles == valid (eval only)
      invalid  — roles == invalid (ablation)
      all      — every column
      none     — empty
    """
    mode = str(mode).lower()
    if mode == "auto":
        return [r["j"] for r in scores if r["selected"]]
    if mode == "strong":
        return [r["j"] for r in scores if r["pass_strength"]]
    if mode == "oracle":
        return [r["j"] for r in scores if r["role"] == "valid"]
    if mode == "invalid":
        return [r["j"] for r in scores if r["role"] == "invalid"]
    if mode == "weak":
        return [r["j"] for r in scores if r["role"] == "weak"]
    if mode == "all":
        return [r["j"] for r in scores]
    if mode == "none":
        return []
    raise ValueError("unknown IV select mode %s" % mode)


# ---------------------------------------------------------------------------
# Estimators
# ---------------------------------------------------------------------------


def ols_ate(Y, T, X):
    """Naive controlled regression Y ~ T + X (biased under confounding)."""
    X = np.asarray(X, dtype=float)
    T = np.asarray(T, dtype=float).reshape(-1, 1)
    Y = np.asarray(Y, dtype=float).reshape(-1)
    design = np.hstack([T, X])
    coef = LinearRegression().fit(design, Y).coef_
    return float(coef[0])


def tlearner_ate(Y, T, X, treash=None):
    """Simple continuous-T proxy: residual-on-residual OLS (still confounded)."""
    Tr = _resid(T, X)
    Yr = _resid(Y, X)
    if float(np.var(Tr)) < 1e-12:
        return float("nan")
    return float(np.dot(Tr, Yr) / (np.dot(Tr, Tr) + 1e-12))


def iv_2sls_ate(Y, T, X, Z):
    """Classic 2SLS: T̂ = P_Z T (with X), then Y ~ T̂ + X."""
    X = np.asarray(X, dtype=float)
    Z = np.asarray(Z, dtype=float)
    T = np.asarray(T, dtype=float).reshape(-1)
    Y = np.asarray(Y, dtype=float).reshape(-1)
    if Z.ndim == 1:
        Z = Z[:, None]
    if Z.shape[1] == 0:
        return float("nan"), {"stage1_f": 0.0, "n_z": 0}
    # stage 1: T ~ X + Z
    stage1_X = np.hstack([X, Z])
    m1 = LinearRegression().fit(stage1_X, T)
    T_hat = m1.predict(stage1_X)
    # joint first-stage F for Z block
    ssr_r = float(np.sum(_resid(T, X) ** 2))
    ssr_u = float(np.sum((T - T_hat) ** 2))
    q = Z.shape[1]
    n = len(T)
    df_u = max(n - X.shape[1] - q - 1, 1)
    stage1_f = ((ssr_r - ssr_u) / max(q, 1)) / (ssr_u / df_u + 1e-12)
    # stage 2: Y ~ T_hat + X
    m2 = LinearRegression().fit(np.hstack([T_hat.reshape(-1, 1), X]), Y)
    return float(m2.coef_[0]), {
        "stage1_f": float(max(stage1_f, 0.0)),
        "n_z": int(q),
        "stage1_r2": float(max(1.0 - ssr_u / (ssr_r + 1e-12), 0.0)),
    }


# ---------------------------------------------------------------------------
# E-value (VanderWeele & Ding) + Rosenbaum Γ sensitivity
# ---------------------------------------------------------------------------


def approx_rr_from_continuous(beta, sd_y):
    """Approximate risk-ratio scale for a continuous effect (Ding & VanderWeele).

    RR ≈ exp(0.91 · |β| / σ_Y). Used only to put E-values on a familiar scale.
    """
    sd = max(float(sd_y), 1e-8)
    return float(np.exp(0.91 * abs(float(beta)) / sd))


def evalue_from_rr(rr):
    """Point E-value for an approximate risk ratio (≥1)."""
    rr = float(rr)
    if rr < 1.0:
        rr = 1.0 / max(rr, 1e-12)
    if rr <= 1.0:
        return 1.0
    return float(rr + np.sqrt(rr * (rr - 1.0)))


def evalue_analysis(Y, T, X, ate_hat=None):
    """E-value for the naive controlled association (how strong U must be).

    Reports:
      - point E-value for the OLS/partial effect
      - E-value for the CI lower bound (null-crossing fragility)
    """
    Y = np.asarray(Y, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    X = np.asarray(X, dtype=float)
    if ate_hat is None:
        ate_hat = ols_ate(Y, T, X)
    # residual SE of Y | T,X for RR approximation + Wald CI
    design = np.hstack([T.reshape(-1, 1), X])
    reg = LinearRegression().fit(design, Y)
    resid = Y - reg.predict(design)
    n, k = design.shape
    sigma = float(np.sqrt(np.sum(resid ** 2) / max(n - k - 1, 1)))
    # SE of β_T from OLS formula
    XtX_inv = np.linalg.pinv(design.T @ design)
    se = float(sigma * np.sqrt(max(XtX_inv[0, 0], 0.0)))
    z = 1.96
    lo = float(ate_hat - z * se)
    hi = float(ate_hat + z * se)
    rr_pt = approx_rr_from_continuous(ate_hat, sigma)
    # CI bound closer to null on the RR scale
    if ate_hat >= 0:
        rr_ci = approx_rr_from_continuous(max(lo, 0.0), sigma) if lo > 0 else 1.0
    else:
        rr_ci = approx_rr_from_continuous(min(hi, 0.0), sigma) if hi < 0 else 1.0
    return {
        "ate_hat": float(ate_hat),
        "se": se,
        "ci95": [lo, hi],
        "sigma_y": sigma,
        "approx_rr": rr_pt,
        "evalue_point": evalue_from_rr(rr_pt),
        "evalue_ci": evalue_from_rr(rr_ci),
        "null_in_ci": bool(lo <= 0.0 <= hi),
    }


def _matched_pair_diffs(Y, T, X, n_pairs=400, seed=SEED):
    """Build residualized high-T vs low-T matched pairs on X (Mahalanobis-ish).

    Continuous Rosenbaum analysis needs paired outcomes; we pair each high-T
    unit to a nearest low-T neighbor in covariate space.
    """
    rng = np.random.default_rng(int(seed))
    Y = np.asarray(Y, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    X = np.asarray(X, dtype=float)
    Yr = _resid(Y, X)
    Tr = _resid(T, X)
    med = float(np.median(Tr))
    hi = np.where(Tr >= med)[0]
    lo = np.where(Tr < med)[0]
    if len(hi) == 0 or len(lo) == 0:
        return np.zeros(0)
    # subsample for speed
    n_pairs = int(min(n_pairs, len(hi), len(lo)))
    hi = rng.choice(hi, size=n_pairs, replace=False)
    # standardize X for nearest neighbor
    sd = X.std(0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    Xs = (X - X.mean(0)) / sd
    diffs = []
    used = set()
    for i in hi:
        d2 = np.sum((Xs[lo] - Xs[i]) ** 2, axis=1)
        order = np.argsort(d2)
        for j in order:
            jj = int(lo[j])
            if jj not in used:
                used.add(jj)
                diffs.append(float(Yr[i] - Yr[jj]))
                break
    return np.asarray(diffs, dtype=float)


def rosenbaum_sign_pvalue(diffs, gamma):
    """Rosenbaum upper-bound p-value for one-sided positive effect (sign test).

    Under Γ, P(positive pair) ≤ p+ = Γ/(1+Γ). Bound = P(Bin(n,p+) ≥ #positives).
    """
    diffs = np.asarray(diffs, dtype=float)
    diffs = diffs[np.isfinite(diffs) & (np.abs(diffs) > 1e-15)]
    n = len(diffs)
    if n == 0:
        return 1.0
    k = int(np.sum(diffs > 0))
    p_plus = float(gamma) / (1.0 + float(gamma))
    # survival function of Binomial
    return float(stats.binom.sf(k - 1, n, p_plus))


def rosenbaum_analysis(Y, T, X, gammas=None, n_pairs=400, seed=SEED, alpha=0.05):
    """Rosenbaum Γ curve for the naive high-vs-low treatment contrast.

    Returns p-value upper bounds vs Γ and the critical Γ* where the bound
    first exceeds ``alpha`` (effect becomes insignificant under that bias).
    """
    gammas = list(gammas or ROSENBAUM_GAMMAS)
    diffs = _matched_pair_diffs(Y, T, X, n_pairs=n_pairs, seed=seed)
    mean_diff = float(np.mean(diffs)) if len(diffs) else float("nan")
    curve = []
    gamma_star = None
    for g in gammas:
        p_up = rosenbaum_sign_pvalue(diffs, g)
        curve.append({"gamma": float(g), "p_upper": p_up})
        if gamma_star is None and p_up > float(alpha):
            gamma_star = float(g)
    robust_beyond = gamma_star is None
    if gamma_star is None and curve:
        # still significant at max Γ — report as lower bound on fragility
        gamma_star = float(curve[-1]["gamma"])
    return {
        "n_pairs": int(len(diffs)),
        "mean_pair_diff": mean_diff,
        "alpha": float(alpha),
        "gamma_star": float(gamma_star) if gamma_star is not None else float("nan"),
        "robust_beyond_max_gamma": bool(robust_beyond),
        "max_gamma_tried": float(gammas[-1]) if gammas else float("nan"),
        "curve": curve,
        "p_at_gamma1": float(curve[0]["p_upper"]) if curve else float("nan"),
        "p_at_max_gamma": float(curve[-1]["p_upper"]) if curve else float("nan"),
    }


def sensitivity_bundle(bundle, z_mode="auto", seed=SEED):
    """E-value on naive association + Rosenbaum Γ + IV point for contrast."""
    scores = score_instruments(bundle)
    naive = ols_ate(bundle.Y, bundle.T, bundle.X)
    ev = evalue_analysis(bundle.Y, bundle.T, bundle.X, ate_hat=naive)
    rb = rosenbaum_analysis(bundle.Y, bundle.T, bundle.X, seed=seed)
    idx = select_instruments(scores, mode=z_mode)
    Zs = bundle.Z[:, idx] if idx else np.zeros((len(bundle.T), 0))
    iv_hat, iv_info = iv_2sls_ate(bundle.Y, bundle.T, bundle.X, Zs)
    return {
        "naive_ate": float(naive),
        "iv_ate": float(iv_hat),
        "true_ate": float(bundle.true_ate),
        "evalue": ev,
        "rosenbaum": rb,
        "iv_info": iv_info,
        "selected": idx,
    }


# ---------------------------------------------------------------------------
# Subset localization — IV strength heterogeneity across segments
# ---------------------------------------------------------------------------


def subset_iv_localization(bundle, z_mode="oracle", min_n=80):
    """Per-segment first-stage F / R² / 2SLS ATE (IV strength heterogeneity).

    Segments = quartiles of X0 (customer segment proxy). Under the DGP,
    high-X0 segments have stronger valid instruments.
    """
    seg = np.asarray(bundle.z_meta["segment"], dtype=int)
    names = bundle.z_meta["segment_names"]
    scores = score_instruments(bundle)
    idx = select_instruments(scores, mode=z_mode)
    if not idx:
        idx = select_instruments(scores, mode="oracle")
    rows = []
    for s, name in enumerate(names):
        mask = seg == s
        n_s = int(mask.sum())
        if n_s < int(min_n):
            rows.append(
                {
                    "segment": int(s),
                    "name": name,
                    "n": n_s,
                    "stage1_f": float("nan"),
                    "stage1_r2": float("nan"),
                    "ate_hat": float("nan"),
                    "abs_err": float("nan"),
                    "mean_seg_scale": float("nan"),
                }
            )
            continue
        X_s = bundle.X[mask]
        T_s = bundle.T[mask]
        Y_s = bundle.Y[mask]
        Z_s = bundle.Z[mask][:, idx]
        ate, info = iv_2sls_ate(Y_s, T_s, X_s, Z_s)
        # also report mean designed strength scaler
        scale = np.asarray(bundle.z_meta["seg_scale"], dtype=float)[mask]
        rows.append(
            {
                "segment": int(s),
                "name": name,
                "n": n_s,
                "stage1_f": float(info.get("stage1_f", float("nan"))),
                "stage1_r2": float(info.get("stage1_r2", float("nan"))),
                "ate_hat": float(ate),
                "abs_err": float(abs(ate - bundle.true_ate)),
                "mean_seg_scale": float(np.mean(scale)),
                "n_z": int(info.get("n_z", 0)),
            }
        )
    fs = np.array([r["stage1_f"] for r in rows], dtype=float)
    finite = fs[np.isfinite(fs)]
    het_ratio = (
        float(np.nanmax(finite) / max(np.nanmin(finite), 1e-8))
        if finite.size >= 2
        else float("nan")
    )
    return {
        "z_mode": z_mode,
        "segments": rows,
        "f_heterogeneity_ratio": het_ratio,
        "f_max": float(np.nanmax(finite)) if finite.size else float("nan"),
        "f_min": float(np.nanmin(finite)) if finite.size else float("nan"),
    }


def run_sensitivity_suite(seeds=None, n=2000, conf_strength=1.0):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 8))
    sens_rows, subset_rows = [], []
    for s in seeds:
        bundle = make_iv_dgp(n=n, conf_strength=conf_strength, seed=int(s))
        sens = sensitivity_bundle(bundle, z_mode="auto", seed=int(s))
        sens_rows.append(
            {
                "seed": int(s),
                "naive_ate": sens["naive_ate"],
                "iv_ate": sens["iv_ate"],
                "true_ate": sens["true_ate"],
                "evalue_point": sens["evalue"]["evalue_point"],
                "evalue_ci": sens["evalue"]["evalue_ci"],
                "approx_rr": sens["evalue"]["approx_rr"],
                "gamma_star": sens["rosenbaum"]["gamma_star"],
                "robust_beyond_max_gamma": int(
                    sens["rosenbaum"]["robust_beyond_max_gamma"]
                ),
                "rosenbaum_p_gamma1": sens["rosenbaum"]["p_at_gamma1"],
                "rosenbaum_p_max": sens["rosenbaum"]["p_at_max_gamma"],
                "rosenbaum_curve": sens["rosenbaum"]["curve"],
            }
        )
        loc = subset_iv_localization(bundle, z_mode="oracle")
        for seg in loc["segments"]:
            subset_rows.append({"seed": int(s), **seg})
    # aggregate sensitivity
    def _agg(key):
        v = np.array([r[key] for r in sens_rows], dtype=float)
        return {
            "mean": float(np.nanmean(v)),
            "sd": float(np.nanstd(v, ddof=1) if len(v) > 1 else 0.0),
        }

    # mean Rosenbaum curve
    gammas = [c["gamma"] for c in sens_rows[0]["rosenbaum_curve"]]
    mean_curve = []
    for i, g in enumerate(gammas):
        ps = [r["rosenbaum_curve"][i]["p_upper"] for r in sens_rows]
        mean_curve.append({"gamma": float(g), "p_upper_mean": float(np.mean(ps))})

    # aggregate subset by segment name
    seg_table = {}
    for name in ("Q1_low", "Q2", "Q3", "Q4_high"):
        sub = [r for r in subset_rows if r["name"] == name]
        seg_table[name] = {
            "n_mean": float(np.mean([r["n"] for r in sub])),
            "stage1_f": _mstd([r["stage1_f"] for r in sub]),
            "stage1_r2": _mstd([r["stage1_r2"] for r in sub]),
            "ate_hat": _mstd([r["ate_hat"] for r in sub]),
            "abs_err": _mstd([r["abs_err"] for r in sub]),
            "mean_seg_scale": _mstd([r["mean_seg_scale"] for r in sub]),
        }
    f_by_seg = [seg_table[n]["stage1_f"]["mean"] for n in ("Q1_low", "Q2", "Q3", "Q4_high")]
    return {
        "seeds": seeds,
        "n": int(n),
        "conf_strength": float(conf_strength),
        "sensitivity": {
            "evalue_point": _agg("evalue_point"),
            "evalue_ci": _agg("evalue_ci"),
            "approx_rr": _agg("approx_rr"),
            "gamma_star": _agg("gamma_star"),
            "robust_beyond_frac": float(
                np.mean([r["robust_beyond_max_gamma"] for r in sens_rows])
            ),
            "rosenbaum_p_max": _agg("rosenbaum_p_max"),
            "naive_ate": _agg("naive_ate"),
            "iv_ate": _agg("iv_ate"),
            "rosenbaum_curve_mean": mean_curve,
            "rows": sens_rows,
        },
        "subset": {
            "by_segment": seg_table,
            "f_heterogeneity_ratio": float(
                max(f_by_seg) / max(min(f_by_seg), 1e-8)
            ),
            "rows": subset_rows,
        },
    }


def _mstd(vals):
    v = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.nanmean(v)),
        "sd": float(np.nanstd(v, ddof=1) if np.sum(np.isfinite(v)) > 1 else 0.0),
    }


# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------


METHODS = (
    "naive_ols",
    "naive_partial",
    "iv_oracle",
    "iv_auto",       # strength + exclusion selection
    "iv_strong",     # F-screen only (may keep invalid if first-stage leaks)
    "iv_all",
    "iv_invalid",    # ablation: exclusion-broken instruments only
    "iv_weak",       # ablation: weak instruments only
)


def run_once(bundle, scores=None):
    scores = scores if scores is not None else score_instruments(bundle)
    out = {}
    out["naive_ols"] = {
        "ate_hat": ols_ate(bundle.Y, bundle.T, bundle.X),
        "n_z": 0,
        "stage1_f": float("nan"),
        "selected": [],
    }
    out["naive_partial"] = {
        "ate_hat": tlearner_ate(bundle.Y, bundle.T, bundle.X),
        "n_z": 0,
        "stage1_f": float("nan"),
        "selected": [],
    }
    for mode, key in (
        ("oracle", "iv_oracle"),
        ("auto", "iv_auto"),
        ("strong", "iv_strong"),
        ("all", "iv_all"),
        ("invalid", "iv_invalid"),
        ("weak", "iv_weak"),
    ):
        idx = select_instruments(scores, mode=mode)
        Zs = bundle.Z[:, idx] if idx else np.zeros((len(bundle.T), 0))
        ate, info = iv_2sls_ate(bundle.Y, bundle.T, bundle.X, Zs)
        roles = [bundle.z_meta["roles"][j] for j in idx]
        out[key] = {
            "ate_hat": ate,
            "n_z": int(info["n_z"]),
            "stage1_f": float(info.get("stage1_f", float("nan"))),
            "stage1_r2": float(info.get("stage1_r2", float("nan"))),
            "selected": idx,
            "selected_roles": roles,
            "n_valid_selected": int(sum(r == "valid" for r in roles)),
            "n_invalid_selected": int(sum(r == "invalid" for r in roles)),
            "n_weak_selected": int(sum(r == "weak" for r in roles)),
        }
    return out, scores


def run_suite(
    seeds=None,
    n=2000,
    conf_strength=1.0,
    methods=None,
):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 8))
    methods = list(methods or METHODS)
    rows = []
    score_book = []
    for s in seeds:
        bundle = make_iv_dgp(n=n, conf_strength=conf_strength, seed=int(s))
        est, scores = run_once(bundle)
        score_book.append({"seed": int(s), "scores": scores})
        for m in methods:
            cell = est[m]
            err = float(cell["ate_hat"] - bundle.true_ate)
            rows.append(
                {
                    "seed": int(s),
                    "method": m,
                    "ate_hat": float(cell["ate_hat"]),
                    "true_ate": float(bundle.true_ate),
                    "abs_err": abs(err),
                    "err": err,
                    "n_z": int(cell.get("n_z", 0)),
                    "stage1_f": float(cell.get("stage1_f", float("nan"))),
                    "n_valid_selected": int(cell.get("n_valid_selected", 0)),
                    "n_invalid_selected": int(cell.get("n_invalid_selected", 0)),
                    "n_weak_selected": int(cell.get("n_weak_selected", 0)),
                }
            )
    table = {}
    for m in methods:
        sub = [r for r in rows if r["method"] == m]
        ae = np.array([r["abs_err"] for r in sub], dtype=float)
        hats = np.array([r["ate_hat"] for r in sub], dtype=float)
        table[m] = {
            "abs_err": {
                "mean": float(np.nanmean(ae)),
                "sd": float(np.nanstd(ae, ddof=1) if len(ae) > 1 else 0.0),
            },
            "ate_hat": {
                "mean": float(np.nanmean(hats)),
                "sd": float(np.nanstd(hats, ddof=1) if len(hats) > 1 else 0.0),
            },
            "mean_n_z": float(np.nanmean([r["n_z"] for r in sub])),
            "mean_n_valid": float(np.nanmean([r["n_valid_selected"] for r in sub])),
            "mean_n_invalid": float(np.nanmean([r["n_invalid_selected"] for r in sub])),
            "n": int(len(sub)),
        }
    # selection quality: precision/recall of auto vs oracle valid set
    sel_stats = _selection_quality(score_book)
    return {
        "seeds": seeds,
        "methods": methods,
        "true_ate": float(make_iv_dgp(seed=SEED).true_ate),
        "conf_strength": float(conf_strength),
        "n": int(n),
        "table": table,
        "rows": rows,
        "selection": sel_stats,
        "score_book": score_book,
    }


def _selection_quality(score_book):
    precs, recs, f1s = [], [], []
    for entry in score_book:
        scores = entry["scores"]
        truth = {r["j"] for r in scores if r["role"] == "valid"}
        pred = {r["j"] for r in scores if r["selected"]}
        tp = len(truth & pred)
        prec = tp / max(len(pred), 1)
        rec = tp / max(len(truth), 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-12)
        precs.append(prec)
        recs.append(rec)
        f1s.append(f1)
    return {
        "precision": {"mean": float(np.mean(precs)), "sd": float(np.std(precs, ddof=1) if len(precs) > 1 else 0.0)},
        "recall": {"mean": float(np.mean(recs)), "sd": float(np.std(recs, ddof=1) if len(recs) > 1 else 0.0)},
        "f1": {"mean": float(np.mean(f1s)), "sd": float(np.std(f1s, ddof=1) if len(f1s) > 1 else 0.0)},
    }


def conf_sweep(conf_grid=None, seeds=None, n=2000):
    """Bias vs confounding strength: naive blows up, IV-auto stays near τ."""
    conf_grid = list(conf_grid or [0.0, 0.5, 1.0, 1.5, 2.0])
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 5))
    out = []
    for c in conf_grid:
        suite = run_suite(seeds=seeds, n=n, conf_strength=float(c))
        out.append(
            {
                "conf_strength": float(c),
                "naive_ols_abs_err": suite["table"]["naive_ols"]["abs_err"]["mean"],
                "iv_oracle_abs_err": suite["table"]["iv_oracle"]["abs_err"]["mean"],
                "iv_auto_abs_err": suite["table"]["iv_auto"]["abs_err"]["mean"],
                "iv_invalid_abs_err": suite["table"]["iv_invalid"]["abs_err"]["mean"],
                "selection_f1": suite["selection"]["f1"]["mean"],
            }
        )
    return out


def plot_suite(suite, sweep, path, sens=None):
    import matplotlib.pyplot as plt

    methods = suite["methods"]
    colors = {
        "naive_ols": "#9AA3AE",
        "naive_partial": "#6B7280",
        "iv_oracle": "#2F6B4F",
        "iv_auto": "#C45C26",
        "iv_strong": "#C4892A",
        "iv_all": "#5B3A8C",
        "iv_invalid": "#A33B24",
        "iv_weak": "#7A4E2D",
    }
    n_panels = 4 if sens is not None else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(3.4 * n_panels, 4.2), squeeze=False)
    axes = axes[0]

    ax = axes[0]
    means = [suite["table"][m]["abs_err"]["mean"] for m in methods]
    sds = [suite["table"][m]["abs_err"]["sd"] for m in methods]
    ypos = np.arange(len(methods))[::-1]
    ax.barh(
        ypos, means, xerr=sds,
        color=[colors.get(m, "#333") for m in methods], height=0.62,
        error_kw={"ecolor": "#555", "lw": 0.9, "capsize": 2},
    )
    ax.axvline(0.0, color="#222", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels(methods, fontsize=8)
    ax.set_xlabel("|ATÊ − τ| ↓")
    ax.set_title("IV absolute error", loc="left", fontsize=10, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.35)

    ax = axes[1]
    cs = [r["conf_strength"] for r in sweep]
    ax.plot(cs, [r["naive_ols_abs_err"] for r in sweep], "o-", color="#9AA3AE", label="naive_ols")
    ax.plot(cs, [r["iv_oracle_abs_err"] for r in sweep], "s-", color="#2F6B4F", label="iv_oracle")
    ax.plot(cs, [r["iv_auto_abs_err"] for r in sweep], "D-", color="#C45C26", label="iv_auto")
    ax.plot(cs, [r["iv_invalid_abs_err"] for r in sweep], "x--", color="#A33B24", label="iv_invalid")
    ax.set_xlabel("confounding strength α")
    ax.set_ylabel("|ATÊ − τ| ↓")
    ax.set_title("Error vs confounding", loc="left", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(True, alpha=0.35)

    if sens is not None:
        ax = axes[2]
        curve = sens["sensitivity"]["rosenbaum_curve_mean"]
        ax.plot(
            [c["gamma"] for c in curve],
            [c["p_upper_mean"] for c in curve],
            "o-", color="#5B3A8C", lw=1.6,
        )
        ax.axhline(0.05, color="#A33B24", ls="--", lw=0.9, label="α=0.05")
        gstar = sens["sensitivity"]["gamma_star"]["mean"]
        ax.axvline(gstar, color="#C45C26", ls=":", lw=1.2, label="Γ*≈%.1f" % gstar)
        ax.set_xlabel("Rosenbaum Γ")
        ax.set_ylabel("p-value upper bound")
        ax.set_title("Rosenbaum sensitivity", loc="left", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7.5, frameon=False)
        ax.grid(True, alpha=0.35)

        ax = axes[3]
        names = ["Q1_low", "Q2", "Q3", "Q4_high"]
        fs = [sens["subset"]["by_segment"][n]["stage1_f"]["mean"] for n in names]
        fe = [sens["subset"]["by_segment"][n]["stage1_f"]["sd"] for n in names]
        xpos = np.arange(len(names))
        ax.bar(xpos, fs, yerr=fe, color="#1F6F8B", width=0.66, capsize=3, ecolor="#555")
        ax.axhline(F_WEAK, color="#A33B24", ls="--", lw=0.9, label="F=10")
        ax.set_xticks(xpos)
        ax.set_xticklabels(names, fontsize=8)
        ax.set_ylabel("first-stage F ↑")
        ax.set_title(
            "Subset IV strength (het≈%.1f×)" % sens["subset"]["f_heterogeneity_ratio"],
            loc="left", fontsize=10, fontweight="bold",
        )
        ax.legend(fontsize=7.5, frameon=False)
        ax.grid(True, axis="y", alpha=0.35)

    tau = suite["true_ate"]
    ev = ""
    if sens is not None:
        ev = "  |  E-value≈%.2f  Γ*≈%.1f" % (
            sens["sensitivity"]["evalue_point"]["mean"],
            sens["sensitivity"]["gamma_star"]["mean"],
        )
    fig.suptitle(
        "Instrumental-variable selection + sensitivity (τ=%.2f)%s" % (tau, ev),
        fontsize=12, fontweight="bold", y=0.98,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def write_readme(suite, sweep, path, sens=None):
    lines = [
        "# Instrumental Variable selection — synthetic argument",
        "",
        "Hidden confounder `U → (T, Y)`. Naive OLS / partialling is biased.",
        "2SLS with **valid** instruments recovers `τ`; invalid / weak Z do not.",
        "",
        "Selection rule (`iv_auto`): first-stage partial **F ≥ 10** and",
        "exclusion excess screen.",
        "",
        "## Absolute error `|ATÊ − τ|` (lower better)",
        "",
        "| method | abs_err | ATÊ | n_Z | valid / invalid picked |",
        "| --- | --- | --- | --- | --- |",
    ]
    for m in suite["methods"]:
        t = suite["table"][m]
        lines.append(
            "| %s | %.3f (%.3f) | %.3f | %.1f | %.1f / %.1f |"
            % (
                m,
                t["abs_err"]["mean"],
                t["abs_err"]["sd"],
                t["ate_hat"]["mean"],
                t["mean_n_z"],
                t["mean_n_valid"],
                t["mean_n_invalid"],
            )
        )
    sel = suite["selection"]
    lines.extend(
        [
            "",
            "## IV selection quality (`iv_auto` vs true valid set)",
            "",
            "| precision | recall | F1 |",
            "| --- | --- | --- |",
            "| %.3f | %.3f | %.3f |"
            % (sel["precision"]["mean"], sel["recall"]["mean"], sel["f1"]["mean"]),
            "",
            "## Confounding sweep",
            "",
            "| α | naive | iv_oracle | iv_auto | iv_invalid | sel F1 |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for r in sweep:
        lines.append(
            "| %.1f | %.3f | %.3f | %.3f | %.3f | %.3f |"
            % (
                r["conf_strength"],
                r["naive_ols_abs_err"],
                r["iv_oracle_abs_err"],
                r["iv_auto_abs_err"],
                r["iv_invalid_abs_err"],
                r["selection_f1"],
            )
        )
    if sens is not None:
        S = sens["sensitivity"]
        lines.extend(
            [
                "",
                "## Sensitivity on the naive association",
                "",
                "| metric | mean (sd) |",
                "| --- | --- |",
                "| E-value (point) | %.3f (%.3f) |"
                % (S["evalue_point"]["mean"], S["evalue_point"]["sd"]),
                "| E-value (CI) | %.3f (%.3f) |"
                % (S["evalue_ci"]["mean"], S["evalue_ci"]["sd"]),
                "| approx RR | %.3f (%.3f) |"
                % (S["approx_rr"]["mean"], S["approx_rr"]["sd"]),
                "| Rosenbaum Γ* (p-bound > 0.05) | %.2f (%.2f) |"
                % (S["gamma_star"]["mean"], S["gamma_star"]["sd"]),
                "| naive ATÊ | %.3f (%.3f) |"
                % (S["naive_ate"]["mean"], S["naive_ate"]["sd"]),
                "| iv ATÊ | %.3f (%.3f) |"
                % (S["iv_ate"]["mean"], S["iv_ate"]["sd"]),
                "",
                "E-value: min strength of an unmeasured confounder (on the RR scale)",
                "needed to explain away the naive effect. Rosenbaum Γ*: odds-ratio",
                "bias in treatment assignment at which the matched-pair sign test",
                "upper-bounds to insignificance.",
                "",
                "## Subset localization — IV strength heterogeneity",
                "",
                "Segments = quartiles of `X0` (designed stronger first stage in Q4).",
                "",
                "| segment | n | stage-1 F | R² | ATÊ | abs_err | scale |",
                "| --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for name in ("Q1_low", "Q2", "Q3", "Q4_high"):
            t = sens["subset"]["by_segment"][name]
            lines.append(
                "| %s | %.0f | %.1f (%.1f) | %.3f | %.3f | %.3f | %.2f |"
                % (
                    name,
                    t["n_mean"],
                    t["stage1_f"]["mean"],
                    t["stage1_f"]["sd"],
                    t["stage1_r2"]["mean"],
                    t["ate_hat"]["mean"],
                    t["abs_err"]["mean"],
                    t["mean_seg_scale"]["mean"],
                )
            )
        lines.append("")
        lines.append(
            "**F heterogeneity ratio** (max/min across segments): **%.2f×**."
            % sens["subset"]["f_heterogeneity_ratio"]
        )
    lines.extend(
        [
            "",
            "## Run",
            "",
            "```bash",
            "python3 scripts/run_iv_selection_synth.py",
            "python3 scripts/run_iv_selection_synth.py --quick",
            "```",
            "",
            "Code: `Python/src/iv_selection_synth.py`",
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

    path.write_text(json.dumps(obj, indent=2, default=_c))
    return path
