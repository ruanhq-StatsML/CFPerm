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

Locked takeaway for the FSDS product case:
  when treatment is confounded, uplift monitors need an IV channel;
  selection must prefer strong + exclusion-plausible Z, not raw correlation.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression

SEED = 2026
F_WEAK = 10.0  # Staiger–Stock rule-of-thumb for first-stage strength


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
    seed=SEED,
):
    """Linear IV world with one hidden confounder.

    T = X β_t + Z_valid γ_v + Z_weak γ_w + α_t U + ε_t
    Y = τ T + X β_y + Z_invalid δ + α_y U + ε_y

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

    eps_t = rng.normal(scale=0.5, size=n)
    eps_y = rng.normal(scale=0.5, size=n)
    T = (
        X @ beta_t
        + Zv @ gamma_v
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
        },
        meta={
            "n": n,
            "p_x": p_x,
            "conf_strength": float(conf_strength),
            "valid_strength": float(valid_strength),
            "invalid_y": float(invalid_y),
            "weak_strength": float(weak_strength),
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


def plot_suite(suite, sweep, path):
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
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

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
    ax.set_title("Synthetic IV: absolute error", loc="left", fontsize=10, fontweight="bold")
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

    tau = suite["true_ate"]
    fig.suptitle(
        "Instrumental-variable selection on synthetic uplift (τ=%.2f)" % tau,
        fontsize=12, fontweight="bold", y=0.98,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def write_readme(suite, sweep, path):
    lines = [
        "# Instrumental Variable selection — synthetic argument",
        "",
        "Hidden confounder `U → (T, Y)`. Naive OLS / partialling is biased.",
        "2SLS with **valid** instruments recovers `τ`; invalid / weak Z do not.",
        "",
        "Selection rule (`iv_auto`): first-stage partial **F ≥ 10** and",
        "exclusion residual score ≤ 0.12.",
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
