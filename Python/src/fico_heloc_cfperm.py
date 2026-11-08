"""FICO HELOC CFPerm experiment: LOCO / PermuCATE / GRF VIMPs + subset localization.

Design (matches the notebook sketch):
  1. Load HELOC, drop rows with special code -9.
  2. Fit logistic RiskPerformance ~ X; take top-3 |coef| features as split axes.
  3. For each split feature, median-partition into exist (W=0) vs new (W=1).
  4. Run CFPerm (permute W, threshold on VIMP) with B trees and n_perm batch perms.
  5. Subset-localize mean differences on rejected / top VIMP features; summarize power.

Power on this real-data design = rejection under the constructed alternative + recovery
of the true median-split feature among rejected / top-ranked VIMPs.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

HELOC_URL = (
    "https://raw.githubusercontent.com/benoitparis/explainable-challenge/"
    "master/heloc_dataset_v1.csv"
)
VIMP_NOTIONS = ("loco", "permucate", "grf")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _as_1d(x) -> np.ndarray:
    return np.asarray(x).reshape(-1)


def _safe_norm(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    s = float(np.sum(np.abs(v)))
    if not np.isfinite(s) or s <= 0:
        return np.zeros_like(v, dtype=float)
    return np.abs(v) / s


def _encode_y(y: pd.Series | np.ndarray) -> np.ndarray:
    """Map RiskPerformance Bad/Good -> {0,1} (Good=1)."""
    s = pd.Series(y)
    if pd.api.types.is_numeric_dtype(s):
        return _as_1d(s).astype(float)
    mapped = s.astype(str).str.strip().str.lower().map({"bad": 0, "good": 1})
    if mapped.isna().any():
        cats = {c: i for i, c in enumerate(pd.unique(s.astype(str)))}
        return np.asarray([cats[v] for v in s.astype(str)], dtype=float)
    return mapped.astype(float).to_numpy()


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_heloc(url: str = HELOC_URL, cache_path: Optional[Path] = None) -> pd.DataFrame:
    if cache_path is not None and Path(cache_path).exists():
        return pd.read_csv(cache_path)
    df = pd.read_csv(url)
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache_path, index=False)
    return df


def clean_heloc(df: pd.DataFrame, missing_code: int = -9) -> pd.DataFrame:
    mask = ~df.isin([missing_code]).any(axis=1)
    return df.loc[mask].copy().reset_index(drop=True)


def top_lr_features(
    df: pd.DataFrame,
    y_col: str = "RiskPerformance",
    top_k: int = 3,
    max_iter: int = 2000,
    seed: int = 0,
) -> List[str]:
    feature_cols = [c for c in df.columns if c != y_col]
    X = df[feature_cols].to_numpy(dtype=float)
    y = _encode_y(df[y_col])
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "lr",
                LogisticRegression(max_iter=max_iter, random_state=seed, solver="lbfgs"),
            ),
        ]
    )
    pipe.fit(X, y)
    coef = np.abs(pipe.named_steps["lr"].coef_).ravel()
    order = np.argsort(-coef)[:top_k]
    return [feature_cols[i] for i in order]


def median_split_batches(
    df: pd.DataFrame,
    feature: str,
    y_col: str = "RiskPerformance",
    max_n_per_arm: Optional[int] = None,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    """Split on feature median → (X, Y, W, feature_names). W=0 exist, W=1 new."""
    feature_cols = [c for c in df.columns if c != y_col]
    med = float(df[feature].median())
    exist = df.loc[df[feature] < med]
    new = df.loc[df[feature] >= med]
    rng = np.random.default_rng(seed)
    if max_n_per_arm is not None:
        if len(exist) > max_n_per_arm:
            idx = rng.choice(len(exist), size=max_n_per_arm, replace=False)
            exist = exist.iloc[idx]
        if len(new) > max_n_per_arm:
            idx = rng.choice(len(new), size=max_n_per_arm, replace=False)
            new = new.iloc[idx]
    exist = exist.reset_index(drop=True)
    new = new.reset_index(drop=True)
    X = np.vstack(
        [
            exist[feature_cols].to_numpy(dtype=float),
            new[feature_cols].to_numpy(dtype=float),
        ]
    )
    Y = np.concatenate([_encode_y(exist[y_col]), _encode_y(new[y_col])])
    W = np.concatenate(
        [np.zeros(len(exist), dtype=int), np.ones(len(new), dtype=int)]
    )
    return X, Y, W, feature_cols


# ---------------------------------------------------------------------------
# VIMP: LOCO R-risk
# ---------------------------------------------------------------------------
def _crossfit_nuisances(
    X, Y, W, *, seed=0, n_splits=3, clip_e=0.01, n_estimators=120
):
    X = np.asarray(X, dtype=float)
    Y = _as_1d(Y).astype(float)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    for tr, te in cv.split(X, W):
        rf = RandomForestRegressor(
            n_estimators=n_estimators,
            min_samples_leaf=5,
            max_features=max(1, int(round(p / 3))),
            random_state=seed,
            n_jobs=1,
        )
        rf.fit(X[tr], Y[tr])
        m_hat[te] = rf.predict(X[te])
        lr = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed),
                ),
            ]
        )
        lr.fit(X[tr], W[tr])
        e_hat[te] = lr.predict_proba(X[te])[:, 1]
    return m_hat, np.clip(e_hat, clip_e, 1.0 - clip_e)


def _fit_tau_ridge(X, Y_tilde, W_tilde, *, seed=0, clip_wtilde=1e-3, alpha=1.0):
    mask = np.abs(W_tilde) > clip_wtilde
    if int(mask.sum()) < max(30, X.shape[1] + 5):

        class _Z:
            def predict(self, X_):
                return np.zeros(X_.shape[0])

        return _Z()
    z = Y_tilde[mask] / W_tilde[mask]
    w = W_tilde[mask] ** 2
    model = Pipeline(
        [("scaler", StandardScaler()), ("ridge", Ridge(alpha=alpha, random_state=seed))]
    )
    model.fit(X[mask], z, ridge__sample_weight=w)
    return model


def _r_risk(Y_tilde, W_tilde, tau_hat) -> float:
    return float(np.mean((_as_1d(Y_tilde) - _as_1d(W_tilde) * _as_1d(tau_hat)) ** 2))


def vimp_loco(
    X, Y, W, *, seed=0, n_splits=3, n_estimators=120, clip_e=0.01, clip_wtilde=1e-3
) -> np.ndarray:
    m_hat, e_hat = _crossfit_nuisances(
        X, Y, W, seed=seed, n_splits=n_splits, clip_e=clip_e, n_estimators=n_estimators
    )
    Y_tilde = _as_1d(Y) - m_hat
    W_tilde = _as_1d(W).astype(float) - e_hat
    tau_full = _fit_tau_ridge(X, Y_tilde, W_tilde, seed=seed, clip_wtilde=clip_wtilde)
    risk_full = _r_risk(Y_tilde, W_tilde, tau_full.predict(X))
    p = X.shape[1]
    imps = np.zeros(p)
    for j in range(p):
        Xj = np.delete(X, j, axis=1)
        tau_j = _fit_tau_ridge(
            Xj, Y_tilde, W_tilde, seed=seed + 17 * (j + 1), clip_wtilde=clip_wtilde
        )
        imps[j] = _r_risk(Y_tilde, W_tilde, tau_j.predict(Xj)) - risk_full
    return _safe_norm(imps)


# ---------------------------------------------------------------------------
# VIMP: GRF (CausalForestDML impurity / split importance)
# ---------------------------------------------------------------------------
def vimp_grf(
    X, Y, W, *, seed=0, n_estimators=150, min_samples_leaf=5, max_depth=None
) -> np.ndarray:
    from econml.dml import CausalForestDML

    X = np.asarray(X, dtype=float)
    Y = _as_1d(Y).astype(float)
    W = _as_1d(W).astype(int)
    p = X.shape[1]
    # econml CausalForest requires n_estimators divisible by subforest_size (4).
    n_estimators = int(max(4, ((int(n_estimators) + 3) // 4) * 4))
    model_y = RandomForestRegressor(
        n_estimators=min(120, n_estimators),
        min_samples_leaf=min_samples_leaf,
        max_features=max(1, int(round(p / 3))),
        random_state=seed,
        n_jobs=1,
    )
    model_t = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed),
            ),
        ]
    )
    cf = CausalForestDML(
        model_y=model_y,
        model_t=model_t,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        max_depth=max_depth,
        discrete_treatment=True,
        cv=3,
        random_state=seed,
    )
    cf.fit(Y=Y, T=W, X=X)
    imp = np.asarray(cf.feature_importances_, dtype=float).ravel()
    if imp.size < p:
        imp = np.pad(imp, (0, p - imp.size))
    return _safe_norm(imp[:p])


# ---------------------------------------------------------------------------
# VIMP: PermuCATE (conditional residual permutation on PO-risk)
# ---------------------------------------------------------------------------
def _po_risk(tau, Y, W, mu0, mu1, pi, clip=1e-3) -> float:
    pi_c = np.clip(pi, clip, 1.0 - clip)
    mu_t = np.where(W == 1, mu1, mu0)
    phi = (Y - mu_t) * (W - pi_c) / (pi_c * (1.0 - pi_c)) + (mu1 - mu0)
    return float(np.mean((phi - tau) ** 2))


def vimp_permucate(
    X,
    Y,
    W,
    *,
    seed=0,
    n_perm_cate: int = 50,
    n_estimators: int = 120,
    test_size: float = 0.4,
    clip_e: float = 0.01,
) -> np.ndarray:
    """Conditional residual permutation importance for CATE / PO-risk."""
    X = np.asarray(X, dtype=float)
    Y = _as_1d(Y).astype(float)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    Xtr, Xte, Ytr, Yte, Wtr, Wte = train_test_split(
        X, Y, W, test_size=test_size, stratify=W, random_state=seed
    )
    # PermuCATE nuisances: RandomForest outcome models (mu0/mu1/tau) + RF for X_j|X_{-j}
    rf0 = RandomForestRegressor(
        n_estimators=n_estimators, min_samples_leaf=5, random_state=seed, n_jobs=1
    )
    rf1 = RandomForestRegressor(
        n_estimators=n_estimators, min_samples_leaf=5, random_state=seed + 1, n_jobs=1
    )
    m0 = Wtr == 0
    m1 = Wtr == 1
    rf0.fit(Xtr[m0], Ytr[m0])
    rf1.fit(Xtr[m1], Ytr[m1])
    lr = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed),
            ),
        ]
    )
    lr.fit(Xtr, Wtr)
    mu0_te = rf0.predict(Xte)
    mu1_te = rf1.predict(Xte)
    pi_tr = np.clip(lr.predict_proba(Xtr)[:, 1], clip_e, 1.0 - clip_e)
    pi_te = np.clip(lr.predict_proba(Xte)[:, 1], clip_e, 1.0 - clip_e)
    mu0_tr = rf0.predict(Xtr)
    mu1_tr = rf1.predict(Xtr)
    phi_tr = (Ytr - np.where(Wtr == 1, mu1_tr, mu0_tr)) * (Wtr - pi_tr) / (
        pi_tr * (1.0 - pi_tr)
    ) + (mu1_tr - mu0_tr)
    tau_model = RandomForestRegressor(
        n_estimators=n_estimators, min_samples_leaf=5, random_state=seed + 2, n_jobs=1
    )
    tau_model.fit(Xtr, phi_tr)
    tau_te = tau_model.predict(Xte)
    base = _po_risk(tau_te, Yte, Wte, mu0_te, mu1_te, pi_te, clip=clip_e)

    # conditional models X_j | X_{-j} also RandomForest (standard PermuCATE)
    nu_models = []
    for j in range(p):
        Xmj = np.delete(Xtr, j, axis=1)
        nu = RandomForestRegressor(
            n_estimators=n_estimators,
            min_samples_leaf=5,
            random_state=seed + 100 + j,
            n_jobs=1,
        )
        nu.fit(Xmj, Xtr[:, j])
        nu_models.append(nu)

    rng = np.random.default_rng(seed)
    imps = np.zeros(p)
    for j in range(p):
        Xmj_te = np.delete(Xte, j, axis=1)
        nu_hat = nu_models[j].predict(Xmj_te)
        resid = Xte[:, j] - nu_hat
        deltas = np.zeros(n_perm_cate)
        for k in range(n_perm_cate):
            Xp = Xte.copy()
            Xp[:, j] = nu_hat + rng.permutation(resid)
            tau_p = tau_model.predict(Xp)
            deltas[k] = _po_risk(tau_p, Yte, Wte, mu0_te, mu1_te, pi_te, clip=clip_e) - base
        imps[j] = float(np.mean(deltas))
    return _safe_norm(imps)


# ---------------------------------------------------------------------------
# CFPerm permutation test
# ---------------------------------------------------------------------------
@dataclass
class CFPermResult:
    vimp: str
    split_feature: str
    rejected: int
    rejected_features: List[str]
    recovered_split: int
    threshold: float
    imp: Dict[str, float]
    pvals: Dict[str, float]
    top_features: List[str]
    n: int
    n0: int
    n1: int
    n_perm: int
    n_estimators: int


def _compute_vimp(vimp: str, X, Y, W, *, seed, n_estimators, n_perm_cate) -> np.ndarray:
    if vimp == "loco":
        return vimp_loco(X, Y, W, seed=seed, n_estimators=n_estimators)
    if vimp == "grf":
        return vimp_grf(X, Y, W, seed=seed, n_estimators=n_estimators)
    if vimp == "permucate":
        return vimp_permucate(
            X, Y, W, seed=seed, n_estimators=n_estimators, n_perm_cate=n_perm_cate
        )
    raise ValueError(f"unknown vimp={vimp}")


def cfperm_test(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    feature_names: Sequence[str],
    *,
    vimp: str = "loco",
    split_feature: str = "",
    n_perm: int = 200,
    n_estimators: int = 150,
    n_perm_cate: int = 50,
    level_feature: float = 0.99,
    level_across_feature: float = 0.95,
    top_k: int = 1,
    seed: int = 0,
    n_jobs: int = 1,
) -> CFPermResult:
    X = np.asarray(X, dtype=float)
    Y = _as_1d(Y).astype(float)
    W = _as_1d(W).astype(int)
    names = list(feature_names)
    p = X.shape[1]
    imp_obs = _compute_vimp(
        vimp, X, Y, W, seed=seed, n_estimators=n_estimators, n_perm_cate=n_perm_cate
    )

    def _one(b: int) -> np.ndarray:
        rng = np.random.default_rng(seed + 10_000 + b)
        Wb = rng.permutation(W)
        return _compute_vimp(
            vimp,
            X,
            Y,
            Wb,
            seed=seed + 10_000 + b,
            n_estimators=n_estimators,
            n_perm_cate=n_perm_cate,
        )

    if n_jobs == 1 or n_perm <= 1:
        perm = np.column_stack([_one(b) for b in range(n_perm)]) if n_perm else np.zeros((p, 0))
    else:
        parts = Parallel(n_jobs=n_jobs, prefer="processes")(
            delayed(_one)(b) for b in range(n_perm)
        )
        perm = np.column_stack(parts) if parts else np.zeros((p, 0))

    if n_perm > 0:
        pvals = (1.0 + np.sum(perm >= imp_obs[:, None], axis=1)) / (1.0 + n_perm)
        q1 = np.quantile(perm, q=level_feature, axis=1)
        thr = float(np.quantile(q1, q=level_across_feature))
    else:
        pvals = np.ones(p)
        thr = float("inf")

    rejected_idx = [j for j in range(p) if imp_obs[j] > thr]
    rejected = int(len(rejected_idx) >= top_k)
    rejected_names = [names[j] for j in rejected_idx]
    order = np.argsort(-imp_obs)
    top_features = [names[j] for j in order[: min(5, p)]]
    recovered = int(bool(split_feature) and split_feature in rejected_names)
    return CFPermResult(
        vimp=vimp,
        split_feature=split_feature,
        rejected=rejected,
        rejected_features=rejected_names,
        recovered_split=recovered,
        threshold=thr,
        imp={names[j]: float(imp_obs[j]) for j in range(p)},
        pvals={names[j]: float(pvals[j]) for j in range(p)},
        top_features=top_features,
        n=int(X.shape[0]),
        n0=int(np.sum(W == 0)),
        n1=int(np.sum(W == 1)),
        n_perm=int(n_perm),
        n_estimators=int(n_estimators),
    )


# ---------------------------------------------------------------------------
# subset localization
# ---------------------------------------------------------------------------
@dataclass
class SubsetLocalRow:
    split_feature: str
    vimp: str
    localize_feature: str
    bin: str
    n0: int
    n1: int
    mean_y0: float
    mean_y1: float
    delta_y: float
    mean_x0: float
    mean_x1: float
    delta_x: float


def subset_localize(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    feature_names: Sequence[str],
    localize_features: Sequence[str],
    *,
    split_feature: str,
    vimp: str,
    n_bins: int = 4,
) -> List[SubsetLocalRow]:
    """Discretize localize_features into quantile bins; report ΔȲ and ΔX̄ by W."""
    names = list(feature_names)
    name_to_j = {n: i for i, n in enumerate(names)}
    Y = _as_1d(Y).astype(float)
    W = _as_1d(W).astype(int)
    rows: List[SubsetLocalRow] = []
    for feat in localize_features:
        if feat not in name_to_j:
            continue
        j = name_to_j[feat]
        col = X[:, j]
        try:
            bins = pd.qcut(col, q=n_bins, duplicates="drop")
        except ValueError:
            bins = pd.cut(col, bins=min(n_bins, max(2, len(np.unique(col)))), duplicates="drop")
        for lab in bins.categories:
            mask = bins == lab
            if int(mask.sum()) < 10:
                continue
            w0 = mask & (W == 0)
            w1 = mask & (W == 1)
            if int(w0.sum()) < 3 or int(w1.sum()) < 3:
                continue
            my0 = float(np.mean(Y[w0]))
            my1 = float(np.mean(Y[w1]))
            mx0 = float(np.mean(col[w0]))
            mx1 = float(np.mean(col[w1]))
            rows.append(
                SubsetLocalRow(
                    split_feature=split_feature,
                    vimp=vimp,
                    localize_feature=feat,
                    bin=str(lab),
                    n0=int(w0.sum()),
                    n1=int(w1.sum()),
                    mean_y0=my0,
                    mean_y1=my1,
                    delta_y=my1 - my0,
                    mean_x0=mx0,
                    mean_x1=mx1,
                    delta_x=mx1 - mx0,
                )
            )
    return rows


def power_summary(results: Sequence[CFPermResult]) -> Dict[str, object]:
    """Aggregate rejection / recovery power across split×vimp cells."""
    by_vimp: Dict[str, Dict[str, float]] = {}
    for v in VIMP_NOTIONS:
        cells = [r for r in results if r.vimp == v]
        if not cells:
            continue
        by_vimp[v] = {
            "n_cells": float(len(cells)),
            "reject_rate": float(np.mean([r.rejected for r in cells])),
            "recover_rate": float(np.mean([r.recovered_split for r in cells])),
            "top1_hit_rate": float(
                np.mean([1.0 if r.top_features and r.top_features[0] == r.split_feature else 0.0 for r in cells])
            ),
            "top3_hit_rate": float(
                np.mean(
                    [
                        1.0
                        if r.split_feature in r.top_features[:3]
                        else 0.0
                        for r in cells
                    ]
                )
            ),
        }
    overall = {
        "n_cells": float(len(results)),
        "reject_rate": float(np.mean([r.rejected for r in results])) if results else 0.0,
        "recover_rate": float(np.mean([r.recovered_split for r in results])) if results else 0.0,
    }
    return {"overall": overall, "by_vimp": by_vimp}


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def run_fico_experiment(
    *,
    url: str = HELOC_URL,
    cache_path: Optional[Path] = None,
    top_k_split: int = 3,
    vimps: Sequence[str] = VIMP_NOTIONS,
    n_perm: int = 200,
    n_estimators: int = 150,
    n_perm_cate: int = 50,
    max_n_per_arm: Optional[int] = 1500,
    level_feature: float = 0.99,
    level_across_feature: float = 0.95,
    seed: int = 2026,
    n_jobs: int = 2,
    n_bins: int = 4,
) -> Dict[str, object]:
    raw = load_heloc(url=url, cache_path=cache_path)
    clean = clean_heloc(raw)
    split_feats = top_lr_features(clean, top_k=top_k_split, seed=seed)
    results: List[CFPermResult] = []
    subset_rows: List[SubsetLocalRow] = []
    meta = {
        "n_raw": int(len(raw)),
        "n_clean": int(len(clean)),
        "split_features": split_feats,
        "n_perm": int(n_perm),
        "n_estimators": int(n_estimators),
        "n_perm_cate": int(n_perm_cate),
        "max_n_per_arm": max_n_per_arm,
        "vimps": list(vimps),
        "level_feature": level_feature,
        "level_across_feature": level_across_feature,
        "seed": seed,
    }
    for i, feat in enumerate(split_feats):
        X, Y, W, names = median_split_batches(
            clean, feat, max_n_per_arm=max_n_per_arm, seed=seed + i
        )
        for vimp in vimps:
            print(
                f"[fico] split={feat} vimp={vimp} n={len(Y)} n_perm={n_perm} B={n_estimators}",
                flush=True,
            )
            res = cfperm_test(
                X,
                Y,
                W,
                names,
                vimp=vimp,
                split_feature=feat,
                n_perm=n_perm,
                n_estimators=n_estimators,
                n_perm_cate=n_perm_cate,
                level_feature=level_feature,
                level_across_feature=level_across_feature,
                seed=seed + 100 * i,
                n_jobs=n_jobs,
            )
            results.append(res)
            # localize on rejected features, else top-3 VIMP
            loc_feats = res.rejected_features or res.top_features[:3]
            if feat not in loc_feats:
                loc_feats = [feat] + list(loc_feats)
            subset_rows.extend(
                subset_localize(
                    X,
                    Y,
                    W,
                    names,
                    loc_feats[:4],
                    split_feature=feat,
                    vimp=vimp,
                    n_bins=n_bins,
                )
            )
    power = power_summary(results)
    return {
        "meta": meta,
        "results": [asdict(r) for r in results],
        "subset_localization": [asdict(r) for r in subset_rows],
        "power": power,
    }


def write_outputs(payload: Dict[str, object], out_dir: Path) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "fico_cfperm_results.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    # flat tables
    pd.DataFrame(payload["results"]).to_csv(out_dir / "cfperm_cells.csv", index=False)
    pd.DataFrame(payload["subset_localization"]).to_csv(
        out_dir / "subset_localization.csv", index=False
    )
    power = payload["power"]
    rows = [{"scope": "overall", **power["overall"]}]
    for v, d in power["by_vimp"].items():
        rows.append({"scope": v, **d})
    pd.DataFrame(rows).to_csv(out_dir / "power_summary.csv", index=False)

    lines = [
        "# FICO HELOC CFPerm — VIMP power + subset localization",
        "",
        f"- Clean n = `{payload['meta']['n_clean']}` (raw `{payload['meta']['n_raw']}`).",
        f"- Split features (top-{len(payload['meta']['split_features'])} LR |coef|): "
        + ", ".join(f"`{c}`" for c in payload["meta"]["split_features"])
        + ".",
        f"- Settings: `B={payload['meta']['n_estimators']}` trees, "
        f"`n_perm={payload['meta']['n_perm']}`, "
        f"`n_perm_cate={payload['meta']['n_perm_cate']}`, "
        f"`max_n_per_arm={payload['meta']['max_n_per_arm']}`.",
        "",
        "## Power",
        "",
        "| scope | reject_rate | recover_rate | top1_hit | top3_hit | n_cells |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    ov = power["overall"]
    lines.append(
        f"| overall | {ov['reject_rate']:.3f} | {ov['recover_rate']:.3f} | — | — | {int(ov['n_cells'])} |"
    )
    for v, d in power["by_vimp"].items():
        lines.append(
            f"| {v} | {d['reject_rate']:.3f} | {d['recover_rate']:.3f} | "
            f"{d['top1_hit_rate']:.3f} | {d['top3_hit_rate']:.3f} | {int(d['n_cells'])} |"
        )
    lines += [
        "",
        "## Cells",
        "",
        "| split | vimp | reject | recovered | top features |",
        "|---|---|---:|---:|---|",
    ]
    for r in payload["results"]:
        lines.append(
            f"| `{r['split_feature']}` | {r['vimp']} | {r['rejected']} | "
            f"{r['recovered_split']} | {', '.join('`'+t+'`' for t in r['top_features'][:3])} |"
        )
    lines += [
        "",
        "## Subset localization",
        "",
        "Quantile bins of rejected / top VIMP features; `delta_y = mean(Y|W=1) - mean(Y|W=0)` "
        "localizes where the batch gap concentrates. See `subset_localization.csv`.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
