
from __future__ import annotations

import numpy_compat  # noqa: F401 — NumPy≥1.25 removed np.bool8; shap+econml

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Literal

import numpy as np
import pandas as pd
import numpy as np
import pandas as pd
from typing import Dict, Optional, Sequence, Literal, Tuple, Any
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from econml.dml import CausalForestDML


# ============================================================
# Utility
# ============================================================
def _as_1d(x) -> np.ndarray:
    return np.asarray(x).reshape(-1)


def _safe_normalize(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    s = float(np.sum(v))
    return v / s if np.isfinite(s) and s > 0 else np.zeros_like(v)


# ============================================================
# Method A: CausalForestDML feature_importances_
# ============================================================
def _fit_cf_importance(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 300,
    min_samples_leaf: int = 10,
    max_depth: Optional[int] = None,
    normalize: bool = True,
) -> np.ndarray:
    """
    Fit a causal forest for batch-effect heterogeneity and return feature_importances_.

    W is the batch indicator (0=existing, 1=new). Under H0 (no shift),
    tau(x)=E[Y|X=x,W=1]-E[Y|X=x,W=0] = 0 for all x.

    Requires: econml + scikit-learn.
    """
    try:
        from econml.dml import CausalForestDML
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.linear_model import LogisticRegression
    except Exception as e:
        raise ImportError(
            "This method requires `econml` and `scikit-learn`.\n"
            "Install: pip install econml scikit-learn"
        ) from e
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W)
    p = X.shape[1]
    uniq = np.unique(W)
    binary01 = (
        uniq.size <= 2
        and np.all(np.isin(uniq, (0, 1)))
        and np.issubdtype(W.dtype, np.integer)
    )
    # Residualized W is float and not {0,1}; use continuous-treatment path + KFold
    # so econml never asks for propensity folds missing a treatment class.
    discrete_treatment = bool(binary01)
    if discrete_treatment:
        W_fit = W.astype(int)
        model_t = LogisticRegression(
            solver="lbfgs",
            max_iter=2000,
            random_state=seed,
        )
        cv_t = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    else:
        W_fit = W.astype(float)
        model_t = RandomForestRegressor(
            n_estimators=120,
            min_samples_leaf=5,
            max_features=max(1, int(round(p / 3))),
            random_state=seed,
            n_jobs=-1,
        )
        cv_t = KFold(n_splits=3, shuffle=True, random_state=seed)
    model_y = RandomForestRegressor(
        n_estimators=200,
        min_samples_leaf=5,
        max_features=max(1, int(round(p / 3))),
        random_state=seed,
        n_jobs=-1,
    )
    cf = CausalForestDML(
        model_y=model_y,
        model_t=model_t,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        max_depth=max_depth,
        discrete_treatment=discrete_treatment,
        cv=cv_t,
        random_state=seed,
    )
    cf.fit(Y=Y, T=W_fit, X=X)
    if not hasattr(cf, "feature_importances_"):
        raise AttributeError(
            "Your econml version does not expose `feature_importances_` on CausalForestDML."
        )
    imp = np.asarray(cf.feature_importances_, dtype=float)
    imp = imp[:p] if imp.shape[0] >= p else np.pad(imp, (0, p - imp.shape[0]))
    return _safe_normalize(imp) if normalize else imp


# ============================================================
# Method B: LOCO R-risk importance (R-learner)
# ============================================================
def _crossfit_nuisances(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_splits: int = 3,
    binary_outcome: bool = False,
    clip_e: float = 0.01,
):
    """
    Cross-fit:
      m_hat(x) = E[Y|X=x]
      e_hat(x) = P(W=1|X=x)

    Returns cross-fitted (m_hat, e_hat) evaluated at all rows.
    """
    try:
        from sklearn.model_selection import StratifiedKFold
        from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
    except Exception as e:
        raise ImportError("This method requires `scikit-learn`. Install: pip install scikit-learn") from e
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    if binary_outcome:
        model_m = RandomForestClassifier(
            n_estimators=200, min_samples_leaf=5,
            max_features=max(1, int(round(np.sqrt(p)))),
            random_state=seed, n_jobs=-1
        )
    else:
        model_m = RandomForestRegressor(
            n_estimators=200, min_samples_leaf=5,
            max_features=max(1, int(round(p / 3))),
            random_state=seed, n_jobs=-1
        )
    model_e = LogisticRegression(solver="lbfgs", max_iter=2000, random_state=seed)
    m_hat = np.zeros(n, dtype=float)
    e_hat = np.zeros(n, dtype=float)
    for tr, te in cv.split(X, W):
        Xm_tr, Xm_te = X[tr], X[te]
        Y_tr, W_tr = Y[tr], W[tr]
        # Fit m
        model_m_fold = model_m.__class__(**model_m.get_params())
        model_m_fold.fit(Xm_tr, Y_tr)
        if hasattr(model_m_fold, "predict_proba"):
            m_hat[te] = model_m_fold.predict_proba(Xm_te)[:, 1]
        else:
            m_hat[te] = model_m_fold.predict(Xm_te)
        # Fit e
        model_e_fold = LogisticRegression(**model_e.get_params())
        model_e_fold.fit(Xm_tr, W_tr)
        e_hat[te] = model_e_fold.predict_proba(Xm_te)[:, 1]
    e_hat = np.clip(e_hat, clip_e, 1.0 - clip_e)
    return m_hat, e_hat


def _fit_tau_rlearner_weighted(
    X: np.ndarray,
    Y_tilde: np.ndarray,
    W_tilde: np.ndarray,
    *,
    seed: int = 0,
    clip_wtilde: float = 1e-3,
):
    """
    Fit tau(x) by weighted regression in the R-learner objective:
      minimize sum_i (Y_tilde_i - W_tilde_i * tau(X_i))^2

    Using z_i = Y_tilde_i / W_tilde_i and weights w_i = W_tilde_i^2
    on rows with |W_tilde_i| > clip_wtilde.
    """
    try:
        from sklearn.base import clone
        from sklearn.linear_model import Ridge
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
    except Exception as e:
        raise ImportError("This method requires `scikit-learn`. Install: pip install scikit-learn") from e
    X = np.asarray(X)
    Y_tilde = _as_1d(Y_tilde)
    W_tilde = _as_1d(W_tilde)
    mask = np.abs(W_tilde) > clip_wtilde
    if np.sum(mask) < max(30, X.shape[1] + 5):
        # Not enough stable rows; return a simple constant-zero model
        class _Zero:
            def predict(self, X_):
                return np.zeros(X_.shape[0], dtype=float)
        return _Zero()

    z = Y_tilde[mask] / W_tilde[mask]
    weights = (W_tilde[mask] ** 2)
    # Default tau-model: standardized Ridge (fast, stable, supports sample_weight)
    tau_model = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("ridge", Ridge(alpha=1.0, random_state=seed)),
        ]
    )
    tau_model.fit(X[mask], z, ridge__sample_weight=weights)
    return tau_model


def r_risk(
    Y_tilde: np.ndarray,
    W_tilde: np.ndarray,
    tau_hat: np.ndarray,
) -> float:
    """
    Empirical R-risk:
      R = mean( (Y_tilde - W_tilde * tau_hat)^2 )
    """
    Y_tilde = _as_1d(Y_tilde)
    W_tilde = _as_1d(W_tilde)
    tau_hat = _as_1d(tau_hat)
    return float(np.mean((Y_tilde - W_tilde * tau_hat) ** 2))


def importance_loco_r_risk(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_splits: int = 3,
    binary_outcome: bool = False,
    clip_e: float = 0.01,
    clip_wtilde: float = 1e-3,
    normalize: bool = False,
) -> np.ndarray:
    """
    LOCO importance for a causal learner via R-risk (R-learner).
    Steps:
      1) Cross-fit nuisances m_hat(x), e_hat(x)
      2) Form residuals: Y_tilde = Y - m_hat(X), W_tilde = W - e_hat(X)
      3) Fit tau_full(x) by weighted R-learner regression
      4) For each feature j, refit tau on X without feature j, compute risk_j
      5) importance_j = risk_j - risk_full
    """
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    m_hat, e_hat = _crossfit_nuisances(
        X, Y, W, seed=seed, n_splits=n_splits, binary_outcome=binary_outcome, clip_e=clip_e
    )
    Y_tilde = Y - m_hat
    W_tilde = W - e_hat

    tau_full_model = _fit_tau_rlearner_weighted(
        X, Y_tilde, W_tilde, seed=seed, clip_wtilde=clip_wtilde
    )
    tau_full = tau_full_model.predict(X)
    risk_full = r_risk(Y_tilde, W_tilde, tau_full)
    p = X.shape[1]
    imps = np.zeros(p, dtype=float)
    for j in range(p):
        X_minus = np.delete(X, j, axis=1)
        tau_j_model = _fit_tau_rlearner_weighted(
            X_minus, Y_tilde, W_tilde, seed=seed + 17 * (j + 1), clip_wtilde=clip_wtilde
        )
        tau_j = tau_j_model.predict(X_minus)
        risk_j = r_risk(Y_tilde, W_tilde, tau_j)
        imps[j] = risk_j - risk_full
    return _safe_normalize(imps) if normalize else imps


# ============================================================
# CFPerm wrapper: choose importance and run permutation p-values
# ============================================================
VimpMethod = Literal["cf_feature_importance", "loco_r_risk"]


@dataclass
class CFPermResult:
    feature_names: list
    imp_observed: pd.Series
    p_values: pd.Series
    q1_upper: pd.Series
    cross_threshold_q2: float
    rejected_features: list
    perm_importances: np.ndarray  # shape (p, n_perm)
    method: str


def cfperm_feature_pvals(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    feature_names: Optional[list] = None,
    method: VimpMethod = "cf_feature_importance",
    n_perm: int = 200,
    q1: float = 0.99,  # feature-wise upper quantile
    q2: float = 0.95,  # cross-feature quantile over q1_upper vector
    seed: int = 0,
    verbose: bool = False,
    # Method A params
    n_estimators: int = 300,
    min_samples_leaf: int = 10,
    max_depth: Optional[int] = None,
    normalize_importance: bool = True,
    # Method B params
    n_splits: int = 3,
    binary_outcome: bool = False,
    clip_e: float = 0.01,
    clip_wtilde: float = 1e-3,
) -> CFPermResult:
    rng = np.random.default_rng(seed)
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    if feature_names is None:
        feature_names = [f"X{j+1}" for j in range(p)]
    feature_names = list(feature_names)
    def compute_importance(X_, Y_, W_, *, seed_):
        if method == "cf_feature_importance":
            return _fit_cf_importance(
                X_, Y_, W_,
                seed=seed_,
                n_estimators=n_estimators,
                min_samples_leaf=min_samples_leaf,
                max_depth=max_depth,
                normalize=normalize_importance,
            )
        elif method == "loco_r_risk":
            # For LOCO R-risk, normalization is usually not required; keep as-is
            return importance_loco_r_risk(
                X_, Y_, W_,
                seed=seed_,
                n_splits=n_splits,
                binary_outcome=binary_outcome,
                clip_e=clip_e,
                clip_wtilde=clip_wtilde,
                normalize=False,
            )
        else:
            raise ValueError(f"Unknown method: {method}")
    # Observed importances
    imp_obs = compute_importance(X, Y, W, seed_=seed)
    # Permutation importances (permute batch labels W and refit)
    perm_imps = np.zeros((p, n_perm), dtype=float)
    for b in range(n_perm):
        Wb = rng.permutation(W)
        perm_imps[:, b] = compute_importance(X, Y, Wb, seed_=seed + 10_000 + b)
        if verbose:
            print(b)
    # One-sided feature-wise p-values: P(Imp_perm >= Imp_obs)
    pvals = np.sum(perm_imps >= imp_obs[:, None], axis=1) / n_perm
    # Feature-wise upper quantile q1
    q1_upper = np.quantile(perm_imps, q=q1, axis=1)
    # Cross-feature threshold q2 over the vector of q1 uppers
    cross_thr = float(np.quantile(q1_upper, q=q2))
    rejected = [feature_names[j] for j in range(p) if imp_obs[j] > cross_thr]
    return CFPermResult(
        feature_names=feature_names,
        imp_observed=pd.Series(imp_obs, index=feature_names, name="importance_observed"),
        p_values=pd.Series(pvals, index=feature_names, name="p_value"),
        q1_upper=pd.Series(q1_upper, index=feature_names, name=f"perm_q{q1}_upper"),
        cross_threshold_q2=cross_thr,
        rejected_features=rejected,
        perm_importances=perm_imps,
        method=method,
    )


def LMGeneration(n = 1000, p = 20, beta = [1,1,1,1,1,1,1], variances = 1, cor = 0.3, eps_noi = 1, mean = 0):
    '''
    Generate a linear DGP(Data-Generating Process) with
    Y = X\beta + \epsilon with AR(1) correation structure with Gaussian Noise.
    '''
    n_signal = len(beta)
    mean_vector = [mean] * n_signal
    cov_mat = np.zeros((n_signal, n_signal))
    for i in range(n_signal):
        for j in range(n_signal):
            cov_mat[i, j] = cor ** abs(i - j) * variances
    cov = np.dot(cov_mat, cov_mat.T)
    X_design = np.random.multivariate_normal(mean_vector, cov, size=n)
    #generate the response:
    Y = np.dot(X_design, beta) + np.random.normal(loc = 0, scale = eps_noi, size = n)
    X_noise = np.zeros((n, p - n_signal))
    for k in range(p - n_signal):
        X_noise[:, k] = np.random.normal(loc = 0, scale = 1, size = n)
    X_design = np.concatenate((X_design, X_noise), axis = 1)
    df_design = pd.DataFrame(np.concatenate((X_design, Y.reshape(-1, 1)), axis = 1))
    df_design.columns = ["X" + str(i) for i in range(p)] + ["Y"]
    return df_design




    
