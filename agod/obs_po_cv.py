"""CV-MSE selection of observation-level PO weight hyperparameters.

On a rejected OOD batch, split rows into K folds; for each candidate
(power, temper, family) fit a cheap regressor with the corresponding
sample weights and pick the setting with lowest mean hold-out MSE.

Power=0 / temper=0 → uniform. Lets calm packs fall back without a
hand-tuned drift gate, while beijing-like batches can pick a soft power.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

from agod.obs_po_weights import obs_po_to_weights

DEFAULT_POWERS: Tuple[float, ...] = (0.0, 0.125, 0.25, 1.0 / 3.0, 0.5)
DEFAULT_TEMPERS: Tuple[float, ...] = (0.0, 0.35, 0.55, 0.75)
DEFAULT_FAMILIES: Tuple[str, ...] = (
    "uniform",
    "qrt",
    "log1p",
    "softmax",
    "hard_support",
    "quantile",
)


def _cheap_rf(seed: int) -> RandomForestRegressor:
    return RandomForestRegressor(
        n_estimators=20,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )


def cv_mse_for_weights(
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    *,
    n_folds: int = 3,
    seed: int = 0,
    model_factory: Optional[Callable[[int], object]] = None,
) -> float:
    """K-fold MSE of a weighted fit (lower is better)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    w = np.asarray(w, float).ravel()
    n = len(y)
    if n < max(2 * n_folds, 8):
        # too small: in-sample residual proxy
        m = (model_factory or _cheap_rf)(seed)
        m.fit(X, y, sample_weight=w)
        pred = m.predict(X)
        return float(np.mean((y - pred) ** 2))

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    losses: List[float] = []
    for fold_i, (tr, te) in enumerate(kf.split(X)):
        m = (model_factory or _cheap_rf)(seed + fold_i)
        m.fit(X[tr], y[tr], sample_weight=w[tr])
        pred = m.predict(X[te])
        losses.append(float(np.mean((y[te] - pred) ** 2)))
    return float(np.mean(losses))


def make_power_weights(
    po: np.ndarray,
    power: float,
    *,
    temper: float = 1.0,
    clip: tuple[float, float] = (0.05, 20.0),
    eps: float = 1e-6,
) -> np.ndarray:
    """w ∝ PO^power, then temper-mix with uniform. power≤0 → uniform."""
    if float(power) <= 0.0 or float(temper) <= 0.0:
        return np.ones(len(np.asarray(po).ravel()), float)
    return obs_po_to_weights(
        po, mode="qrt", power=float(power), temper=float(temper), clip=clip, eps=eps
    )


def cv_select_power(
    X: np.ndarray,
    y: np.ndarray,
    po: np.ndarray,
    *,
    powers: Sequence[float] = DEFAULT_POWERS,
    tempers: Sequence[float] = DEFAULT_TEMPERS,
    n_folds: int = 3,
    seed: int = 0,
    temper_cap: float | None = None,
    model_factory: Optional[Callable[[int], object]] = None,
) -> Dict[str, object]:
    """Grid-search (power, temper) by CV-MSE on the OOD batch.

    ``temper_cap`` (e.g. adaptive_temper(drift)) upper-bounds λ so CV
    cannot overshoot the drift budget; None → full temper grid.
    """
    rows: List[dict] = []
    best = None
    for p in powers:
        for lam in tempers:
            if temper_cap is not None:
                lam = min(float(lam), float(temper_cap))
            w = make_power_weights(po, float(p), temper=float(lam))
            mse = cv_mse_for_weights(
                X, y, w, n_folds=n_folds, seed=seed, model_factory=model_factory
            )
            row = {"power": float(p), "temper": float(lam), "cv_mse": mse}
            rows.append(row)
            if best is None or mse < best["cv_mse"]:
                best = row

    assert best is not None
    w_best = make_power_weights(po, best["power"], temper=best["temper"])
    return {
        "power": best["power"],
        "temper": best["temper"],
        "cv_mse": best["cv_mse"],
        "weights": w_best,
        "grid": rows,
    }


def cv_select_family(
    X: np.ndarray,
    y: np.ndarray,
    po: np.ndarray,
    *,
    families: Sequence[str] = DEFAULT_FAMILIES,
    temper: float = 0.55,
    n_folds: int = 3,
    seed: int = 0,
    topk_frac: float = 0.2,
    model_factory: Optional[Callable[[int], object]] = None,
) -> Dict[str, object]:
    """Pick a discrete weight family by CV-MSE (temper shared for soft maps)."""
    rows: List[dict] = []
    best = None
    lam = float(temper)
    for fam in families:
        if fam == "uniform" or lam <= 0.0:
            w = np.ones(len(np.asarray(po).ravel()), float)
            fam_eff = "uniform"
        else:
            fam_eff = fam
            w = obs_po_to_weights(
                po,
                mode=fam,  # type: ignore[arg-type]
                temper=lam,
                topk_frac=topk_frac,
                boost_max=3.0,
            )
        mse = cv_mse_for_weights(
            X, y, w, n_folds=n_folds, seed=seed, model_factory=model_factory
        )
        row = {"family": fam_eff, "temper": lam if fam_eff != "uniform" else 0.0, "cv_mse": mse}
        rows.append(row)
        if best is None or mse < best["cv_mse"]:
            best = {**row, "weights": w}

    assert best is not None
    return {
        "family": best["family"],
        "temper": best["temper"],
        "cv_mse": best["cv_mse"],
        "weights": best["weights"],
        "grid": rows,
    }
