"""Sliding-window candidate feature bank for monitoring.

If compute allows, keep rolling statistics on a window and score them
against the same frozen D_ref as PO-risk / MMD. Online PCA is cheap.

These are candidate monitoring features (and later clever covariates),
not a unique decomposition of the drop.

MMD口径 stays MMD²(X_window, X_ref), never pairwise vs history.
No online-bootstrap.
"""
from __future__ import annotations

from collections import deque
from typing import Mapping

import numpy as np
from sklearn.decomposition import PCA, IncrementalPCA

from logo_modality import as_groups, brier_or_mse, fit_serving, serving_mu
from online_rfperm import FrozenRFPerm
from streaming_po_risk import mmd_vs_reference, rbf_bandwidth, streaming_po_and_mse


def _as_2d(X) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


def subspace_gap(P_ref: np.ndarray, P_now: np.ndarray) -> float:
    """1 − mean canonical overlap of two orthonormal bases. 0 = same subspace."""
    A = np.asarray(P_ref, dtype=float)
    B = np.asarray(P_now, dtype=float)
    k = min(A.shape[1], B.shape[1])
    if k == 0:
        return 0.0
    M = A[:, :k].T @ B[:, :k]
    s = np.linalg.svd(M, compute_uv=False)
    return float(1.0 - np.mean(np.clip(s, 0.0, 1.0)))


def recon_error(pca: PCA, X) -> float:
    X = _as_2d(X)
    Z = pca.transform(X)
    Xh = pca.inverse_transform(Z)
    return float(np.mean((X - Xh) ** 2))


class SlidingWindowBank:
    """Frozen D_ref probes + a sliding window of the last `window` rows.

    Cheap features run every step. PO-risk is optional (heavier).
    Per-group PCA is the cheap stand-in for LOGO localization on P(X).
    """

    def __init__(
        self,
        X_ref,
        Y_ref,
        *,
        window: int = 200,
        n_components: int = 4,
        groups: Mapping | None = None,
        seed: int = 2026,
        with_po: bool = False,
    ):
        self.X_ref = _as_2d(X_ref)
        self.Y_ref = np.asarray(Y_ref, dtype=float).ravel()
        self.window = max(int(window), 8)
        self.n_components = int(n_components)
        self.seed = int(seed)
        self.with_po = bool(with_po)
        self.groups = as_groups(groups) if groups else None

        k = min(self.n_components, self.X_ref.shape[1], max(self.X_ref.shape[0] - 1, 1))
        self.pca_ref = PCA(n_components=k, random_state=seed)
        self.pca_ref.fit(self.X_ref)
        self.recon_ref = recon_error(self.pca_ref, self.X_ref)
        self.mean_ref = self.X_ref.mean(axis=0)
        self.pca_live = IncrementalPCA(n_components=k)
        self.pca_live.partial_fit(self.X_ref[: max(k + 8, min(len(self.X_ref), 256))])

        self.sigma = rbf_bandwidth(self.X_ref, seed=seed)
        self.probe = FrozenRFPerm(self.X_ref, self.Y_ref, seed=seed)
        self.serve, self.binary = fit_serving(self.X_ref, self.Y_ref, seed=seed)
        self.brier_ref = brier_or_mse(
            self.Y_ref, serving_mu(self.serve, self.X_ref, self.binary), self.binary
        )

        self.group_pca = {}
        if self.groups:
            for name, idx in self.groups.items():
                Xg = self.X_ref[:, idx]
                kg = min(2, max(1, Xg.shape[1] - 1))
                pca = PCA(n_components=kg, random_state=seed)
                pca.fit(Xg)
                self.group_pca[name] = {
                    "pca": pca,
                    "recon_ref": recon_error(pca, Xg),
                    "idx": idx,
                }

        self.buf_x: deque[np.ndarray] = deque()
        self.buf_y: deque[np.ndarray] = deque()
        self.n_buf = 0
        self.last: dict = {}

    def _trim(self) -> None:
        while self.n_buf > self.window and self.buf_x:
            x0 = self.buf_x.popleft()
            self.buf_y.popleft()
            self.n_buf -= len(x0)

    def _window_xy(self):
        Xw = np.vstack(self.buf_x)
        Yw = np.concatenate(self.buf_y)
        return Xw, Yw

    def step(self, X_new, Y_new) -> dict:
        X_new = _as_2d(X_new)
        Y_new = np.asarray(Y_new, dtype=float).ravel()
        self.buf_x.append(X_new)
        self.buf_y.append(Y_new)
        self.n_buf += len(X_new)
        self._trim()
        Xw, Yw = self._window_xy()

        k_fit = max(self.pca_live.n_components, 8)
        self.pca_live.partial_fit(X_new[: max(len(X_new), k_fit)])

        recon = recon_error(self.pca_ref, Xw)
        score_now = self.pca_ref.transform(Xw).mean(axis=0)
        score_ref = self.pca_ref.transform(self.X_ref).mean(axis=0)
        feat = {
            "n_window": int(len(Xw)),
            "mean_l2": float(np.linalg.norm(Xw.mean(axis=0) - self.mean_ref)),
            "pca_recon": float(recon),
            "pca_recon_excess": float(recon - self.recon_ref),
            "pca_score_l2": float(np.linalg.norm(score_now - score_ref)),
            "pca_subspace_gap": subspace_gap(self.pca_ref.components_.T, self.pca_live.components_.T),
            "mmd_vs_ref": mmd_vs_reference(self.X_ref, Xw, sigma=self.sigma, seed=self.seed),
            "brier": brier_or_mse(Yw, serving_mu(self.serve, Xw, self.binary), self.binary),
            "brier_excess": None,
            "rfperm_T": None,
            "po_risk": None,
        }
        feat["brier_excess"] = float(feat["brier"] - self.brier_ref)
        rf = self.probe.step(Xw, Yw)
        feat["rfperm_T"] = float(rf["rfperm_T"])
        feat["rfperm_hop"] = bool(rf["rfperm_hop"])
        if self.with_po:
            po, _ = streaming_po_and_mse(self.X_ref, self.Y_ref, Xw, Yw, seed=self.seed)
            feat["po_risk"] = float(po)
        if self.group_pca:
            recon_g = {}
            for name, g in self.group_pca.items():
                Xg = Xw[:, g["idx"]]
                recon_g[name] = recon_error(g["pca"], Xg) - float(g["recon_ref"])
            feat["pca_recon_group"] = recon_g
            feat["pca_recon_group_argmax"] = max(recon_g, key=recon_g.get)
        self.last = feat
        return feat

    def vector(self, feat: dict | None = None) -> np.ndarray:
        """Flat candidate vector for a shadow / monitoring model."""
        f = self.last if feat is None else feat
        keys = (
            "mean_l2",
            "pca_recon_excess",
            "pca_score_l2",
            "pca_subspace_gap",
            "mmd_vs_ref",
            "brier_excess",
            "rfperm_T",
        )
        vals = [float(f[k]) for k in keys]
        if f.get("pca_recon_group"):
            for name in sorted(f["pca_recon_group"]):
                vals.append(float(f["pca_recon_group"][name]))
        if f.get("po_risk") is not None:
            vals.append(float(f["po_risk"]))
        return np.asarray(vals, dtype=float)
