"""Observation matrix X for boards / stats-ToT: edge features ∥ stream features.

Contract (aligned with review)
------------------------------
- X = ``hstack(edge, stream)`` — same row index / time order
- Y = downstream target (click / convert / traffic / …)
- ToT *intermediate* = policy Thought state (adapt / α / freeze + scorecard),
  **not** a predicted Ŷ

Edge block: UI / graph-edge tabular feats (e.g. TencentGR ``e_*``/``u_*``/``i_*``).
Stream block: time-ordered pack feats (metro / stocks / pm25 / waymo_proxy / …).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class FeatBlocks:
    """Named horizontal blocks that make up observation X."""

    X: Array
    y: Optional[Array]
    names: List[str]
    edge_dim: int
    stream_dim: int
    n: int

    @property
    def edge_slice(self) -> slice:
        return slice(0, self.edge_dim)

    @property
    def stream_slice(self) -> slice:
        return slice(self.edge_dim, self.edge_dim + self.stream_dim)

    def to_meta(self) -> Dict[str, Any]:
        return {
            "n": int(self.n),
            "d": int(self.X.shape[1]) if self.X.ndim == 2 else 0,
            "edge_dim": int(self.edge_dim),
            "stream_dim": int(self.stream_dim),
            "layout": "edge||stream",
            "y_is_target": True,
            "intermediate_is": "policy_thought_not_yhat",
        }


def _as_2d(a: Array, *, name: str) -> Array:
    x = np.asarray(a, dtype=np.float64)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    if x.ndim != 2:
        raise ValueError(f"{name} must be 1d or 2d, got shape {x.shape}")
    return x


def _default_names(prefix: str, d: int, given: Optional[Sequence[str]]) -> List[str]:
    if given is None:
        return [f"{prefix}{j}" for j in range(d)]
    names = [str(s) for s in given]
    if len(names) != d:
        raise ValueError(
            f"{prefix} names length {len(names)} != feature dim {d}"
        )
    return names


def align_rows(
    X_edge: Array,
    X_stream: Array,
    *,
    y: Optional[Array] = None,
    n: Optional[int] = None,
) -> Tuple[Array, Array, Optional[Array]]:
    """Truncate to a common leading row count (time / index order).

    Boards and ToT assume the same observation index across blocks; we do
    not invent join keys here — callers pre-align by ``e_last_ts`` / pack time.
    """
    xe = _as_2d(X_edge, name="X_edge")
    xs = _as_2d(X_stream, name="X_stream")
    n_common = min(xe.shape[0], xs.shape[0])
    if y is not None:
        yy = np.asarray(y).reshape(-1)
        n_common = min(n_common, yy.shape[0])
    else:
        yy = None
    if n is not None:
        n_common = min(n_common, int(n))
    if n_common <= 0:
        raise ValueError("no overlapping rows to concat edge||stream")
    xe, xs = xe[:n_common], xs[:n_common]
    if yy is not None:
        yy = yy[:n_common]
    return xe, xs, yy


def concat_edge_stream(
    X_edge: Array,
    X_stream: Array,
    *,
    y: Optional[Array] = None,
    edge_names: Optional[Sequence[str]] = None,
    stream_names: Optional[Sequence[str]] = None,
    n: Optional[int] = None,
    standardize: bool = False,
) -> FeatBlocks:
    """Build ``X = [edge | stream]`` with prefixed feature names.

    Parameters
    ----------
    X_edge, X_stream
        Feature matrices (n×d). Row-aligned in time / board index order.
    y
        Optional target; carried through truncated to the same ``n``.
    standardize
        If True, z-score each column after concat (W1-style board hygiene).
    """
    xe, xs, yy = align_rows(X_edge, X_stream, y=y, n=n)
    de, ds = int(xe.shape[1]), int(xs.shape[1])
    en = _default_names("edge:", de, edge_names)
    sn = _default_names("stream:", ds, stream_names)
    # Avoid double-prefix if caller already tagged
    en = [c if c.startswith("edge:") else f"edge:{c}" for c in en]
    sn = [c if c.startswith("stream:") else f"stream:{c}" for c in sn]
    X = np.hstack([xe, xs]).astype(np.float64, copy=False)
    if standardize and X.size:
        mu = X.mean(axis=0)
        sd = X.std(axis=0)
        X = (X - mu) / (sd + 1e-8)
    return FeatBlocks(
        X=X,
        y=yy,
        names=en + sn,
        edge_dim=de,
        stream_dim=ds,
        n=int(X.shape[0]),
    )


def split_edge_stream(
    blocks: FeatBlocks,
) -> Tuple[Array, Array]:
    """Inverse view: recover edge / stream matrices from a FeatBlocks."""
    return blocks.X[:, blocks.edge_slice], blocks.X[:, blocks.stream_slice]


def tot_observation_spec() -> Dict[str, str]:
    """Human-readable X/Y/intermediate contract for stats-ToT."""
    return {
        "X": "edge_features || stream_features (hstack, shared row index)",
        "Y": "downstream target on that row (not a ToT intermediate)",
        "intermediate": (
            "policy Thought = {adapt, alpha, freeze, rank_eff, mse_eff, "
            "expected_flops, burn} — not Ŷ"
        ),
        "V": "(mse_eff tier, mse_eff, rank_eff, -flops)",
    }
