"""Reference PO-risk vs post-RFPerm re-fit PO-risk.

Reference (frozen):  PO_ref_i = |Y_i − f_ref(X_i)|
  f_ref = OnlineRFPerm burn-in model (never updated).

Probe (rolling):     PO_probe_i = |Y_i − probe(X_i)|
  probe re-fit each step on previous batch (stale vs OOD).

Refit (recent/OOD):  PO_refit from μ0 fit on recent control R, scored on O.
  Triggered when OnlineRFPerm rejects (re-train the PO-learner).
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from agod.po_iptw import instance_po_risk, po_iptw_weights
from agod.po_refit import build_recent_ood_windows, refit_po_on_windows


def reference_po_from_fref(f_ref, X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Frozen OnlineRFPerm reference residual PO on a batch."""
    pred = np.asarray(f_ref.predict(X), float)
    return np.abs(np.asarray(y, float).ravel() - pred)


def probe_po(
    probe,
    X: np.ndarray,
    y: np.ndarray,
    *,
    batch_po: float | None = None,
    mix: float = 0.0,
) -> np.ndarray:
    """Rolling-probe residual PO (optionally blended with batch PO lift)."""
    pred = np.asarray(probe.predict(X), float)
    return instance_po_risk(y, pred, batch_po=batch_po, mix=mix)


def refit_po_current(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_recent: int = 1,
    window_mode: str = "recent_ood",
    blend_mu_gap: float = 0.25,
) -> np.ndarray:
    """Re-train PO-learner on recent control; score PO on current OOD batch."""
    windows = build_recent_ood_windows(
        stream, t, n_recent=n_recent, window_mode=window_mode  # type: ignore[arg-type]
    )
    return refit_po_on_windows(windows, seed=seed, blend_mu_gap=blend_mu_gap)


def po_quality_vs_truth(po: np.ndarray, truth_hard: np.ndarray) -> Dict[str, float]:
    """How well PO ranks truly hard rows (|residual| of a holdout oracle)."""
    po = np.asarray(po, float).ravel()
    truth = np.asarray(truth_hard, float).ravel()
    if len(po) < 8:
        return {
            "spearman": float("nan"),
            "pearson": float("nan"),
            "topk_overlap": float("nan"),
            "po_mean": float(np.mean(po)) if len(po) else float("nan"),
        }
    r_po = po.argsort().argsort().astype(float)
    r_tr = truth.argsort().argsort().astype(float)
    spearman = float(np.corrcoef(r_po, r_tr)[0, 1])
    pearson = float(np.corrcoef(po, truth)[0, 1])
    k = max(1, len(po) // 5)
    top_po = set(np.argsort(po)[-k:])
    top_tr = set(np.argsort(truth)[-k:])
    overlap = float(len(top_po & top_tr) / k)
    return {
        "spearman": spearman,
        "pearson": pearson,
        "topk_overlap": overlap,
        "po_mean": float(np.mean(po)),
    }


def weights_from_po(po: np.ndarray, mode: str = "sqrt") -> np.ndarray:
    return po_iptw_weights(po, mode=mode)  # type: ignore[arg-type]
