"""Adaptation opportunity for novel-category cold-start.

Core claim
----------
Attribution + separation are *batch-local* actuators. The product question is
not "did we separate this batch?" but

    O_t = opportunity to adapt to new categories at time t
        = how much expected novel-acc lift is available *now*,
          given that P_{t+1} is unknown.

Decomposition (all factors in [0, 1])
------------------------------------
    O_local  = D · L · R · M     # act on this batch
    O_persist = O_local · Stab   # reuse S / offset next batch (usually <<)

  D  detectability   — LIVE/REF error rise (RFPerm gate)
  L  localizability  — VIMP mass on attributed subspace S
  R  separability    — membership gap / room to push off known prototypes
  M  actable mass    — enough high-mem / private-live rows to spend
  Stab footprint stability vs previous batch (unknown next-batch discount)

Policy
------
  O_local  ≥ τ_sep     → run representation separation on S
  O_local  ≥ τ_exp and M high → allow head expand
  O_persist < τ_keep   → decay / drop persisted S-offset (do not bet on t+1)

Benchmarks for this track (external)
------------------------------------
  C-GCD: CIFAR-100, ImageNet-100, Tiny-ImageNet, CUB
         (MetaGCD / Happy: per-session Old / New / All)
  GCD/NCD: CIFAR-10/100, ImageNet-100, CUB, Cars, Herbarium
  PDA/OSDA/UniDA: Office-Home, VisDA, DomainNet (H-score)
  Our proxies: synthetic C-GCD + Office-Home PDA + tabular boards
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np

from weighted_cold_start import _membership_conc


def _clip01(x: float) -> float:
    return float(np.clip(x, 0.0, 1.0))


def detectability(error_ratio: float, *, thresh: float = 1.15) -> float:
    """D: how clearly LIVE is worse than REF. Soft step around the gate."""
    r = float(error_ratio)
    # 1.0 → 0; thresh → ~0.5; 1.5·thresh → ~1
    return _clip01((r - 1.0) / max(float(thresh) - 1.0 + 1e-6, 1e-6) * 0.5)


def localizability(vimp: np.ndarray, feat_idx: np.ndarray) -> float:
    """L: mass of attributed subspace S (sharp footprint ⇒ high)."""
    v = np.asarray(vimp, float)
    feat = np.asarray(feat_idx, int)
    if v.size == 0 or feat.size == 0:
        return 0.0
    return _clip01(float(v[feat].sum()))


def separability(
    membership: np.ndarray,
    *,
    X_live: Optional[np.ndarray] = None,
    prototypes_S: Optional[np.ndarray] = None,
    feat_idx: Optional[np.ndarray] = None,
    mem_thr: float = 0.55,
    margin: float = 1.25,
) -> Dict[str, float]:
    """R: can we still push high-mem rows off known centers on S?

    Uses membership concentration as a cheap proxy; if prototypes are given,
    also scores residual margin violations.
    """
    mem = np.asarray(membership, float).ravel()
    conc = _membership_conc(mem)
    high = mem >= float(mem_thr)
    if not high.any():
        high = mem >= float(np.quantile(mem, 0.75))
    # membership gap proxy: top vs bottom half
    mid = int(len(mem) // 2) or 1
    top = float(np.sort(mem)[-mid:].mean())
    bot = float(np.sort(mem)[:mid].mean())
    gap_mem = _clip01((top - bot) / 0.6)

    resid = 0.5
    if (
        X_live is not None and prototypes_S is not None and feat_idx is not None
        and high.any() and np.asarray(prototypes_S).size
    ):
        Xs = np.asarray(X_live, float)[:, np.asarray(feat_idx, int)]
        proto = np.asarray(prototypes_S, float)
        d = []
        for i in np.flatnonzero(high):
            d.append(float(np.linalg.norm(Xs[i] - proto, axis=1).min()))
        d = np.asarray(d, float)
        # fraction still inside / near margin → room to separate
        frac_in = float(np.mean(d < float(margin)))
        # also reward large existing gap (already separable geometry)
        gap_geo = _clip01(float(np.median(d)) / max(2.0 * float(margin), 1e-6))
        resid = _clip01(0.55 * frac_in + 0.45 * gap_geo)

    R = _clip01(0.40 * conc + 0.35 * gap_mem + 0.25 * resid)
    return {
        "R": R,
        "mem_conc": float(conc),
        "gap_mem": float(gap_mem),
        "resid": float(resid),
    }


def actable_mass(
    membership: np.ndarray,
    *,
    n_private: int = 0,
    expand_min: int = 16,
    mem_thr: float = 0.55,
) -> Dict[str, float]:
    """M: enough rows to spend on sep / expand without instant overfit."""
    mem = np.asarray(membership, float).ravel()
    n_high = int((mem >= float(mem_thr)).sum())
    n_eff = max(n_high, int(n_private))
    # sep can fire with few rows; expand needs closer to N*
    M_sep = _clip01(n_eff / max(4.0, 0.35 * float(expand_min)))
    M_exp = _clip01(n_eff / max(float(expand_min), 1.0))
    M = _clip01(0.55 * M_sep + 0.45 * M_exp)
    return {
        "M": M,
        "M_sep": float(M_sep),
        "M_exp": float(M_exp),
        "n_high": n_high,
        "n_private": int(n_private),
        "n_eff": n_eff,
    }


def footprint_stability_score(
    vimp_top_mass: float,
    prev_vimp_top_mass: Optional[float],
    membership: np.ndarray,
    *,
    prev_feat: Optional[np.ndarray] = None,
    feat_idx: Optional[np.ndarray] = None,
    mass_tol: float = 0.10,
) -> float:
    """Stab ∈ [0,1]: how safe is reusing last batch's S / offset?

    Unknown next-batch ⇒ default low unless mass + feature support overlap.
    """
    if prev_vimp_top_mass is None:
        return 0.35  # first batch: no persist claim
    mass_ok = abs(float(vimp_top_mass) - float(prev_vimp_top_mass)) <= float(mass_tol)
    conc = _membership_conc(membership)
    overlap = 0.0
    if prev_feat is not None and feat_idx is not None:
        a = set(int(x) for x in np.asarray(prev_feat).ravel())
        b = set(int(x) for x in np.asarray(feat_idx).ravel())
        overlap = len(a & b) / max(len(a | b), 1)
    raw = (0.40 * float(mass_ok) + 0.35 * overlap + 0.25 * conc)
    return _clip01(raw)


def adaptation_opportunity(
    *,
    error_ratio: float,
    vimp: np.ndarray,
    feat_idx: np.ndarray,
    membership: np.ndarray,
    n_private: int = 0,
    expand_min: int = 16,
    detect_thresh: float = 1.15,
    mem_thr: float = 0.55,
    margin: float = 1.25,
    X_live: Optional[np.ndarray] = None,
    prototypes_S: Optional[np.ndarray] = None,
    prev_vimp_top_mass: Optional[float] = None,
    prev_feat: Optional[np.ndarray] = None,
    tau_sep: float = 0.18,
    tau_exp: float = 0.28,
    tau_keep: float = 0.45,
    force: bool = False,
) -> Dict:
    """Return O_local / O_persist and the act policy for this batch."""
    D = detectability(error_ratio, thresh=detect_thresh)
    if force:
        D = max(D, 0.55)
    L = localizability(vimp, feat_idx)
    sep = separability(
        membership, X_live=X_live, prototypes_S=prototypes_S,
        feat_idx=feat_idx, mem_thr=mem_thr, margin=margin,
    )
    mass = actable_mass(
        membership, n_private=n_private, expand_min=expand_min, mem_thr=mem_thr,
    )
    vimp_top = float(np.asarray(vimp, float)[np.asarray(feat_idx, int)].sum()) if len(feat_idx) else 0.0
    Stab = footprint_stability_score(
        vimp_top, prev_vimp_top_mass, membership,
        prev_feat=prev_feat, feat_idx=feat_idx,
    )

    R, M = float(sep["R"]), float(mass["M"])
    # geometric mean keeps a single weak factor from being hidden by others
    O_local = float((max(D, 1e-6) * max(L, 1e-6) * max(R, 1e-6) * max(M, 1e-6)) ** 0.25)
    O_persist = _clip01(O_local * Stab)

    do_sep = bool(force or O_local >= float(tau_sep))
    do_exp = bool(do_sep and (O_local >= float(tau_exp)) and (mass["M_exp"] >= 0.45 or n_private >= max(4, expand_min // 2)))
    keep_offset = bool(O_persist >= float(tau_keep))

    if not do_sep:
        action = "hold"
    elif do_exp:
        action = "sep+expand"
    else:
        action = "sep_only"

    return {
        "O_local": O_local,
        "O_persist": O_persist,
        "D": D,
        "L": L,
        "R": R,
        "M": M,
        "Stab": Stab,
        "factors": {"sep": sep, "mass": mass},
        "vimp_top_mass": vimp_top,
        "tau": {"sep": float(tau_sep), "exp": float(tau_exp), "keep": float(tau_keep)},
        "do_sep": do_sep,
        "do_expand_ok": do_exp,
        "keep_offset": keep_offset,
        "action": action,
        # one-liner for logs / PO
        "read": (
            f"O_local={O_local:.2f} (D={D:.2f} L={L:.2f} R={R:.2f} M={M:.2f}) "
            f"O_persist={O_persist:.2f} → {action}"
        ),
    }


def opportunity_from_bound_labels(
    y_live: Sequence,
    private_live: Sequence,
) -> int:
    y = np.asarray(y_live).ravel()
    priv = list(private_live or [])
    if not priv:
        return 0
    return int(np.isin(y, priv).sum())
