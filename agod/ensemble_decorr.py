"""Ensemble decorrelation under high modality correlation.

When pairwise grad cosines are high, modalities are *not* independent voters —
they form a near-rank-1 ensemble. Skipping alone is blunt; the control plane
should **decompose roles** and **decorrelate** step sizes:

  Gram G[i,j] = cos(g_i, g_j)           # alignment geometry
  effective_rank(G) ↓  ⇒  collinear ensemble
  leader      = argmax_m α_m · (1−ρ_m)_+  # owns shared update direction
  diversifier = high residual uniqueness after removing leader direction
  redundant   = collinear with leader / low residual → damp adapt LR

LR mapping (FWD always on; only BWD step sizes reshape):
  LR_m ∝ α_m · role_gain_m
  role_gain: leader↑, diversifier mid (residual), redundant→β floor

Justification: paying proj-BWD on a redundant tower retrains the *same*
direction; diversifiers buy orthogonal residual mass — that is the only
ensemble member worth keeping when correlation is high.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .compute_ratio import modality_redundancy, pair_cos_matrix, unique_adapt_mass


ROLE_NAMES = ("leader", "diversifier", "redundant")


def gram_effective_rank(pair_cos: Mapping[str, float], mods: Sequence[str]) -> dict:
    """Spectral summary of the pairwise cosine Gram matrix.

    effective_rank = exp(H(λ/Σλ)) ∈ [1, |M|]  (Roy & Vetterli).
    Near 1 ⇒ one shared direction dominates ⇒ decorrelate / role-split.
    """
    mods = list(mods)
    g = pair_cos_matrix(pair_cos, mods)
    # symmetrize + clip for numeric PSD-ish
    g = 0.5 * (g + g.T)
    np.fill_diagonal(g, 1.0)
    evals = np.linalg.eigvalsh(g)
    evals = np.clip(evals, 0.0, None)
    s = float(evals.sum())
    if s <= 1e-12:
        return {
            "evals": evals.tolist(),
            "effective_rank": float(len(mods)),
            "participation_ratio": float(len(mods)),
            "top_frac": 0.0,
        }
    p = evals / s
    p = p[p > 1e-12]
    h = float(-np.sum(p * np.log(p)))
    erank = float(np.exp(h))
    # participation ratio (alternative rank proxy)
    pr = float((s * s) / max(float(np.sum(evals * evals)), 1e-12))
    top_frac = float(evals.max() / s)
    return {
        "evals": [float(x) for x in evals],
        "effective_rank": erank,
        "participation_ratio": pr,
        "top_frac": top_frac,
        "n_mods": len(mods),
    }


def residual_uniqueness(
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    leader: str,
) -> dict[str, float]:
    """How much of each modality is *not* explained by the leader direction.

    r_m = (1 − max(cos(g_m, g_leader), 0))_+
    leader itself gets r=0 (it *is* the shared direction).
    """
    mods = list(mods)
    mat = pair_cos_matrix(pair_cos, mods)
    idx = {m: i for i, m in enumerate(mods)}
    li = idx[leader]
    out = {}
    for m in mods:
        if m == leader:
            out[m] = 0.0
            continue
        c = float(mat[idx[m], li])
        out[m] = float(max(1.0 - max(c, 0.0), 0.0))
    return out


def assign_ensemble_roles(
    alpha: Mapping[str, float],
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
    erank_collinear: float = 1.55,
) -> dict:
    """Decompose modalities into ensemble roles under correlation geometry.

    Trigger: mean redundancy ≥ rho_high **or** effective_rank ≤ erank_collinear.
    When *not* triggered (low correlation / full-rank ensemble), every modality
    is treated as an independent voter → role=diversifier, no damp.

    Returns roles, unique/residual masses, spectral stats, and whether
    decorrelation mode is active.
    """
    mods = list(mods)
    rho = modality_redundancy(pair_cos, mods)
    uniq = unique_adapt_mass(alpha, rho, mods, eta=eta)
    spec = gram_effective_rank(pair_cos, mods)
    mean_rho = float(np.mean([rho[m] for m in mods]))
    decorr = bool(mean_rho >= float(rho_high) or spec["effective_rank"] <= float(erank_collinear))

    # leader = highest unique mass (attribution after redundancy discount)
    leader = max(mods, key=lambda m: uniq[m])
    resid = residual_uniqueness(pair_cos, mods, leader)
    # residual *mass* for budgeting: α · residual uniqueness
    resid_mass_raw = {m: float(alpha[m]) * float(resid[m]) for m in mods}
    s = sum(resid_mass_raw.values())
    resid_mass = (
        {m: resid_mass_raw[m] / s for m in mods}
        if s > 1e-12
        else {m: 0.0 for m in mods}
    )

    roles: dict[str, str] = {}
    if not decorr:
        roles = {m: "diversifier" for m in mods}
    else:
        roles[leader] = "leader"
        for m in mods:
            if m == leader:
                continue
            # diversifier if residual uniqueness clears threshold
            if resid[m] >= float(residual_div):
                roles[m] = "diversifier"
            else:
                roles[m] = "redundant"

    return {
        "decorr_active": decorr,
        "leader": leader,
        "roles": roles,
        "redundancy": rho,
        "unique_mass": uniq,
        "residual": resid,
        "residual_mass": resid_mass,
        "mean_redundancy": mean_rho,
        "spectral": spec,
        "rho_high": float(rho_high),
        "residual_div": float(residual_div),
        "erank_collinear": float(erank_collinear),
        "eta": float(eta),
    }


def role_lr_gains(
    roles: Mapping[str, str],
    mods: Sequence[str],
    *,
    residual: Mapping[str, float] | None = None,
    leader_gain: float = 1.35,
    diversifier_gain: float = 1.00,
    redundant_gain: float = 0.25,
    residual_boost: float = 0.50,
) -> dict[str, float]:
    """Map ensemble roles → multiplicative LR gains (clipped ≥0).

    Diversifiers get an extra residual_boost · r_m so more-orthogonal members
    step harder than near-collinear ones still labeled diversifier.
    """
    mods = list(mods)
    residual = residual or {m: 0.0 for m in mods}
    out = {}
    for m in mods:
        role = roles.get(m, "diversifier")
        if role == "leader":
            out[m] = float(leader_gain)
        elif role == "redundant":
            out[m] = float(redundant_gain)
        else:
            out[m] = float(diversifier_gain + residual_boost * float(residual.get(m, 0.0)))
    return out


def soft_decorr_lr(
    alpha: Mapping[str, float],
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
    leader_gain: float = 1.35,
    diversifier_gain: float = 1.00,
    redundant_gain: float = 0.25,
    residual_boost: float = 0.50,
    gain: Mapping[str, float] | None = None,
) -> dict:
    """Soft α→LR with ensemble decorrelation gains when correlation is high.

    When decorr is inactive (low ρ / high effective rank), reduces to soft LR
    (role gains ≈ 1 with mild residual boost only if roles are all diversifier
    and residual_boost is applied — here we force gain=1 when inactive).
    """
    from .lr_controller import alpha_to_lr  # local import avoids cycle

    mods = list(mods)
    decomp = assign_ensemble_roles(
        alpha,
        pair_cos,
        mods,
        eta=eta,
        rho_high=rho_high,
        residual_div=residual_div,
    )
    if decomp["decorr_active"]:
        rg = role_lr_gains(
            decomp["roles"],
            mods,
            residual=decomp["residual"],
            leader_gain=leader_gain,
            diversifier_gain=diversifier_gain,
            redundant_gain=redundant_gain,
            residual_boost=residual_boost,
        )
    else:
        rg = {m: 1.0 for m in mods}

    # fold optional external gain
    gain = gain or {m: 1.0 for m in mods}
    combined = {m: float(rg[m]) * float(gain[m]) for m in mods}
    lr = alpha_to_lr(alpha, mods, beta=beta, gain=combined)
    return {
        "lr": lr,
        "role_gain": rg,
        "decomp": decomp,
    }


def characterize_ensemble_traj(
    traj: Sequence[Mapping],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
) -> dict:
    """Summarize role assignments / spectral rank over a window trajectory."""
    mods = list(mods)
    if not traj:
        return {}
    eranks, mean_rhos, decorr_flags = [], [], []
    role_counts = {r: 0 for r in ROLE_NAMES}
    leaders: dict[str, int] = {}
    for r in traj:
        alpha = r.get("alpha") or {}
        pair = r.get("pair_cos") or {}
        if not pair and "mean_pair_cos" in r:
            mpc = float(r["mean_pair_cos"])
            pair = {
                f"{a}|{b}": mpc
                for i, a in enumerate(mods)
                for b in mods[i + 1 :]
            }
        decomp = assign_ensemble_roles(
            alpha, pair, mods, eta=eta, rho_high=rho_high, residual_div=residual_div
        )
        eranks.append(decomp["spectral"]["effective_rank"])
        mean_rhos.append(decomp["mean_redundancy"])
        decorr_flags.append(1.0 if decomp["decorr_active"] else 0.0)
        leaders[decomp["leader"]] = leaders.get(decomp["leader"], 0) + 1
        for m in mods:
            role_counts[decomp["roles"][m]] = role_counts.get(decomp["roles"][m], 0) + 1
    n = max(len(traj), 1)
    return {
        "n_windows": len(traj),
        "mean_redundancy": float(np.nanmean(mean_rhos)),
        "mean_effective_rank": float(np.nanmean(eranks)),
        "frac_decorr_active": float(np.nanmean(decorr_flags)),
        "leader_counts": leaders,
        "role_counts": role_counts,
        "role_frac": {k: float(v) / (n * len(mods)) for k, v in role_counts.items()},
        "rho_high": float(rho_high),
        "residual_div": float(residual_div),
        "eta": float(eta),
    }


# ---------------------------------------------------------------------------
# Correlation buckets + R-run trajectory diagnostics (online decorr eval)
# ---------------------------------------------------------------------------

BUCKET_NAMES = ("low", "mid", "high")


def corr_bucket(
    mean_rho: float,
    effective_rank: float,
    *,
    rho_high: float = 0.55,
    rho_mid: float = 0.35,
    erank_collinear: float = 1.55,
    erank_mid: float = 2.20,
) -> str:
    """Map a window's (ρ̄, erank) → low / mid / high correlation bucket.

    high: trigger-aligned (ρ̄≥ρ_high OR erank≤erank_collinear)
    low:  clearly independent (ρ̄<ρ_mid AND erank>erank_mid)
    mid:  everything else (borderline / mixed)
    """
    if float(mean_rho) >= float(rho_high) or float(effective_rank) <= float(erank_collinear):
        return "high"
    if float(mean_rho) < float(rho_mid) and float(effective_rank) > float(erank_mid):
        return "low"
    return "mid"


def _pair_from_row(row: Mapping, mods: Sequence[str]) -> dict[str, float]:
    mods = list(mods)
    pair = dict(row.get("pair_cos") or {})
    if pair:
        return {k: float(v) for k, v in pair.items()}
    mpc = float(row.get("mean_pair_cos", row.get("mean_redundancy", 0.0)) or 0.0)
    return {
        f"{a}|{b}": mpc
        for i, a in enumerate(mods)
        for b in mods[i + 1 :]
    }


def expand_role_windows(
    traj: Sequence[Mapping],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
    erank_collinear: float = 1.55,
) -> list[dict]:
    """Per-window role / bucket / LR expansion (for R-run & swap stats)."""
    mods = list(mods)
    out = []
    for i, row in enumerate(traj):
        alpha = {m: float((row.get("alpha") or {}).get(m, 0.0)) for m in mods}
        # tolerate missing α → uniform
        s = sum(alpha.values())
        if s <= 1e-12:
            alpha = {m: 1.0 / len(mods) for m in mods}
        pair = _pair_from_row(row, mods)
        packed = soft_decorr_lr(
            alpha,
            pair,
            mods,
            eta=eta,
            rho_high=rho_high,
            residual_div=residual_div,
        )
        d = packed["decomp"]
        erank = float(d["spectral"]["effective_rank"])
        mean_rho = float(d["mean_redundancy"])
        bucket = corr_bucket(
            mean_rho,
            erank,
            rho_high=rho_high,
            erank_collinear=erank_collinear,
        )
        leader_alpha = max(mods, key=lambda m: alpha[m])
        leader_u = d["leader"]
        lr = packed["lr"]
        lr_vals = [float(lr[m]) for m in mods]
        out.append(
            {
                "t": int(row.get("t", i)),
                "bucket": bucket,
                "decorr_active": bool(d["decorr_active"]),
                "mean_redundancy": mean_rho,
                "effective_rank": erank,
                "top_frac": float(d["spectral"].get("top_frac", float("nan"))),
                "alpha": alpha,
                "unique_mass": d["unique_mass"],
                "residual": d["residual"],
                "roles": d["roles"],
                "leader_alpha": leader_alpha,
                "leader_u": leader_u,
                "leader_swap": bool(leader_alpha != leader_u),
                "role_gain": packed["role_gain"],
                "lr": {m: float(lr[m]) for m in mods},
                "lr_ratio": float(max(lr_vals) / max(min(lr_vals), 1e-9)),
                "acc_lift": (
                    float(row["acc_lift"])
                    if row.get("acc_lift") is not None
                    else None
                ),
            }
        )
    return out


def extract_r_runs(
    windows: Sequence[Mapping],
    mods: Sequence[str],
    *,
    min_streak: int = 2,
) -> list[dict]:
    """Extract redundant streaks (R-runs) per modality.

    An R-run is a contiguous span where role==redundant AND decorr_active,
    of length ≥ min_streak. Low-corr windows (decorr off) break streaks.
    """
    mods = list(mods)
    runs: list[dict] = []
    for m in mods:
        start = None
        length = 0
        for i, w in enumerate(windows):
            is_r = bool(w.get("decorr_active")) and (w.get("roles") or {}).get(m) == "redundant"
            if is_r:
                if start is None:
                    start = i
                    length = 1
                else:
                    length += 1
            else:
                if start is not None and length >= int(min_streak):
                    span = list(windows[start : start + length])
                    lifts = [x["acc_lift"] for x in span if x.get("acc_lift") is not None]
                    runs.append(
                        {
                            "modality": m,
                            "t_start": int(span[0]["t"]),
                            "t_end": int(span[-1]["t"]),
                            "length": int(length),
                            "mean_rho": float(np.mean([x["mean_redundancy"] for x in span])),
                            "mean_erank": float(np.mean([x["effective_rank"] for x in span])),
                            "mean_lr": float(np.mean([x["lr"][m] for x in span])),
                            "mean_residual": float(np.mean([x["residual"][m] for x in span])),
                            "mean_acc_lift": float(np.mean(lifts)) if lifts else None,
                            "buckets": [x["bucket"] for x in span],
                        }
                    )
                start = None
                length = 0
        if start is not None and length >= int(min_streak):
            span = list(windows[start : start + length])
            lifts = [x["acc_lift"] for x in span if x.get("acc_lift") is not None]
            runs.append(
                {
                    "modality": m,
                    "t_start": int(span[0]["t"]),
                    "t_end": int(span[-1]["t"]),
                    "length": int(length),
                    "mean_rho": float(np.mean([x["mean_redundancy"] for x in span])),
                    "mean_erank": float(np.mean([x["effective_rank"] for x in span])),
                    "mean_lr": float(np.mean([x["lr"][m] for x in span])),
                    "mean_residual": float(np.mean([x["residual"][m] for x in span])),
                    "mean_acc_lift": float(np.mean(lifts)) if lifts else None,
                    "buckets": [x["bucket"] for x in span],
                }
            )
    return runs


def _role_transition_matrix(windows: Sequence[Mapping], mods: Sequence[str]) -> dict:
    """P(role_{t+1}|role_t) pooled over modalities (decorr-on steps only)."""
    labels = list(ROLE_NAMES)
    counts = {a: {b: 0 for b in labels} for a in labels}
    for m in mods:
        prev = None
        for w in windows:
            if not w.get("decorr_active"):
                prev = None
                continue
            cur = (w.get("roles") or {}).get(m)
            if cur not in counts:
                prev = None
                continue
            if prev in counts:
                counts[prev][cur] += 1
            prev = cur
    out = {}
    for a in labels:
        tot = sum(counts[a].values())
        out[a] = {
            b: (float(counts[a][b]) / tot if tot else 0.0) for b in labels
        }
        out[a]["_n"] = int(tot)
    return out


def characterize_buckets_and_rruns(
    traj: Sequence[Mapping],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
    erank_collinear: float = 1.55,
    min_streak: int = 2,
) -> dict:
    """Bucket windows by correlation regime + R-run / leader-swap diagnostics.

    Evaluation protocol:
      - low bucket: soft_decorr should match soft (γ=1); R-run freq ≈ 0
      - high bucket: decorr on; report R-run rate, swap rate, LR_L/LR_R
      - mid bucket: borderline — report separately, do not pool into high
    """
    mods = list(mods)
    windows = expand_role_windows(
        traj,
        mods,
        eta=eta,
        rho_high=rho_high,
        residual_div=residual_div,
        erank_collinear=erank_collinear,
    )
    if not windows:
        return {}

    buckets: dict[str, list] = {b: [] for b in BUCKET_NAMES}
    for w in windows:
        buckets[w["bucket"]].append(w)

    def _bucket_stats(rows: list[dict]) -> dict:
        if not rows:
            return {
                "n_windows": 0,
                "frac_windows": 0.0,
                "mean_rho": float("nan"),
                "mean_erank": float("nan"),
                "frac_decorr_active": float("nan"),
                "frac_leader_swap": float("nan"),
                "role_frac": {r: float("nan") for r in ROLE_NAMES},
                "mean_lr_ratio": float("nan"),
                "mean_acc_lift": float("nan"),
            }
        n = len(rows)
        role_counts = {r: 0 for r in ROLE_NAMES}
        for w in rows:
            for m in mods:
                role_counts[w["roles"][m]] += 1
        lifts = [w["acc_lift"] for w in rows if w.get("acc_lift") is not None]
        return {
            "n_windows": n,
            "frac_windows": float(n / len(windows)),
            "mean_rho": float(np.mean([w["mean_redundancy"] for w in rows])),
            "mean_erank": float(np.mean([w["effective_rank"] for w in rows])),
            "frac_decorr_active": float(np.mean([1.0 if w["decorr_active"] else 0.0 for w in rows])),
            "frac_leader_swap": float(np.mean([1.0 if w["leader_swap"] else 0.0 for w in rows])),
            "role_frac": {
                r: float(role_counts[r]) / (n * len(mods)) for r in ROLE_NAMES
            },
            "mean_lr_ratio": float(np.mean([w["lr_ratio"] for w in rows])),
            "mean_acc_lift": float(np.mean(lifts)) if lifts else float("nan"),
        }

    bucket_stats = {b: _bucket_stats(buckets[b]) for b in BUCKET_NAMES}
    r_runs = extract_r_runs(windows, mods, min_streak=min_streak)

    # R-run frequency: runs per decorr-on window (and per high-bucket window)
    n_decorr = sum(1 for w in windows if w["decorr_active"])
    n_high = len(buckets["high"])
    n_r_windows = 0
    for w in windows:
        if not w["decorr_active"]:
            continue
        if any(w["roles"][m] == "redundant" for m in mods):
            n_r_windows += 1

    run_lengths = [r["length"] for r in r_runs]
    rrun_stats = {
        "min_streak": int(min_streak),
        "n_runs": len(r_runs),
        "runs_per_window": float(len(r_runs) / max(len(windows), 1)),
        "runs_per_decorr_window": float(len(r_runs) / max(n_decorr, 1)),
        "runs_per_high_bucket_window": float(len(r_runs) / max(n_high, 1)),
        "frac_windows_with_any_R": float(n_r_windows / max(n_decorr, 1)) if n_decorr else 0.0,
        "mean_run_length": float(np.mean(run_lengths)) if run_lengths else 0.0,
        "max_run_length": int(max(run_lengths)) if run_lengths else 0,
        "by_modality": {
            m: {
                "n_runs": sum(1 for r in r_runs if r["modality"] == m),
                "mean_length": float(
                    np.mean([r["length"] for r in r_runs if r["modality"] == m])
                    if any(r["modality"] == m for r in r_runs)
                    else 0.0
                ),
            }
            for m in mods
        },
        "runs": r_runs,
    }

    # length≥1 R events (single-window hits) vs sticky runs (length≥min_streak)
    hits = extract_r_runs(windows, mods, min_streak=1)
    sticky = [s for s in hits if s["length"] >= min_streak]
    rrun_stats["n_r_events_len1"] = len(hits)
    rrun_stats["n_r_mod_windows"] = int(sum(h["length"] for h in hits))

    return {
        "n_windows": len(windows),
        "mods": mods,
        "bucket_edges": {
            "rho_high": float(rho_high),
            "rho_mid": 0.35,
            "erank_collinear": float(erank_collinear),
            "erank_mid": 2.20,
        },
        "buckets": bucket_stats,
        "r_runs": rrun_stats,
        "n_sticky_r_runs": len(sticky),
        "role_transitions_decorr_on": _role_transition_matrix(windows, mods),
        "overall": {
            "mean_rho": float(np.mean([w["mean_redundancy"] for w in windows])),
            "mean_erank": float(np.mean([w["effective_rank"] for w in windows])),
            "frac_decorr_active": float(
                np.mean([1.0 if w["decorr_active"] else 0.0 for w in windows])
            ),
            "frac_leader_swap": float(
                np.mean([1.0 if w["leader_swap"] else 0.0 for w in windows])
            ),
        },
        "windows": windows,
    }
