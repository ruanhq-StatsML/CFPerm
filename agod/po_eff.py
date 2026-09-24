"""PO-risk adaptation efficiency: quality / relative FLOPs.

Pivot away from adjacent-chunk transfer probes.  Here the question is:

  On OnlineRFPerm reject batches, which PO scorer (frozen ref / rolling
  probe / recent-window refit) buys the most *hard-row ranking* or
  *sig-only MSE drop* per unit of adaptation compute?

FLOPs model (order-of-magnitude; not wall-clock)
-----------------------------------------------
- ``ref``:     uses frozen burn-in f_ref — **0** re-fits at reject
- ``probe``:   re-fits a probe each stream step (always-on cost)
- ``refit``:   re-fits μ0 on recent control **only on reject**

relative_flops ≈ n_fits · batch_size · n_trees · depth_proxy

Efficiency
----------
- rank_eff  = ΔSpearman vs ref  /  relative_flops_million
- mse_eff   = (MSE_unif_sig − MSE_mode_sig) / relative_flops_million

Positive mse_eff ⇒ cheaper error reduction; negative ⇒ paid FLOPs to
*worsen* sig MSE (common when ranking wins but IPTW hurts).
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np

# RF fit proxy: n_samples * n_trees * depth (same order as board probes)
_DEFAULT_TREES = 30
_DEFAULT_DEPTH = 6


def fit_flops(
    n_samples: int,
    *,
    n_trees: int = _DEFAULT_TREES,
    depth: int = _DEFAULT_DEPTH,
) -> float:
    return float(max(int(n_samples), 1) * n_trees * depth)


def mode_relative_flops(
    mode: str,
    *,
    n_batches: int,
    n_reject: int,
    batch_size: int,
    n_control: int = 1,
) -> float:
    """Adaptation FLOPs proxy for one PO mode over a stream (realized rejects)."""
    fit = fit_flops(batch_size)
    m = mode.lower().replace("_po", "")
    if m in ("ref", "uniform", "dre"):
        # ref: no re-fit; uniform/dre: no PO-learner re-fit counted here
        return 1.0  # ε floor so ratios stay finite
    if m == "probe":
        # probe updates every step after burn-in ≈ n_batches fits
        return float(max(n_batches, 1) * fit)
    if m == "refit":
        # one fit on ~n_control batches of size batch_size, per reject
        win = fit_flops(batch_size * max(int(n_control), 1))
        return float(max(n_reject, 1) * win)
    return float("nan")


def expected_adapt_flops(
    mode: str,
    *,
    gate_duty: float,
    n_batches: int,
    batch_size: int,
    n_control: int = 1,
) -> float:
    """Ex-ante adaptation budget: E[flops] = f(duty), not realized n_reject.

    Statistical justify
    -------------------
    Realized ``n_reject`` is a post-hoc count.  Planning / comparing gates
    needs **expected** cost under a Bernoulli reject rate ``duty``:

      E[flops_refit]  = n_batches · duty · fit(window)
      E[flops_probe]  = n_batches · fit(batch)          # duty-invariant
      E[flops_ref]    ≈ ε

    ``budget_ratio = E[refit] / E[probe] ≈ duty · (n_control)`` when fit
    scales with window size — low duty ⇒ refit is the cheap PO path.
    """
    duty = float(np.clip(gate_duty, 0.0, 1.0))
    fit = fit_flops(batch_size)
    m = mode.lower().replace("_po", "")
    if m in ("ref", "uniform", "dre"):
        return 1.0
    if m == "probe":
        return float(max(n_batches, 1) * fit)
    if m == "refit":
        win = fit_flops(batch_size * max(int(n_control), 1))
        return float(max(n_batches, 1) * duty * win)
    return float("nan")


def budget_ratio_refit_vs_probe(
    gate_duty: float,
    *,
    n_control: int = 1,
) -> float:
    """E[refit]/E[probe] ≈ duty · n_control (same batch_size scale)."""
    duty = float(np.clip(gate_duty, 0.0, 1.0))
    return float(duty * max(int(n_control), 1))


def rank_efficiency(
    spearman: float,
    spearman_ref: float,
    relative_flops: float,
    *,
    scale: float = 1e6,
) -> float:
    if not np.isfinite(spearman) or not np.isfinite(spearman_ref):
        return float("nan")
    if not np.isfinite(relative_flops) or relative_flops <= 0:
        return float("nan")
    return float((spearman - spearman_ref) / (relative_flops / scale))


def mse_efficiency(
    mse_mode: float,
    mse_uniform: float,
    relative_flops: float,
    *,
    scale: float = 1e6,
) -> float:
    """Positive ⇒ mode reduces sig MSE per million FLOPs vs uniform."""
    if not np.isfinite(mse_mode) or not np.isfinite(mse_uniform):
        return float("nan")
    if not np.isfinite(relative_flops) or relative_flops <= 0:
        return float("nan")
    return float((mse_uniform - mse_mode) / (relative_flops / scale))


def scorecard_from_dataset_block(
    name: str,
    block: Mapping[str, Any],
    *,
    n_batches: int,
    batch_size: int,
    n_control: int = 1,
) -> Dict[str, Any]:
    """Build per-mode efficiency rows from one ``summary.json`` dataset block."""
    results = block.get("results") or {}
    # quality lives under any mode that recorded po_quality (usually uniform)
    pq = None
    for _m, payload in results.items():
        if isinstance(payload, dict) and payload.get("po_quality"):
            pq = payload["po_quality"]
            break
    if not pq:
        return {"dataset": name, "ok": False, "reason": "no po_quality"}

    n_reject = int(pq.get("n_reject") or 0)
    ref_sp = float((pq.get("ref") or {}).get("spearman", float("nan")))
    rows: List[Dict[str, Any]] = []
    unif = results.get("uniform") or {}
    unif_mse = float(unif.get("mse_mean_sig", float("nan")))
    gate_duty = float(unif.get("gate_duty", float("nan")))
    if not np.isfinite(gate_duty) and n_batches > 0:
        # fallback: realized reject rate
        gate_duty = float(n_reject) / float(max(n_batches, 1))

    for label, key in (
        ("ref", "ref"),
        ("probe", "probe"),
        ("refit", "refit"),
    ):
        q = pq.get(key) or {}
        sp = float(q.get("spearman", float("nan")))
        flops = mode_relative_flops(
            label,
            n_batches=n_batches,
            n_reject=n_reject,
            batch_size=batch_size,
            n_control=n_control,
        )
        e_flops = expected_adapt_flops(
            label,
            gate_duty=gate_duty if np.isfinite(gate_duty) else 0.0,
            n_batches=n_batches,
            batch_size=batch_size,
            n_control=n_control,
        )
        mode_key = f"{label}_po" if label != "ref" else "ref_po"
        if label == "ref":
            mode_key = "ref_po"
        mse = float((results.get(mode_key) or {}).get("mse_mean_sig", float("nan")))
        rows.append(
            {
                "mode": label,
                "spearman": sp,
                "delta_spearman_vs_ref": (
                    float(sp - ref_sp) if np.isfinite(sp) and np.isfinite(ref_sp) else float("nan")
                ),
                "mse_mean_sig": mse,
                "relative_flops": flops,
                "expected_flops": e_flops,
                "rank_eff": rank_efficiency(sp, ref_sp, flops),
                "mse_eff": mse_efficiency(mse, unif_mse, flops),
                "rank_eff_expected": rank_efficiency(sp, ref_sp, e_flops),
                "mse_eff_expected": mse_efficiency(mse, unif_mse, e_flops),
            }
        )

    # Best rank_eff among modes that actually spend FLOPs (probe/refit)
    spend = [r for r in rows if r["mode"] in ("probe", "refit") and np.isfinite(r["rank_eff"])]
    best_rank = max(spend, key=lambda r: r["rank_eff"])["mode"] if spend else None
    spend_m = [r for r in rows if r["mode"] in ("probe", "refit") and np.isfinite(r["mse_eff"])]
    best_mse = max(spend_m, key=lambda r: r["mse_eff"])["mode"] if spend_m else None
    spend_e = [
        r
        for r in rows
        if r["mode"] in ("probe", "refit") and np.isfinite(r["rank_eff_expected"])
    ]
    best_rank_E = (
        max(spend_e, key=lambda r: r["rank_eff_expected"])["mode"] if spend_e else None
    )

    return {
        "dataset": name,
        "ok": True,
        "n_reject": n_reject,
        "n_batches": n_batches,
        "batch_size": batch_size,
        "gate_duty": gate_duty,
        "budget_ratio_refit_vs_probe": budget_ratio_refit_vs_probe(
            gate_duty if np.isfinite(gate_duty) else 0.0, n_control=n_control
        ),
        "uniform_mse_sig": unif_mse,
        "ref_spearman": ref_sp,
        "modes": rows,
        "best_rank_eff_mode": best_rank,
        "best_mse_eff_mode": best_mse,
        "best_rank_eff_expected_mode": best_rank_E,
        "reading": _reading(rows, best_rank, best_mse, gate_duty),
    }


def _reading(
    rows: Sequence[Dict[str, Any]],
    best_rank: Optional[str],
    best_mse: Optional[str],
    gate_duty: float = float("nan"),
) -> str:
    by = {r["mode"]: r for r in rows}
    refit = by.get("refit") or {}
    probe = by.get("probe") or {}
    duty_s = f"duty={gate_duty:.2f}" if np.isfinite(gate_duty) else "duty=?"
    if best_rank == "probe" and (probe.get("mse_eff") or 0) < 0:
        return f"{duty_s}: probe buys ranking cheaply but IPTW hurts sig MSE"
    if best_rank == "refit" and (refit.get("rank_eff") or 0) > (probe.get("rank_eff") or -1e9):
        return f"{duty_s}: refit wins rank_eff (reject-only vs always-on probe)"
    if best_mse == "refit" and (refit.get("mse_eff") or 0) > 0:
        return f"{duty_s}: refit reduces sig MSE per reject-FLOP vs uniform"
    if best_mse is None or all((r.get("mse_eff") or 0) <= 0 for r in rows if r["mode"] != "ref"):
        return f"{duty_s}: adaptation spends FLOPs without sig-MSE win — prefer ref/uniform"
    return f"{duty_s}: rank_eff→{best_rank}; mse_eff→{best_mse}"


def scorecard_from_summary(summary: Mapping[str, Any]) -> Dict[str, Any]:
    n_batches = int(summary.get("n_batches") or 40)
    batch_size = int(summary.get("batch_size") or 100)
    n_control = int(summary.get("n_control") or 1)
    cards = []
    for name, block in (summary.get("datasets") or {}).items():
        if not isinstance(block, dict):
            continue
        cards.append(
            scorecard_from_dataset_block(
                name,
                block,
                n_batches=n_batches,
                batch_size=batch_size,
                n_control=n_control,
            )
        )
    ok = [c for c in cards if c.get("ok")]
    # Aggregate: count best modes
    from collections import Counter

    br = Counter(c["best_rank_eff_mode"] for c in ok if c.get("best_rank_eff_mode"))
    bm = Counter(c["best_mse_eff_mode"] for c in ok if c.get("best_mse_eff_mode"))
    be = Counter(
        c["best_rank_eff_expected_mode"]
        for c in ok
        if c.get("best_rank_eff_expected_mode")
    )
    mean_duty = float(
        np.mean([c["gate_duty"] for c in ok if np.isfinite(c.get("gate_duty", np.nan))])
    ) if ok else float("nan")
    return {
        "n_datasets": len(ok),
        "batch_size": batch_size,
        "n_batches": n_batches,
        "n_control": n_control,
        "mean_gate_duty": mean_duty,
        "mean_budget_ratio_refit_vs_probe": (
            float(budget_ratio_refit_vs_probe(mean_duty, n_control=n_control))
            if np.isfinite(mean_duty)
            else float("nan")
        ),
        "best_rank_eff_counts": dict(br),
        "best_mse_eff_counts": dict(bm),
        "best_rank_eff_expected_counts": dict(be),
        "cards": cards,
        "headline": (
            f"mean duty={mean_duty:.3f} ⇒ E[refit]/E[probe]≈"
            f"{budget_ratio_refit_vs_probe(mean_duty, n_control=n_control):.3f}; "
            f"rank_eff wins {dict(br)}; expected-rank wins {dict(be)}; "
            f"mse_eff wins {dict(bm)}."
        ),
    }
