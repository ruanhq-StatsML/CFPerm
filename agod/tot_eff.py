"""Tree-of-Thoughts over PO/AGOD efficiency policies.

Statistical empowerment (not LLM self-scores)
---------------------------------------------
- Expand: adapt ∈ {ref, probe, refit} × alpha ∈ {None, 0.5, 1/3} × freeze?
- Value:  rank_eff / mse_eff / expected_flops (from scorecards)
- Prune:  soft_burn decisions + budget + conditional Pareto

Thought = a policy node; V = scorecard metrics; prune = burn/Pareto/budget.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from agod.po_eff import (
    budget_ratio_refit_vs_probe,
    expected_adapt_flops,
    prefer_refit_on_budget,
)
from agod.soft_burn import burn_soft_weights


@dataclass
class Thought:
    """One node in the efficiency ToT."""

    adapt: str  # ref | probe | refit
    alpha: Optional[float] = None  # None=uniform, 0.5=√, 1/3=∛
    freeze: Optional[str] = None  # None | early | low_share
    rank_eff: float = float("nan")
    mse_eff: float = float("nan")
    rel_vs_uniform: float = float("nan")
    expected_flops: float = float("nan")
    burn: Optional[Dict[str, Any]] = None
    note: str = ""
    children: List["Thought"] = field(default_factory=list)

    @property
    def label(self) -> str:
        a = "unif" if self.alpha is None else (f"α={self.alpha:.2f}")
        f = f"+{self.freeze}" if self.freeze else ""
        return f"{self.adapt}|{a}{f}"


def _lex_key(t: Thought) -> Tuple[float, float, float]:
    """Prefer downstream win, then rank efficiency, then cheaper."""
    mse = t.mse_eff if np.isfinite(t.mse_eff) else -1e18
    # soft preference: positive mse_eff first
    mse_tier = 1.0 if mse > 0 else 0.0
    rank = t.rank_eff if np.isfinite(t.rank_eff) else -1e18
    cost = -(t.expected_flops if np.isfinite(t.expected_flops) else 1e18)
    return (mse_tier, mse, rank + 1e-12 * cost)


def expand_adapt_layer(
    *,
    gate_duty: float,
    n_batches: int,
    batch_size: int,
    n_control: int,
    rank_effs: Dict[str, float],
    mse_effs: Dict[str, float],
) -> List[Thought]:
    """Root children: choose adaptation mode under duty budget."""
    out: List[Thought] = []
    for adapt in ("ref", "probe", "refit"):
        flops = expected_adapt_flops(
            adapt,
            gate_duty=gate_duty,
            n_batches=n_batches,
            batch_size=batch_size,
            n_control=n_control,
        )
        out.append(
            Thought(
                adapt=adapt,
                rank_eff=float(rank_effs.get(adapt, float("nan"))),
                mse_eff=float(mse_effs.get(adapt, float("nan"))),
                expected_flops=flops,
                note=(
                    "prefer_refit_budget"
                    if adapt == "refit"
                    and prefer_refit_on_budget(gate_duty, n_control=n_control)
                    else ""
                ),
            )
        )
    # budget prune: drop probe if refit preferred and probe rank not better
    if prefer_refit_on_budget(gate_duty, n_control=n_control):
        refit = next(t for t in out if t.adapt == "refit")
        pruned = []
        for t in out:
            if t.adapt == "probe" and not (
                np.isfinite(t.rank_eff)
                and np.isfinite(refit.rank_eff)
                and t.rank_eff > refit.rank_eff
            ):
                continue  # prune always-on probe under low duty
            pruned.append(t)
        out = pruned
    return sorted(out, key=_lex_key, reverse=True)


def expand_weight_layer(
    parent: Thought,
    *,
    rel_sqrt: float,
    rel_cbrt: float,
    mse_eff_sqrt: float = float("nan"),
    mse_eff_cbrt: float = float("nan"),
) -> List[Thought]:
    """Children of an adapt node: uniform / √ / ∛ with soft-burn prune."""
    kids: List[Thought] = []
    # uniform always legal
    kids.append(
        Thought(
            adapt=parent.adapt,
            alpha=None,
            rank_eff=parent.rank_eff,
            mse_eff=0.0,  # baseline
            rel_vs_uniform=1.0,
            expected_flops=parent.expected_flops,
            burn={"decision": "KEEP_UNIFORM", "burn": False},
            note="safe leaf",
        )
    )
    for alpha, rel, me in (
        (0.5, rel_sqrt, mse_eff_sqrt),
        (1.0 / 3.0, rel_cbrt, mse_eff_cbrt),
    ):
        decision = burn_soft_weights(
            rel_gated_sqrt=rel_sqrt,
            rel_gated_cbrt=rel_cbrt,
            soft_win_cbrt_le_sqrt=bool(
                np.isfinite(rel_cbrt) and np.isfinite(rel_sqrt) and rel_cbrt <= rel_sqrt
            ),
            gate_already_on=True,
        )
        # Only keep BURN child matching this alpha; SOFTEN_ONLY keeps ∛ only
        keep = False
        if decision["burn"]:
            if decision["decision"] == "BURN_SQRT" and abs(alpha - 0.5) < 1e-9:
                keep = True
            if decision["decision"] == "BURN_CBRT" and abs(alpha - 1.0 / 3.0) < 1e-9:
                keep = True
        elif decision["decision"] == "SOFTEN_ONLY" and abs(alpha - 1.0 / 3.0) < 1e-9:
            keep = True  # leaf for damage control, not a win claim
        if not keep:
            continue
        kids.append(
            Thought(
                adapt=parent.adapt,
                alpha=alpha,
                rank_eff=parent.rank_eff,
                mse_eff=float(me),
                rel_vs_uniform=float(rel),
                expected_flops=parent.expected_flops,  # α ≠ FLOPs
                burn=decision,
                note=decision.get("reason", ""),
            )
        )
    parent.children = kids
    return kids


def beam_search_policy_tot(
    *,
    gate_duty: float,
    n_batches: int = 40,
    batch_size: int = 100,
    n_control: int = 1,
    rank_effs: Optional[Dict[str, float]] = None,
    mse_effs: Optional[Dict[str, float]] = None,
    rel_sqrt: float = 1.1,
    rel_cbrt: float = 1.05,
    mse_eff_sqrt: float = float("nan"),
    mse_eff_cbrt: float = float("nan"),
    beam_k: int = 3,
) -> Dict[str, Any]:
    """Two-level ToT: adapt → weight; return beam + best leaf."""
    rank_effs = rank_effs or {"ref": 0.0, "probe": 0.2, "refit": 1.0}
    mse_effs = mse_effs or {"ref": 0.0, "probe": -0.1, "refit": -0.05}
    level1 = expand_adapt_layer(
        gate_duty=gate_duty,
        n_batches=n_batches,
        batch_size=batch_size,
        n_control=n_control,
        rank_effs=rank_effs,
        mse_effs=mse_effs,
    )[:beam_k]
    leaves: List[Thought] = []
    for node in level1:
        kids = expand_weight_layer(
            node,
            rel_sqrt=rel_sqrt,
            rel_cbrt=rel_cbrt,
            mse_eff_sqrt=mse_eff_sqrt,
            mse_eff_cbrt=mse_eff_cbrt,
        )
        leaves.extend(kids)
    leaves_sorted = sorted(leaves, key=_lex_key, reverse=True)
    best = leaves_sorted[0] if leaves_sorted else None
    return {
        "budget_ratio_refit_vs_probe": budget_ratio_refit_vs_probe(
            gate_duty, n_control=n_control
        ),
        "beam_adapt": [t.label for t in level1],
        "leaves": [
            {
                "label": t.label,
                "rank_eff": t.rank_eff,
                "mse_eff": t.mse_eff,
                "rel": t.rel_vs_uniform,
                "flops": t.expected_flops,
                "burn": (t.burn or {}).get("decision"),
                "note": t.note,
            }
            for t in leaves_sorted
        ],
        "best": None
        if best is None
        else {
            "label": best.label,
            "burn": (best.burn or {}).get("decision"),
            "note": best.note,
        },
        "reading": (
            f"best={best.label if best else None}; "
            f"duty={gate_duty:.2f}; "
            "V=(mse_eff tier, mse_eff, rank_eff); α≠FLOPs"
        ),
    }


def tot_from_pack_metrics(
    *,
    gate_duty: float,
    rank_refit: float,
    rank_probe: float,
    mse_refit: float,
    mse_probe: float,
    rel_sqrt: float,
    rel_cbrt: float,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Convenience: plug scorecard numbers into ToT."""
    return beam_search_policy_tot(
        gate_duty=gate_duty,
        rank_effs={"ref": 0.0, "probe": rank_probe, "refit": rank_refit},
        mse_effs={"ref": 0.0, "probe": mse_probe, "refit": mse_refit},
        rel_sqrt=rel_sqrt,
        rel_cbrt=rel_cbrt,
        mse_eff_sqrt=mse_probe,  # rough: weight layer inherits adapt mse signal
        mse_eff_cbrt=mse_refit,
        **kwargs,
    )
