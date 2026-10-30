"""Product next-step logic for PO-risk Var(u) — straightforward triage.

Pipeline (one direction only)
-----------------------------
  1. Fit φ, τ̂  (PermuCATE RF or CausalForestDML)
  2. Unit risk   u = (φ − τ̂)²
  3. Mag×Var     μ = E[u],  s = Std(sub-env mean u)     ← shape
  4. Var(u)      CF leaf-size-debiased (default)         ← model uncertainty
  5. Y-flip      (text integrity boards only)            ← outcome break
  6. Actuator    typed medicine from the gate below
  7. Honesty     report bias_ratio = B_leaf / V_naive

Gate (product) — Var(u) confirms, does not solo-freeze
-------------------------------------------------------
  Primary sensors: Mag×Var (μ, s, J) + Y-flip (integrity).
  Var(u) CF-debiased is the *honesty / confirmation* channel:
    - removes leaf-size bias that would masquerade as malignant
    - relative V↑ vs REF supports a freeze only *with* shape/outcome evidence

  quiet / mild μ, no Y-flip, no Mag×Var spike     → hold / conformal
  Y-flip>0  OR  (Mag×Var spike)  OR  (V↑_rel AND Mag×Var/Y evidence)
                                              → VAI freeze

Default estimator for step 4 is always CF-debiased + leaf merge.
BCF only when a Bayesian pipeline is already in place.
Bootstrap is fallback / sanity check — never the product default.

Why this is the next step after Mag×Var
---------------------------------------
Mag×Var separates benign vs malignant *shape* of u.
Var(u) from the forest removes a different failure mode: leaf-size bias in
Var(τ̂) that inflates empirical Var(u) and can *masquerade as malignant*
on a quiet/benign LIVE. Debias is the product honesty check — not a
standalone absolute threshold on V.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Optional


@dataclass
class TriageDecision:
    """Typed actuator call from Mag×Var + Var(u) + optional Y-flip."""

    call: str  # quiet | benign | malignant
    actuator: str  # hold | conformal_calibrate | VAI_freeze
    reason: str
    mu: float
    s: float
    V_debiased: float
    V_naive: float
    bias_ratio: float
    y_flip: float
    J: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# MVP thresholds (calibrate on a benign probe in production)
THR = {
    "mu_mild": 0.05,
    "mu_high": 0.55,
    "s_spike": 0.05,
    "V_rel": 1.50,  # V_live > V_rel * V_ref (relative only)
    "y_flip": 0.02,
    "J_mal": 0.65,
    "bias_warn": 0.40,  # honesty: leaf bias is material
}


def triage(
    *,
    mu: float,
    s: float,
    V_debiased: float,
    V_naive: float,
    J: float = 0.0,
    y_flip: float = 0.0,
    V_ref: Optional[float] = None,
    thr: Optional[dict] = None,
) -> TriageDecision:
    """Straightforward next-step gate. One call → one actuator.

    Var(u) never freezes alone: needs Mag×Var spike and/or Y-flip.
    """
    t = dict(THR)
    if thr:
        t.update(thr)

    bias_ratio = float(
        max(float(V_naive) - float(V_debiased), 0.0) / max(float(V_naive), 1e-12)
    )

    s_spiked = s >= t["s_spike"]
    y_break = y_flip >= t["y_flip"]
    j_mal = J >= t["J_mal"]
    mu_high = mu >= t["mu_high"]
    mu_mild = mu >= t["mu_mild"]
    # Mag×Var malignant = co-move / J / high μ — NOT s alone
    # (synonym boards can wiggle s without integrity break)
    shape_mal = j_mal or mu_high or (s_spiked and mu_mild)

    # relative V spike vs REF — confirmation only
    V_rel_up = False
    if V_ref is not None:
        V_rel_up = float(V_debiased) > float(t["V_rel"]) * max(float(V_ref), 1e-8) + 1e-6

    if y_break or shape_mal or (V_rel_up and (y_break or shape_mal or mu_mild)):
        # V_rel_up alone is NOT enough; require shape/outcome, or mild μ with V↑
        # Re-check: freeze if integrity or Mag×Var malignant; V confirms.
        if y_break or shape_mal:
            call, actuator = "malignant", "VAI_freeze"
            bits = []
            if y_break:
                bits.append("Y-flip")
            if s_spiked:
                bits.append("s↑")
            if j_mal:
                bits.append("J↑")
            if mu_high:
                bits.append("μ↑↑")
            if V_rel_up:
                bits.append("V↑_rel")
            reason = "malignant: " + "+".join(bits)
        elif V_rel_up and mu_mild:
            # mild mean lift + relative V↑ without Y-flip: still conformal (benign)
            call, actuator = "benign", "conformal_calibrate"
            reason = "benign: mild μ + V↑_rel but no Y-flip / Mag×Var spike"
        else:
            call, actuator = "quiet", "hold"
            reason = "quiet: no integrity / Mag×Var spike"
    elif mu_mild:
        call, actuator = "benign", "conformal_calibrate"
        reason = "benign: mild μ, no Y-flip, Mag×Var not spiked"
    else:
        call, actuator = "quiet", "hold"
        reason = "quiet: μ and Mag×Var near REF"

    if bias_ratio >= t["bias_warn"] and call != "malignant":
        reason += f"; honesty: leaf bias_ratio={bias_ratio:.2f} (debias mattered)"

    return TriageDecision(
        call=call,
        actuator=actuator,
        reason=reason,
        mu=float(mu),
        s=float(s),
        V_debiased=float(V_debiased),
        V_naive=float(V_naive),
        bias_ratio=float(bias_ratio),
        y_flip=float(y_flip),
        J=float(J),
    )


def estimator_priority() -> list:
    """Product default order — do not invent a fourth route."""
    return [
        {
            "name": "causal_forest_leaf_debiased",
            "role": "default",
            "when": "LIVE monitoring; leaf merge min n=30 (50 if n≥200)",
        },
        {
            "name": "bcf_posterior_rpy2",
            "role": "bayesian",
            "when": "use rpy2 + R package bcf when Bayesian path is required",
        },
        {
            "name": "bootstrap",
            "role": "fallback",
            "when": "intuitive sanity check; small B + Winsorize u",
        },
    ]
