"""Turn the clever covariate into the next planning prompt.

H is computed only from a finished outcome. The prompt is for the following
step. It does not add a UCB term.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from online_llm_stack import clever_covariate


class RAPCleverCovariateSkill:
    """rap_clever_covariate_guide, version 0.1.0.

    The planner sees H and the adjustment. It does not see a UCB bonus.
    """

    name = "rap_clever_covariate_guide"
    version = "0.1.0"

    def __init__(self, threshold: float = 0.5):
        self.threshold = float(threshold)

    def invoke(
        self,
        state: str,
        prior: float,
        outcome: float,
        history: Optional[Sequence[str]] = None,
    ) -> dict:
        scored = clever_covariate(float(outcome), float(prior))
        # The shared formula uses 0.5. A caller threshold only renames the hint.
        h = float(scored["H"])
        if h > self.threshold:
            hint = "aggressive"
        elif h < -self.threshold:
            hint = "conservative"
        else:
            hint = "keep"
        history_text = "\n".join(history) if history else "None"
        e = float(scored["prior"])
        prompt = (
            "[Clever Covariate Signal]\n"
            f"Previous prior: {e:.3f}\n"
            f"Observed outcome: {float(outcome):.0f}\n"
            f"Clever covariate H: {h:.3f}\n"
            f"Online anomaly |H|: {abs(h):.3f}\n"
            f"Adjustment hint: {hint}\n"
            "\n"
            "- If hint = aggressive: the outcome was much better than expected.\n"
            "  Explore more boldly, consider higher-risk/higher-reward branches.\n"
            "- If hint = conservative: the outcome was much worse than expected.\n"
            "  Be more cautious, prefer safe and proven branches, reduce exploration.\n"
            "- If hint = keep: the outcome matched expectation.\n"
            "  Continue with the current strategy.\n"
            "\n"
            "[Current State]\n"
            f"{state}\n"
            "\n"
            "[Recent History]\n"
            f"{history_text}\n"
            "\n"
            "Based on the above, generate the next planning step or candidate actions.\n"
        )
        return {
            "prompt": prompt,
            "H": h,
            "anomaly": float(abs(h)),
            "hint": hint,
            "prior": e,
            "outcome": float(outcome),
        }
