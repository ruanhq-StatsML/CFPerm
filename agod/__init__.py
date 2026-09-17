"""Streaming gates used by the manuscript LLM-audit prototype.

Two related but distinct monitors:

- ``OnlineARBootstrap`` (Palm & Nagler): frozen-ref excess
  ``Δ_t = s_t − μ_ref``, fire when the online AR-bootstrap CI is
  entirely above 0.
- last-two ``hop_fires`` (OnlineRFPerm): consecutive OOS ratio
  ``e_now / e_prev ≥ γ`` with an error floor.
"""

from agod.online_ar_bootstrap import BETA, OnlineARBootstrap, run_delta_bootstrap

__all__ = ["BETA", "OnlineARBootstrap", "run_delta_bootstrap"]
