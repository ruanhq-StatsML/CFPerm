"""AGOD lightweight package init (PO-IPTW + rolling PO-learner refit)."""

from agod.online_rfperm import (
    fit_online_probe,
    hop_fires,
    run_rfperm_stream,
    score_probe,
)
from agod.po_iptw import instance_po_risk, po_iptw_weights
from agod.po_refit import (
    refit_po_weights,
    run_adaptive_stream,
    run_oracle_switch,
    run_resid_stream,
    run_refit_stream,
    run_switch_stream,
)

__all__ = [
    "fit_online_probe",
    "hop_fires",
    "instance_po_risk",
    "po_iptw_weights",
    "refit_po_weights",
    "run_adaptive_stream",
    "run_oracle_switch",
    "run_resid_stream",
    "run_refit_stream",
    "run_rfperm_stream",
    "run_switch_stream",
    "score_probe",
]
