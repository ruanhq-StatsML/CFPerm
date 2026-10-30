"""po-risk-var — PO-risk Var(u) product package.

Default estimator: Causal Forest leaf-size-debiased + leaf-wise merge.
Next-step gate: Mag×Var + Var(u) + optional Y-flip → typed actuator.
"""
from po_risk_var.apps_product import run_product_apps_suite
from po_risk_var.estimators import estimate_var_u, monitor_batch
from po_risk_var.product import TriageDecision, estimator_priority, triage
from po_risk_var.suite import run_justify_suite

__version__ = "0.1.1"
__all__ = [
    "estimate_var_u",
    "monitor_batch",
    "triage",
    "TriageDecision",
    "estimator_priority",
    "run_justify_suite",
    "run_product_apps_suite",
    "__version__",
]
