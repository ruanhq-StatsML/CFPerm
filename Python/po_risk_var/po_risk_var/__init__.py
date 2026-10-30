"""po-risk-var — PO-risk Var(u) product package.

Default estimator: Causal Forest leaf-size-debiased + leaf-wise merge.
Next-step gate: Mag×Var + Var(u) + optional Y-flip → typed actuator.
"""
from po_risk_var.apps_product import run_product_apps_suite
from po_risk_var.biz_value import (
    estimator_incremental_value,
    product_apps_incremental_value,
    run_biz_value_suite,
)
from po_risk_var.compare_methods import run_compare_suite
from po_risk_var.estimators import (
    estimate_var_u,
    merge_leaf_u,
    monitor_batch,
    tree_parents,
)
from po_risk_var.product import TriageDecision, estimator_priority, triage
from po_risk_var.suite import run_justify_suite

__version__ = "0.1.6"
__all__ = [
    "estimate_var_u",
    "monitor_batch",
    "tree_parents",
    "merge_leaf_u",
    "triage",
    "TriageDecision",
    "estimator_priority",
    "estimator_incremental_value",
    "product_apps_incremental_value",
    "run_justify_suite",
    "run_product_apps_suite",
    "run_biz_value_suite",
    "run_compare_suite",
    "__version__",
]
