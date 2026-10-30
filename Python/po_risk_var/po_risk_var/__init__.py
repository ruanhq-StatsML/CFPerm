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
    leaf_assignment_forest,
    leaf_assignment_one_tree,
    merge_leaves_forest,
    merge_leaves_one_tree,
    monitor_batch,
    var_u_from_leaf_merged,
)
from po_risk_var.product import TriageDecision, estimator_priority, triage
from po_risk_var.suite import run_justify_suite

__version__ = "0.1.5"
__all__ = [
    "estimate_var_u",
    "monitor_batch",
    "leaf_assignment_one_tree",
    "leaf_assignment_forest",
    "merge_leaves_one_tree",
    "merge_leaves_forest",
    "var_u_from_leaf_merged",
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
