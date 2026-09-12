"""AGOD: Attribution-Guided Online adaptation controllers.

Reusable control-plane pieces for multimodal online LR routing:

- ``mmd``: unbiased RBF MMD² (covariate / P(X) sensor)
- ``shift``: concept vs covariate decomposition
- ``lr_controller``: Softmax α → per-modality LR multipliers (actuator)

Amazon smoke CLI: ``python3 scripts/run_agod_amazon_mmd_lr.py``
"""

from .adapter import (
    ADAPTER_LAYERS,
    GATE_POLICIES,
    AdapterConfig,
    AdapterDecision,
    DynamicGate,
    LayeredAdapter,
    alpha_entropy,
    characterize_gate_traj,
    characterize_soft_lr_traj,
    flops_rel_proj,
    n_star_for_lift,
    soft_lr_dispersion,
)
from .lr_controller import (
    EMARouter,
    SCHEDULER_NAMES,
    aligned_cos_sim,
    alpha_to_lr,
    common_dim_grad_signature,
    cos_sim,
    damp_lr,
    equal_lr,
    flatten_grads,
    intensity_gain,
    modality_grad_cosine,
    schedule_modality_lr,
    soft_budget_lr,
    soft_cosine_lr,
    soft_entropy_lr,
    soft_gradcos_lr,
    soft_warmup_lr,
    softmax_scores,
    z_norm,
)
from .mmd import rbf_mmd2, whiten_pair
from .shift import decompose_hybrid, decompose_mmd, decompose_rf

__all__ = [
    "ADAPTER_LAYERS",
    "GATE_POLICIES",
    "AdapterConfig",
    "AdapterDecision",
    "DynamicGate",
    "LayeredAdapter",
    "characterize_gate_traj",
    "characterize_soft_lr_traj",
    "soft_lr_dispersion",
    "alpha_entropy",
    "n_star_for_lift",
    "flops_rel_proj",
    "rbf_mmd2",
    "whiten_pair",
    "decompose_mmd",
    "decompose_rf",
    "decompose_hybrid",
    "z_norm",
    "softmax_scores",
    "alpha_to_lr",
    "intensity_gain",
    "EMARouter",
    "SCHEDULER_NAMES",
    "schedule_modality_lr",
    "equal_lr",
    "soft_cosine_lr",
    "soft_budget_lr",
    "soft_warmup_lr",
    "damp_lr",
    "soft_entropy_lr",
    "soft_gradcos_lr",
    "cos_sim",
    "flatten_grads",
    "common_dim_grad_signature",
    "aligned_cos_sim",
    "modality_grad_cosine",
]

__version__ = "0.1.0"
