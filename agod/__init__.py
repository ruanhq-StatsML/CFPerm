"""AGOD: Attribution-Guided Online adaptation controllers.

Reusable control-plane pieces for multimodal online LR routing:

- ``mmd``: unbiased RBF MMD² (covariate / P(X) sensor)
- ``shift``: concept vs covariate decomposition
- ``lr_controller``: Softmax α → per-modality LR multipliers (actuator)

Amazon smoke CLI: ``python3 scripts/run_agod_amazon_mmd_lr.py``
"""

from .lr_controller import EMARouter, alpha_to_lr, intensity_gain, softmax_scores, z_norm
from .mmd import rbf_mmd2, whiten_pair
from .shift import decompose_hybrid, decompose_mmd, decompose_rf

__all__ = [
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
]

__version__ = "0.1.0"
