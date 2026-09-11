"""Attribution-Guided Online Distillation (AGOD).

Reframes CFPerm / RF-Domain modality-specific gap decomposition as a
clever covariate that routes distillation capacity toward drifting
modalities. The update signal is a distribution-shift statistic, not
prediction error alone.
"""

from .chronoberg import (
    CHRONOBERG_WINDOWS,
    MODALITIES,
    ChronoBergStream,
    SyntheticChronoBergConfig,
    TimeWindowBatch,
)
from .distill import Student, Teacher, agod_loss
from .experiment import ExperimentResult, run_baseline_comparison
from .metrics import attribution_consistency, forgetting_rate, recall_at_k
from .msg import MSGResult, MSGState, compute_msg_state
from .online import AGODConfig, AGODTrainer, Baseline
from .routing import gate_weights, softmax_weights, smoothness_penalty

__all__ = [
    "AGODConfig",
    "AGODTrainer",
    "Baseline",
    "CHRONOBERG_WINDOWS",
    "ChronoBergStream",
    "ExperimentResult",
    "MODALITIES",
    "MSGResult",
    "MSGState",
    "Student",
    "SyntheticChronoBergConfig",
    "Teacher",
    "TimeWindowBatch",
    "agod_loss",
    "attribution_consistency",
    "compute_msg_state",
    "forgetting_rate",
    "gate_weights",
    "recall_at_k",
    "run_baseline_comparison",
    "smoothness_penalty",
    "softmax_weights",
]

__version__ = "0.1.0"
