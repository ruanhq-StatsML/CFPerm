"""Online AGOD loop: monitor → route → distill → gate → fallback."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .chronoberg import MODALITIES, ChronoBergStream, TimeWindowBatch
from .distill import Student, Teacher, agod_loss, sgd_step
from .metrics import cosine_alignment, recall_at_k
from .msg import compute_msg_state
from .routing import route_from_state


@dataclass
class AGODConfig:
    tau: float = 0.35
    gamma: float = 1.35
    lam: float = 0.25
    momentum: float = 0.35
    theta_low: float = 0.55
    theta_high: float = 0.80
    decay: float = 0.45
    lr: float = 0.12
    steps_per_window: int = 8
    pretrain_steps: int = 55
    teacher_dim: int = 12
    kl_temperature: float = 0.07
    seed: int = 2026


@dataclass
class StepLog:
    year: int
    t_index: int
    baseline: str
    alpha: Dict[str, float]
    auc: Dict[str, float]
    po_risk: Dict[str, float]
    gap: Dict[str, float]
    localize: Sequence[str]
    omega: float
    loss: float
    recall: float
    drift_recall: float
    alignment: Dict[str, float]
    fallback: Dict[str, bool]


@dataclass
class AGODTrainer:
    teacher: Teacher
    student: Student
    config: AGODConfig
    baseline: str = "B3"
    modalities: Sequence[str] = MODALITIES
    alpha_prev: Optional[np.ndarray] = None
    logs: List[StepLog] = field(default_factory=list)

    def pretrain(self, reference: TimeWindowBatch) -> None:
        n = len(self.modalities)
        alpha = {m: 1.0 / n for m in self.modalities}
        t_mod = self.teacher.embed_batch(reference)
        for _ in range(self.config.pretrain_steps):
            sgd_step(
                self.student,
                reference,
                t_mod,
                alpha,
                lr=self.config.lr,
                kl_temperature=self.config.kl_temperature,
                modalities=self.modalities,
            )

    def step(
        self,
        reference: TimeWindowBatch,
        current: TimeWindowBatch,
        eval_batch: Optional[TimeWindowBatch] = None,
    ) -> StepLog:
        state = compute_msg_state(
            reference,
            current,
            gamma=self.config.gamma,
            seed=self.config.seed + 1000 * current.t_index,
        )
        alpha_map, localize, omega = route_from_state(
            state,
            tau=self.config.tau,
            theta_low=self.config.theta_low,
            theta_high=self.config.theta_high,
            decay=self.config.decay,
            momentum=self.config.momentum,
            alpha_prev=self.alpha_prev,
            modalities=self.modalities,
            baseline=self.baseline,
        )
        t_mod = self.teacher.embed_batch(current)
        for _ in range(self.config.steps_per_window):
            sgd_step(
                self.student,
                current,
                t_mod,
                alpha_map,
                lr=self.config.lr,
                kl_temperature=self.config.kl_temperature,
                modalities=self.modalities,
            )
        eval_on = eval_batch or current
        s_mod = self.student.embed_batch(eval_on)
        t_mod = self.teacher.embed_batch(eval_on)
        breakdown = agod_loss(
            self.student.embed_batch(current),
            self.teacher.embed_batch(current),
            alpha_map,
            omega=omega,
            lam=self.config.lam,
            kl_temperature=self.config.kl_temperature,
            modalities=self.modalities,
        )
        s_fused = self.student.fused(s_mod, self.modalities)
        t_fused = self.teacher.fused(t_mod, self.modalities)
        rec = recall_at_k(s_fused, t_fused, k=5)
        if eval_on.drifted_audio:
            drift_recall = recall_at_k(s_mod["audio"], t_mod["audio"], k=5)
        else:
            drift_recall = rec
        align = {
            m: cosine_alignment(s_mod[m], t_mod[m]) for m in self.modalities
        }
        self.alpha_prev = np.array([alpha_map[m] for m in self.modalities], dtype=np.float64)
        log = StepLog(
            year=current.year,
            t_index=current.t_index,
            baseline=self.baseline,
            alpha=alpha_map,
            auc={m: state.details[m].auc for m in self.modalities},
            po_risk={m: state.details[m].po_risk for m in self.modalities},
            gap=state.gaps,
            localize=localize,
            omega=omega,
            loss=breakdown.total,
            recall=rec,
            drift_recall=drift_recall,
            alignment=align,
            fallback={m: state.details[m].used_fallback for m in self.modalities},
        )
        self.logs.append(log)
        return log


def run_stream(
    stream: ChronoBergStream,
    baseline: str,
    config: Optional[AGODConfig] = None,
) -> AGODTrainer:
    config = config or AGODConfig(seed=stream.config.seed)
    rng = np.random.default_rng(config.seed)
    teacher = Teacher.from_dims(stream.dims, config.teacher_dim, rng)
    student = Student.from_teacher(teacher, rng)
    trainer = AGODTrainer(
        teacher=teacher,
        student=student,
        config=config,
        baseline=baseline,
    )
    reference = stream.reference()
    trainer.pretrain(reference)
    for batch in stream.iter_online():
        eval_batch = stream.make_window(batch.t_index, split="eval")
        trainer.step(reference, batch, eval_batch=eval_batch)
    return trainer


Baseline = str
