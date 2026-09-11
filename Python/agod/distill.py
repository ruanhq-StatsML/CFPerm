"""Frozen teacher + linear student with global KL matching and local alignment."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, MutableMapping, Sequence

import numpy as np

from .chronoberg import MODALITIES, TimeWindowBatch


def _l2_normalize(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    nrm = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.clip(nrm, eps, None)


def pairwise_softmax(z: np.ndarray, temperature: float = 0.07) -> np.ndarray:
    sim = (z @ z.T) / max(float(temperature), 1e-6)
    sim = sim - np.max(sim, axis=1, keepdims=True)
    e = np.exp(np.clip(sim, -40.0, 40.0))
    return e / np.clip(e.sum(axis=1, keepdims=True), 1e-8, None)


def kl_batch(p: np.ndarray, q: np.ndarray, eps: float = 1e-8) -> float:
    p = np.clip(p, eps, 1.0)
    q = np.clip(q, eps, 1.0)
    return float(np.mean(np.sum(p * (np.log(p) - np.log(q)), axis=1)))


@dataclass
class Teacher:
    """Frozen random projections standing in for CLIP / CLAP / EmbeddingGemma."""

    projections: Dict[str, np.ndarray]
    out_dim: int

    @classmethod
    def from_dims(
        cls,
        dims: Mapping[str, int],
        out_dim: int,
        rng: np.random.Generator,
        modalities: Sequence[str] = MODALITIES,
    ) -> "Teacher":
        proj = {}
        for m in modalities:
            w = rng.normal(size=(int(dims[m]), out_dim))
            w /= np.linalg.norm(w, axis=0, keepdims=True) + 1e-8
            proj[m] = w
        return cls(projections=proj, out_dim=out_dim)

    def embed_modality(self, modality: str, x: np.ndarray) -> np.ndarray:
        return _l2_normalize(np.asarray(x, dtype=np.float64) @ self.projections[modality])

    def embed_batch(self, batch: TimeWindowBatch) -> Dict[str, np.ndarray]:
        return {m: self.embed_modality(m, batch.X[m]) for m in batch.X}

    def fused(self, embeds: Mapping[str, np.ndarray], modalities: Sequence[str] = MODALITIES) -> np.ndarray:
        stacked = np.stack([embeds[m] for m in modalities], axis=0)
        return _l2_normalize(stacked.mean(axis=0))


@dataclass
class Student:
    """Rank-deficient projector per modality. Capacity is intentionally limited."""

    weights: Dict[str, np.ndarray]
    compress: Dict[str, np.ndarray]
    out_dim: int

    @classmethod
    def from_teacher(
        cls,
        teacher: Teacher,
        rng: np.random.Generator,
        scale: float = 0.08,
        rank_ratio: float = 0.40,
    ) -> "Student":
        compress: Dict[str, np.ndarray] = {}
        weights: Dict[str, np.ndarray] = {}
        for m, w_t in teacher.projections.items():
            d_in = int(w_t.shape[0])
            d_obs = max(3, int(round(d_in * rank_ratio)))
            c = rng.normal(size=(d_in, d_obs))
            c /= np.linalg.norm(c, axis=0, keepdims=True) + 1e-8
            compress[m] = c
            weights[m] = scale * rng.normal(size=(d_obs, teacher.out_dim)) / np.sqrt(d_obs)
        return cls(weights=weights, compress=compress, out_dim=teacher.out_dim)

    def clone(self) -> "Student":
        return Student(
            weights={m: w.copy() for m, w in self.weights.items()},
            compress={m: c.copy() for m, c in self.compress.items()},
            out_dim=self.out_dim,
        )

    def observe(self, modality: str, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, dtype=np.float64) @ self.compress[modality]

    def embed_modality(self, modality: str, x: np.ndarray) -> np.ndarray:
        return self.observe(modality, x) @ self.weights[modality]

    def embed_batch(self, batch: TimeWindowBatch) -> Dict[str, np.ndarray]:
        return {m: self.embed_modality(m, batch.X[m]) for m in self.weights}

    def fused(
        self,
        embeds: Mapping[str, np.ndarray],
        modalities: Sequence[str] = MODALITIES,
    ) -> np.ndarray:
        stacked = np.stack([embeds[m] for m in modalities], axis=0)
        return _l2_normalize(stacked.mean(axis=0))

    def parameters(self) -> MutableMapping[str, np.ndarray]:
        return self.weights


@dataclass
class LossBreakdown:
    total: float
    global_kl: float
    local: Dict[str, float]
    omega: float


def agod_loss(
    student_mod: Mapping[str, np.ndarray],
    teacher_mod: Mapping[str, np.ndarray],
    alpha: Mapping[str, float],
    *,
    omega: float = 0.0,
    lam: float = 0.0,
    kl_temperature: float = 0.07,
    modalities: Sequence[str] = MODALITIES,
) -> LossBreakdown:
    s_fused = _l2_normalize(np.stack([student_mod[m] for m in modalities], axis=0).mean(axis=0))
    t_fused = _l2_normalize(np.stack([teacher_mod[m] for m in modalities], axis=0).mean(axis=0))
    global_kl = kl_batch(pairwise_softmax(t_fused, kl_temperature), pairwise_softmax(s_fused, kl_temperature))
    local = {}
    weighted = 0.0
    for m in modalities:
        s = _l2_normalize(student_mod[m])
        t = _l2_normalize(teacher_mod[m])
        local[m] = float(np.mean(np.sum((s - t) ** 2, axis=1)))
        weighted += float(alpha[m]) * local[m]
    total = global_kl + weighted + float(lam) * float(omega)
    return LossBreakdown(total=total, global_kl=global_kl, local=local, omega=float(omega))


def _grad_mse_normalized(x: np.ndarray, s_raw: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Gradient of mean ||s_hat - t||^2 w.r.t. projector W, s_raw = X @ W, s_hat = s_raw / ||.||."""
    s = _l2_normalize(s_raw)
    t = _l2_normalize(t)
    n, d = s.shape
    # d||s-t||^2 / ds_hat = 2(s-t); chain through L2 norm.
    diff = 2.0 * (s - t) / n
    dots = np.sum(diff * s, axis=1, keepdims=True)
    nrm = np.clip(np.linalg.norm(s_raw, axis=1, keepdims=True), 1e-8, None)
    d_raw = (diff - s * dots) / nrm
    return x.T @ d_raw


def _grad_kl_fused(
    x_mods: Mapping[str, np.ndarray],
    s_raw: Mapping[str, np.ndarray],
    t_fused: np.ndarray,
    *,
    temperature: float,
    modalities: Sequence[str],
) -> Dict[str, np.ndarray]:
    """Finite-difference-free gradient of KL(teacher_sim || student_sim) through mean fusion."""
    s_stack = np.stack([s_raw[m] for m in modalities], axis=0)
    fused_raw = s_stack.mean(axis=0)
    s_fused = _l2_normalize(fused_raw)
    t_fused = _l2_normalize(t_fused)
    p = pairwise_softmax(t_fused, temperature)
    q = pairwise_softmax(s_fused, temperature)
    # dKL/dsim_s via softmax Jacobian: (q - p) / tau, using teacher as target.
    n = s_fused.shape[0]
    grad_sim = (q - p) / max(temperature, 1e-6) / n
    grad_fused = (grad_sim + grad_sim.T) @ s_fused
    nrm = np.clip(np.linalg.norm(fused_raw, axis=1, keepdims=True), 1e-8, None)
    dots = np.sum(grad_fused * s_fused, axis=1, keepdims=True)
    d_raw_fused = (grad_fused - s_fused * dots) / nrm
    d_raw_each = d_raw_fused / len(modalities)
    grads = {}
    for m in modalities:
        grads[m] = x_mods[m].T @ d_raw_each
    return grads


def sgd_step(
    student: Student,
    batch: TimeWindowBatch,
    teacher_mod: Mapping[str, np.ndarray],
    alpha: Mapping[str, float],
    *,
    lr: float = 0.08,
    kl_temperature: float = 0.07,
    modalities: Sequence[str] = MODALITIES,
) -> None:
    x_mods = {m: np.asarray(batch.X[m], dtype=np.float64) for m in modalities}
    obs = {m: student.observe(m, x_mods[m]) for m in modalities}
    s_raw = {m: obs[m] @ student.weights[m] for m in modalities}
    t_fused = _l2_normalize(np.stack([teacher_mod[m] for m in modalities], axis=0).mean(axis=0))
    g_kl = _grad_kl_fused(obs, s_raw, t_fused, temperature=kl_temperature, modalities=modalities)
    for m in modalities:
        g_local = _grad_mse_normalized(obs[m], s_raw[m], teacher_mod[m])
        grad = g_kl[m] + float(alpha[m]) * g_local
        student.weights[m] -= lr * grad
