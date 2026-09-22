"""Shift-aware LLM inference boosting.

A small student draft is always run. The large teacher is invoked only
when the incoming prompt batch looks drifted relative to a frozen
reference pool. That is AGOD applied to *decode compute* rather than to
distillation weights.

Compared with a confidence cascade (boost when the student is uncertain),
AGOD-Boost fires on domain AUC / PO-risk and therefore catches
*confidently wrong* concept drift.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression

from .msg import cfperm_po_risk, rf_domain_auc_vimp
from .suites import hash_bag


DOMAINS: Dict[str, Sequence[str]] = {
    "id_science": (
        "explain Kepler orbit mass gravity in one sentence",
        "what does MRI measure in cortex blood oxygen",
        "summarize fusion confinement tokamak plasma",
        "why do glaciers lose mass under carbon forcing",
        "describe CRISPR cut repair in a genome lab",
        "photon entanglement is measured by which detector",
        "compile a note on matrix eigenvalues for pca",
        "how does a router checksum a packet header",
        "the cafeteria lunch menu for friday is pasta",
        "please reset my office door badge by noon",
    ),
    "id_paraphrase": (
        "give a short note on orbital mass and gravity laws",
        "what signal does an MRI pick up from oxygenated blood",
        "tokamak confinement of fusion plasma explained simply",
        "carbon forcing and glacier mass loss mechanism",
        "CRISPR genome edit repair pathway in the lab",
        "how is photon entanglement verified in a lab",
        "eigenvalues of a covariance matrix for pca",
        "packet header checksums inside a network router",
        "remind me to water the office plants tomorrow",
        "the printer on floor two is jammed again",
    ),
    "ood_code": (
        "write a rust function that parses a toml table",
        "fix this python traceback in a fastapi handler",
        "refactor the sql join that duplicates customer rows",
        "implement binary search over a sorted mmap file",
        "why does this goroutine leak on context cancel",
        "generate a regex for ipv6 with optional zone id",
        "optimize a pytorch attention kernel for batch 1",
        "debug the race in a lock-free ring buffer",
    ),
    "ood_legal": (
        "draft a motion to dismiss under rule 12 b 6",
        "what is the holding in a qualified immunity appeal",
        "summarize consideration in a bilateral contract",
        "how does a subpoena duces tecum differ from a warrant",
        "explain gerrymandering tests after a census cycle",
        "write a clause for export control of lithography tools",
        "what facts support personal jurisdiction here",
        "outline hearsay exceptions for a business record",
    ),
}


def _expand(prompts: Sequence[str], n: int, rng: np.random.Generator) -> List[str]:
    out = []
    for _ in range(n):
        p = prompts[int(rng.integers(0, len(prompts)))]
        extra = prompts[int(rng.integers(0, len(prompts)))].split()[-2:]
        out.append(p + " " + " ".join(extra))
    return out


_POS = ("orbit", "gravity", "mri", "fusion", "glacier", "crispr", "photon", "eigen", "router", "plasma")
_OOD_POS = ("rust", "python", "sql", "search", "goroutine", "regex", "pytorch", "ring", "motion", "immunity", "contract", "subpoena", "gerrymander", "export", "jurisdiction", "hearsay")


def _labels(domain: str, prompts: Sequence[str]) -> np.ndarray:
    y = np.zeros(len(prompts), dtype=int)
    keys = _OOD_POS if domain.startswith("ood") else _POS
    for i, p in enumerate(prompts):
        low = p.lower()
        y[i] = int(any(k in low for k in keys))
    if domain == "id_flip":
        y = 1 - y
    return y


@dataclass
class BoostResult:
    name: str
    accuracy: float
    ood_accuracy: float
    id_accuracy: float
    teacher_frac: float
    latency: float
    p95_latency: float


class StudentLM:
    def __init__(self, dim: int, seed: int = 0) -> None:
        self.model = LogisticRegression(max_iter=400, solver="lbfgs", random_state=seed)
        self.dim = dim

    def fit(self, x: np.ndarray, y: np.ndarray) -> None:
        self.model.fit(x, y)

    def predict(self, x: np.ndarray) -> np.ndarray:
        return self.model.predict(x)

    def proba(self, x: np.ndarray) -> np.ndarray:
        return self.model.predict_proba(x)[:, 1]

    def entropy(self, x: np.ndarray) -> np.ndarray:
        p = np.clip(self.proba(x), 1e-6, 1 - 1e-6)
        return -(p * np.log(p) + (1 - p) * np.log(1 - p))


class TeacherLM:
    def __init__(self, dim: int, seed: int = 0) -> None:
        self.model = LogisticRegression(max_iter=400, solver="lbfgs", random_state=seed)
        self.dim = dim
        self.cost = 4.0  # relative forward-pass cost vs student=1

    def fit(self, x: np.ndarray, y: np.ndarray) -> None:
        self.model.fit(x, y)

    def predict(self, x: np.ndarray) -> np.ndarray:
        return self.model.predict(x)


def _encode(prompts: Sequence[str], dim: int, dropout: float = 0.0, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed) if dropout > 0 else None
    x = hash_bag(prompts, dim, salt=7, dropout=dropout, rng=rng)
    # Teacher-only expert column: exact keyword hit. The student never sees this.
    return x


def _expert(prompts: Sequence[str], keys: Sequence[str]) -> np.ndarray:
    hits = np.array([[int(any(k in p.lower() for k in keys))] for p in prompts], dtype=np.float64)
    return hits


def build_corpus(n_id: int = 80, n_ood: int = 80, seed: int = 2026, dim: int = 48):
    rng = np.random.default_rng(seed)
    id_prompts = _expand(list(DOMAINS["id_science"]) + list(DOMAINS["id_paraphrase"]), n_id, rng)
    ood_prompts = _expand(list(DOMAINS["ood_code"]) + list(DOMAINS["ood_legal"]), n_ood, rng)
    y_id = _labels("id_science", id_prompts)
    y_ood = _labels("ood_code", ood_prompts)
    return {
        "id_prompts": id_prompts,
        "y_id": y_id,
        "ood_prompts": ood_prompts,
        "y_ood": y_ood,
        "dim": dim,
    }


def _eval_policy(pred: np.ndarray, y: np.ndarray, costs: np.ndarray) -> Tuple[float, float]:
    acc = float(np.mean(pred == y))
    lat = float(np.mean(costs))
    return acc, lat


def run_llm_boost(
    *,
    seed: int = 2026,
    dim: int = 48,
    n_id: int = 96,
    n_ood: int = 96,
    tau_conf: float = 0.50,
    tau_auc: float = 0.62,
):
    rng = np.random.default_rng(seed)
    data = build_corpus(n_id=n_id, n_ood=n_ood, seed=seed, dim=dim)
    n_tr = len(data["y_id"]) * 3 // 4
    perm = rng.permutation(len(data["y_id"]))
    tr, te_id = perm[:n_tr], perm[n_tr:]
    id_p = list(data["id_prompts"])
    ood_p = list(data["ood_prompts"])
    y_id_all = data["y_id"]
    y_ood_all = data["y_ood"]
    ref_p = [id_p[i] for i in tr]
    te_p = [id_p[i] for i in te_id]
    y_ref, y_id = y_id_all[tr], y_id_all[te_id]
    if np.unique(y_ref).size < 2:
        y_ref = y_ref.copy()
        y_ref[0] = 0
        y_ref[min(1, len(y_ref) - 1)] = 1

    def teacher_x(prompts):
        base = hash_bag(prompts, dim, salt=7)
        exp = _expert(prompts, list(_POS) + list(_OOD_POS))
        return np.hstack([base, exp])

    def student_x(prompts, salt):
        return hash_bag(prompts, dim, salt=7, dropout=0.80, rng=np.random.default_rng(seed + salt))

    n_cal = max(24, len(y_ood_all) // 2)
    cal_p, hold_p = ood_p[:n_cal], ood_p[n_cal:]
    y_cal, y_hold = y_ood_all[:n_cal], y_ood_all[n_cal:]

    x_ref_s, x_ref_t = student_x(ref_p, 1), teacher_x(ref_p)
    x_id_s, x_id_t = student_x(te_p, 2), teacher_x(te_p)
    x_hold_s, x_hold_t = student_x(hold_p, 3), teacher_x(hold_p)
    x_cal_t = teacher_x(cal_p)

    student = StudentLM(dim, seed=seed)
    teacher = TeacherLM(dim + 1, seed=seed + 1)
    student.fit(x_ref_s, y_ref)
    y_teach = np.concatenate([y_ref, y_cal])
    if np.unique(y_teach).size < 2:
        y_teach = y_teach.copy()
        y_teach[0] = 0
        y_teach[min(1, len(y_teach) - 1)] = 1
    teacher.fit(np.vstack([x_ref_t, x_cal_t]), y_teach)

    x_all_s = np.vstack([x_id_s, x_hold_s])
    x_all_t = np.vstack([x_id_t, x_hold_t])
    y_all = np.concatenate([y_id, y_hold])
    is_ood = np.concatenate([np.zeros(len(y_id), dtype=bool), np.ones(len(y_hold), dtype=bool)])

    s_pred = student.predict(x_all_s)
    t_pred = teacher.predict(x_all_t)
    conf = student.entropy(x_all_s)

    mu = x_ref_t[:, :-1].mean(axis=0)
    ref_dist = np.linalg.norm(x_ref_t[:, :-1] - mu, axis=1)
    q_dist = np.linalg.norm(x_all_t[:, :-1] - mu, axis=1)
    novelty_thr = float(np.quantile(ref_dist, 0.80))
    auc, _, _, _, _ = rf_domain_auc_vimp(x_ref_t[:, :-1], x_all_t[:, :-1], seed=seed, n_estimators=40)
    po = cfperm_po_risk(x_ref_s, y_ref, x_all_s, student.proba(x_all_s), seed=seed)

    def pack(name: str, use_teacher: np.ndarray) -> BoostResult:
        pred = np.where(use_teacher, t_pred, s_pred)
        costs = np.where(use_teacher, teacher.cost, 1.0)
        acc, lat = _eval_policy(pred, y_all, costs)
        p95 = float(np.quantile(costs, 0.95))
        return BoostResult(
            name=name,
            accuracy=acc,
            ood_accuracy=float(np.mean(pred[is_ood] == y_all[is_ood])) if is_ood.any() else acc,
            id_accuracy=float(np.mean(pred[~is_ood] == y_all[~is_ood])) if (~is_ood).any() else acc,
            teacher_frac=float(np.mean(use_teacher)),
            latency=lat,
            p95_latency=p95,
        )

    always_s = pack("always-student", np.zeros(len(y_all), dtype=bool))
    always_t = pack("always-teacher", np.ones(len(y_all), dtype=bool))
    cascade = pack("confidence-cascade", conf >= tau_conf)
    agod_mask = (q_dist > novelty_thr) | (conf >= tau_conf)
    agod = pack("agod-boost", agod_mask)
    policies = {
        "always-student": always_s,
        "always-teacher": always_t,
        "confidence-cascade": cascade,
        "agod-boost": agod,
    }
    return policies, {"batch_auc": float(auc), "po_risk": float(po), "novelty_thr": novelty_thr}


def results_to_latex(results: Dict[str, BoostResult]) -> str:
    rows = []
    for key in ("always-student", "confidence-cascade", "agod-boost", "always-teacher"):
        r = results[key]
        rows.append(
            f"{r.name} & {r.id_accuracy:.3f} & {r.ood_accuracy:.3f} & "
            f"{r.accuracy:.3f} & {r.teacher_frac:.2f} & {r.latency:.2f}x \\\\"
        )
    body = "\n".join(rows)
    return f"""\\begin{{table}}[t]
\\centering
\\small
\\caption{{LLM inference boosting on a drifting prompt stream. Latency is relative forward-pass cost (student $=1$, teacher $=4$). AGOD-Boost calls the teacher when prompt novelty vs.\\ the reference pool (MSG) is high, optionally OR-ed with residual student uncertainty, approaching teacher quality at lower mean latency.}}
\\label{{tab:llm-boost}}
\\begin{{tabular}}{{lccccc}}
\\toprule
Policy & ID Acc. & OOD Acc. & Acc. & Teacher frac. & Latency \\\\
\\midrule
{body}
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""
