"""Gap-guided LLM inference: two-sample localization as a reading budget.

Not TMLE. Not causal. W is a batch/domain indicator for a two-sample problem
on X. ê estimates P(W=1 | X). π and s_m(x) say which human block B_m localizes
that two-sample difference. They never identify an effect of W on Y.

The LLM is a reader with a token/tool budget. Three distinct maps:

  same expert     pick one m from π only (not from this x); every query is
                  read by the same specialist f_m(X_{B_m})
  context pack    serialize enabled blocks in read-order into one string;
                  one forward pass on concat_m∈E X_{B_m}  (not an opinion pool)
  instance gate   m(x) = argmax s_m(x)  — different experts for different x;
                  this is *not* “the same expert”

π is a shrinkage prior over channel identity C ∈ {1..M}. Abstain means the
posterior over C is too flat to justify *dropping* a block: pack all, or HITL.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Sequence

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

from clever_covariate_gap import (
    CLIP_E,
    EPS,
    ModalitySpec,
    _as_1d,
    _rf_clf,
    clever_z,
    decompose_modality_gap,
    fit_vimp,
    logit,
    relu_normalize,
    vimp_mass_share,
)

KNOWN_TOOLS: Dict[str, dict] = {
    "text": {
        "tool_id": "read_tokens",
        "description": "Read raw / hashed lexical tokens for this example.",
        "cost": 3.0,
        "action": "dump the text channel into the prompt window",
    },
    "valence": {
        "tool_id": "vad_valence",
        "description": "Look up valence lexicon statistics (mean/max/p90).",
        "cost": 1.0,
        "action": "call the valence lexicon tool before arguing from wording",
    },
    "arousal": {
        "tool_id": "vad_arousal",
        "description": "Look up arousal lexicon statistics.",
        "cost": 1.0,
        "action": "inspect arousal scores",
    },
    "dominance": {
        "tool_id": "vad_dominance",
        "description": "Look up dominance lexicon statistics.",
        "cost": 1.0,
        "action": "inspect dominance scores",
    },
    "premise": {
        "tool_id": "read_premise",
        "description": "Re-read the premise sentence only.",
        "cost": 1.5,
        "action": "load the premise into context",
    },
    "hypothesis": {
        "tool_id": "read_hypothesis",
        "description": "Re-read the hypothesis sentence only.",
        "cost": 1.5,
        "action": "load the hypothesis into context",
    },
    "overlap": {
        "tool_id": "lexical_overlap",
        "description": "Compute lexical / token overlap between the two sentences.",
        "cost": 1.0,
        "action": "run the overlap checker before NLI reasoning",
    },
    "hours": {
        "tool_id": "inspect_hours",
        "description": "Inspect hours-worked and related labor-supply fields.",
        "cost": 1.0,
        "action": "surface hours-worked before other tabular fields",
    },
    "capital": {
        "tool_id": "inspect_capital",
        "description": "Inspect capital-gain / capital-loss.",
        "cost": 1.0,
        "action": "surface capital fields",
    },
    "demography": {
        "tool_id": "inspect_demography",
        "description": "Inspect age / education / race block (not the batch key).",
        "cost": 1.5,
        "action": "load demography columns",
    },
    "work": {
        "tool_id": "inspect_work",
        "description": "Inspect occupation / workclass / hours-adjacent job fields.",
        "cost": 1.5,
        "action": "load work columns",
    },
    "style": {
        "tool_id": "style_markers",
        "description": "Inspect punctuation / length / surface-style markers.",
        "cost": 1.0,
        "action": "load style features",
    },
    "polarity": {
        "tool_id": "polarity_lexicon",
        "description": "Inspect polarity-lexicon atoms (with jitter; not a point mass).",
        "cost": 1.0,
        "action": "load polarity scores before marketplace-specific wording",
    },
}


@dataclass
class ToolChannel:
    name: str
    tool_id: str
    description: str
    cost: float = 1.0
    action: str = ""


@dataclass
class RouterConfig:
    """Deterministic gates. Tune on a held-out domain pair, then freeze."""

    lam: float = 0.4
    tau_hi: float = 0.42
    tau_lo: float = 0.18
    k_max: int = 2
    cost_budget: float = 4.0
    floor_enable: float = 0.08
    entropy_abstain: float | None = None
    clip: float = CLIP_E


@dataclass
class RoutingDecision:
    read_order: List[str]
    tools_enabled: List[str]
    tools_disabled: List[str]
    abstain: bool
    ask_missing: List[str]
    critic_must_cite: str | None
    blend: Dict[str, float]
    pi: Dict[str, float]
    instance_share: Dict[str, float]
    entropy: float
    rationale: str
    packed_blocks: List[str] = field(default_factory=list)
    pack_mode: str = "same_expert"
    packet: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = asdict(self)
        return d


@dataclass
class FrozenGapModels:
    """Deployed Stage-1 scorers. Fit once on the domain pair; score online."""

    spec: ModalitySpec
    names: List[str]
    pi: np.ndarray
    pi_vimp: np.ndarray
    prevalence: float
    full_clf: object
    block_clfs: List[object]
    lomo_clfs: List[object | None]
    has_lomo: bool
    clip: float = CLIP_E
    n_estimators: int = 80
    seed: int = 0

    def pack_pi(self) -> Dict[str, float]:
        return {n: round(float(v), 4) for n, v in zip(self.names, self.pi)}


@dataclass
class CriticVerdict:
    ok: bool
    cited: List[str]
    must_cite: str | None
    mismatch: bool
    reask_prompt: str


def simplex_entropy(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float).reshape(-1)
    p = np.clip(p, EPS, None)
    p = p / float(p.sum())
    return float(-np.sum(p * np.log(p)))


def entropy_gate(M: int, config: "RouterConfig") -> float:
    if config.entropy_abstain is not None:
        return float(config.entropy_abstain)
    return 0.92 * float(np.log(max(int(M), 2)))


def population_pack(
    names: Sequence[str],
    pi: np.ndarray,
    *,
    config: "RouterConfig | None" = None,
    k: int | None = None,
) -> tuple[List[str], str]:
    """Same template for every query: E = E(π), never E(x).

    pack_all     — π too flat; dropping a block is unjustified
    same_expert  — |E|=1, m̂ = argmax π  (the “same specialist”)
    concat_k     — |E|=k>1, still one reader on concatenated columns
    """
    config = config or RouterConfig()
    names = list(names)
    pi = relu_normalize(np.asarray(pi, dtype=float).reshape(-1))
    M = len(names)
    ent = simplex_entropy(pi)
    pmax = float(pi.max())
    order = [names[i] for i in np.argsort(-pi)]
    if ent >= entropy_gate(M, config) or pmax < config.tau_lo:
        return list(names), "pack_all"
    if pmax >= config.tau_hi:
        return [order[0]], "same_expert"
    kk = int(config.k_max if k is None else k)
    kk = max(1, min(kk, M))
    return order[:kk], "concat_k"


def packed_in_read_order(read_order: Sequence[str], enabled_names: Sequence[str]) -> List[str]:
    enabled = set(enabled_names)
    return [n for n in read_order if n in enabled]


def column_index(spec: ModalitySpec, enabled_names: Sequence[str]) -> np.ndarray:
    idx: List[int] = []
    for name in enabled_names:
        sl = spec.slices[spec.names.index(name)]
        idx.extend(range(sl.start, sl.stop))
    return np.asarray(idx, dtype=int)


def pack_blocks(X: np.ndarray, spec: ModalitySpec, enabled_names: Sequence[str]) -> np.ndarray:
    """Column-concat of enabled blocks in the given order (context-packing analogue)."""
    X = np.asarray(X, dtype=float)
    idx = column_index(spec, enabled_names)
    if idx.size == 0:
        return np.ones((X.shape[0], 1), dtype=float)
    return X[:, idx]


def pack_user_context(block_texts: Dict[str, str], decision: RoutingDecision) -> str:
    """Serialize only packed_blocks. Omitted channels are absent, not empty headers."""
    if decision.abstain or decision.pack_mode == "hitl_abstain":
        ask = ", ".join(decision.ask_missing) or "a more localized channel"
        return (
            "[no subset packed] The two-sample gap is not localized. "
            f"Do not concatenate a subset of blocks. Request: {ask}."
        )
    parts = []
    for name in decision.packed_blocks:
        if name not in block_texts:
            continue
        parts.append(f"### channel: {name}\n{block_texts[name]}")
    return "\n\n".join(parts)


def holdout_packed_auc(
    Xtr: np.ndarray,
    Wtr: np.ndarray,
    Xte: np.ndarray,
    Wte: np.ndarray,
    spec: ModalitySpec,
    enabled_names: Sequence[str],
    *,
    seed: int = 0,
    n_estimators: int = 50,
) -> tuple[float, int]:
    """One reader on the packed columns — not a mixture of block-wise scores."""
    Xptr = pack_blocks(Xtr, spec, enabled_names)
    Xpte = pack_blocks(Xte, spec, enabled_names)
    Wtr = _as_1d(Wtr).astype(int)
    clf = _rf_clf(Xptr.shape[1], len(Wtr), seed, n_estimators=n_estimators)
    clf.fit(Xptr, Wtr)
    pred = _predict_pos(clf, Xpte, float(Wtr.mean()), CLIP_E)
    return _auc(Wte, pred), int(Xptr.shape[1])


def make_diffuse_equal_shift(
    *,
    n: int = 500,
    d_text: int = 24,
    d_vad: int = 4,
    mean_shift: float = 0.95,
    seed: int = 0,
) -> tuple:
    """Each block gets the same one-coordinate shift — π should not concentrate."""
    rng = np.random.default_rng(seed)
    W = rng.integers(0, 2, size=n)
    names = ["text", "valence", "arousal", "dominance"]
    dims = [d_text, d_vad, d_vad, d_vad]
    blocks = []
    for d in dims:
        Xj = rng.normal(0.0, 1.0, size=(n, d))
        Xj[W == 1, 0] += mean_shift
        blocks.append(Xj)
    X = np.hstack(blocks)
    Y = rng.normal(0.0, 1.0, size=n)
    slices, start = [], 0
    for d in dims:
        slices.append(slice(start, start + d))
        start += d
    return X, W, Y, ModalitySpec(names=names, slices=slices)


def blend_shares(pi: np.ndarray, instance_share: np.ndarray, *, lam: float = 0.4) -> np.ndarray:
    """r = λπ + (1−λ)s, then renormalize. λ=1 is π-only (same expert for all x)."""
    pi = np.asarray(pi, dtype=float).reshape(1, -1)
    s = np.asarray(instance_share, dtype=float)
    if s.ndim == 1:
        s = s.reshape(1, -1)
    lam = float(np.clip(lam, 0.0, 1.0))
    r = lam * pi + (1.0 - lam) * s
    return r / (r.sum(axis=1, keepdims=True) + EPS)


def catalog_for_spec(spec: ModalitySpec) -> Dict[str, ToolChannel]:
    out: Dict[str, ToolChannel] = {}
    for name in spec.names:
        meta = KNOWN_TOOLS.get(
            name,
            {
                "tool_id": f"inspect_{name}",
                "description": f"Inspect the '{name}' channel for this example.",
                "cost": 1.0,
                "action": f"load the {name} block",
            },
        )
        out[name] = ToolChannel(name=name, **meta)
    return out


def _predict_pos(clf, X: np.ndarray, prevalence: float, clip: float) -> np.ndarray:
    n = int(np.asarray(X).shape[0])
    if clf is None:
        return np.full(n, float(np.clip(prevalence, clip, 1.0 - clip)))
    proba = clf.predict_proba(X)
    classes = list(clf.classes_)
    if 1 not in classes:
        return np.full(n, float(np.clip(prevalence, clip, 1.0 - clip)))
    return np.clip(proba[:, classes.index(1)].astype(float), clip, 1.0 - clip)


def fit_frozen_gap(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    *,
    seed: int = 0,
    n_estimators: int = 80,
    light: bool = False,
    pi: np.ndarray | None = None,
) -> FrozenGapModels:
    """Fit deployable ê / ê_m / ê_{-m}. π comes from OOF decompose on the same train fold."""
    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    n, d = X.shape
    names = list(spec.names)
    prevalence = float(W.mean())
    clip = CLIP_E

    if pi is None:
        gap = decompose_modality_gap(
            X, W, spec, seed=seed, n_splits=3, n_estimators=n_estimators, light=light,
        )
        pi = np.asarray(gap.pi_consensus, dtype=float)
        pi_vimp = np.asarray(gap.pi_vimp, dtype=float)
    else:
        pi = relu_normalize(np.asarray(pi, dtype=float).reshape(-1))
        vimp = fit_vimp(X, W, seed=seed, n_estimators=n_estimators)
        share = vimp_mass_share(vimp, spec)
        pi_vimp = np.array([share[n] for n in names], dtype=float)

    full = _rf_clf(d, n, seed, n_estimators=n_estimators)
    full.fit(X, W)

    block_clfs: List[object] = []
    for j, sl in enumerate(spec.slices):
        Xj = X[:, sl]
        clf = _rf_clf(Xj.shape[1], n, seed + 10 * (j + 1), n_estimators=n_estimators)
        clf.fit(Xj, W)
        block_clfs.append(clf)

    lomo_clfs: List[object | None] = []
    has_lomo = not light
    if has_lomo:
        for j, sl in enumerate(spec.slices):
            keep = [i for i in range(d) if i not in range(sl.start, sl.stop)]
            if not keep:
                lomo_clfs.append(None)
                continue
            Xl = X[:, keep]
            clf = _rf_clf(Xl.shape[1], n, seed + 40 * (j + 1), n_estimators=n_estimators)
            clf.fit(Xl, W)
            lomo_clfs.append(clf)
    else:
        lomo_clfs = [None] * len(names)

    return FrozenGapModels(
        spec=spec,
        names=names,
        pi=pi,
        pi_vimp=pi_vimp,
        prevalence=prevalence,
        full_clf=full,
        block_clfs=block_clfs,
        lomo_clfs=lomo_clfs,
        has_lomo=has_lomo,
        clip=clip,
        n_estimators=n_estimators,
        seed=seed,
    )


def predict_gap_scores(models: FrozenGapModels, X: np.ndarray) -> dict:
    """Online scores for a batch of queries. No W required."""
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    M = len(models.names)
    clip = models.clip
    e_full = _predict_pos(models.full_clf, X, models.prevalence, clip)
    e_block = np.zeros((n, M), dtype=float)
    instance_gap = np.zeros((n, M), dtype=float)
    for j, sl in enumerate(models.spec.slices):
        e_block[:, j] = _predict_pos(models.block_clfs[j], X[:, sl], models.prevalence, clip)
        if models.has_lomo:
            keep = [i for i in range(X.shape[1]) if i not in range(sl.start, sl.stop)]
            Xl = X[:, keep] if keep else np.ones((n, 1))
            e_lomo = _predict_pos(models.lomo_clfs[j], Xl, models.prevalence, clip)
            instance_gap[:, j] = e_full - e_lomo
    if models.has_lomo:
        abs_gap = np.abs(instance_gap)
        instance_share = abs_gap / (abs_gap.sum(axis=1, keepdims=True) + EPS)
    else:
        instance_share = np.repeat(models.pi.reshape(1, -1), n, axis=0)
        instance_gap[:] = 0.0
    z = logit(e_block, clip=clip) * models.pi.reshape(1, -1)
    return {
        "e_full": e_full,
        "e_block": e_block,
        "instance_gap": instance_gap,
        "instance_share": instance_share,
        "z": z,
    }


def route_from_shares(
    names: Sequence[str],
    pi: np.ndarray,
    instance_share: np.ndarray,
    catalog: Dict[str, ToolChannel],
    config: RouterConfig | None = None,
    *,
    e_hat: float | None = None,
    z_row: np.ndarray | None = None,
) -> RoutingDecision:
    config = config or RouterConfig()
    names = list(names)
    pi = relu_normalize(np.asarray(pi, dtype=float).reshape(-1))
    s = relu_normalize(np.asarray(instance_share, dtype=float).reshape(-1))
    r = blend_shares(pi, s, lam=config.lam).reshape(-1)
    M = len(names)
    order = list(np.argsort(-r))
    read_order = [names[i] for i in order]
    top_i = int(np.argmax(r))
    top = names[top_i]
    rmax = float(r[top_i])
    ent = simplex_entropy(r)
    ent_cut = (
        float(config.entropy_abstain)
        if config.entropy_abstain is not None
        else 0.92 * float(np.log(max(M, 2)))
    )
    abstain = bool(rmax < config.tau_lo or ent >= ent_cut)
    k_cap = 1 if (not abstain and rmax >= config.tau_hi) else int(config.k_max)

    ranked = sorted(
        (
            (
                float(r[i]) / max(float(catalog[names[i]].cost), EPS),
                int(i),
                names[i],
                catalog[names[i]],
            )
            for i in range(M)
        ),
        reverse=True,
    )
    enabled: List[ToolChannel] = []
    spent = 0.0
    if not abstain:
        for _score, i, name, ch in ranked:
            if float(r[i]) < config.floor_enable:
                continue
            if len(enabled) >= k_cap:
                break
            if spent + ch.cost > config.cost_budget + 1e-9:
                continue
            enabled.append(ch)
            spent += ch.cost
        if not enabled:
            enabled = [catalog[top]]
            spent = float(enabled[0].cost)

    enabled_ids = [c.tool_id for c in enabled]
    enabled_names = [c.name for c in enabled]
    disabled = [catalog[n].tool_id for n in names if n not in enabled_names]
    ask = [top] if abstain else []
    must = None if abstain else top

    pop_blocks, pop_mode = population_pack(names, pi, config=config)
    if abstain:
        pack_mode = "hitl_abstain"
        packed_blocks: List[str] = []
        rationale = (
            f"Neither π nor s localizes a channel (H(r)={ent:.2f}, r_max={rmax:.2f}). "
            f"Do not drop blocks. HITL: ask for '{top}'. If a score is required, pack_all={list(names)}."
        )
    elif rmax >= config.tau_hi:
        pack_mode = "same_expert"
        packed_blocks = packed_in_read_order(read_order, enabled_names)
        rationale = (
            f"Same expert '{top}' (r={rmax:.2f} ≥ τ_hi). Pack only that block; CoT must cite {top}."
        )
    else:
        pack_mode = "concat_k"
        packed_blocks = packed_in_read_order(read_order, enabled_names)
        rationale = (
            f"Concat {packed_blocks} in read-order (one reader, k sections). "
            f"Population template would be {pop_mode}:{pop_blocks}."
        )

    pack_pi = {n: round(float(v), 4) for n, v in zip(names, pi)}
    pack_s = {n: round(float(v), 4) for n, v in zip(names, s)}
    pack_r = {n: round(float(v), 4) for n, v in zip(names, r)}
    packet = {
        "population_shift_pi": pack_pi,
        "this_example_s": pack_s,
        "blend_r": pack_r,
        "read_order": read_order,
        "tools_enabled": enabled_ids,
        "tools_disabled": disabled,
        "packed_blocks": packed_blocks,
        "pack_mode": pack_mode,
        "population_pack": {"mode": pop_mode, "blocks": pop_blocks},
        "abstain": abstain,
        "ask_missing": ask,
        "critic_must_cite": must,
        "e_hat": None if e_hat is None else round(float(e_hat), 4),
        "Z": None
        if z_row is None
        else {n: round(float(v), 4) for n, v in zip(names, np.asarray(z_row).reshape(-1))},
        "rationale": rationale,
    }
    return RoutingDecision(
        read_order=read_order,
        tools_enabled=enabled_ids,
        tools_disabled=disabled,
        abstain=abstain,
        ask_missing=ask,
        critic_must_cite=must,
        blend=pack_r,
        pi=pack_pi,
        instance_share=pack_s,
        entropy=round(ent, 4),
        rationale=rationale,
        packed_blocks=packed_blocks,
        pack_mode=pack_mode,
        packet=packet,
    )


class GapGuidedRouter:
    def __init__(
        self,
        models: FrozenGapModels,
        *,
        catalog: Dict[str, ToolChannel] | None = None,
        config: RouterConfig | None = None,
    ):
        self.models = models
        self.catalog = catalog or catalog_for_spec(models.spec)
        self.config = config or RouterConfig()
        if not models.has_lomo:
            self.config = RouterConfig(**{**asdict(self.config), "lam": 1.0})

    def route_row(self, scores: dict, i: int) -> RoutingDecision:
        return route_from_shares(
            self.models.names,
            self.models.pi,
            scores["instance_share"][i],
            self.catalog,
            self.config,
            e_hat=float(scores["e_full"][i]),
            z_row=scores["z"][i],
        )

    def route(self, X: np.ndarray) -> List[RoutingDecision]:
        scores = predict_gap_scores(self.models, X)
        return [self.route_row(scores, i) for i in range(len(X))]


def render_system_prompt(decision: RoutingDecision) -> str:
    """Copy-paste system prompt. The model never sees W; only π, s, r, tools."""
    must = decision.critic_must_cite or "(none — abstain)"
    lines = [
        "You are answering under a measured distribution shift between two batches.",
        "A frozen domain model decomposed the shift into human channels B_m.",
        "Use the packet below as a *hard routing prior*, not as a hint you may ignore.",
        "",
        "Population gap shares π (prior over which channel drifted):",
        "  " + json.dumps(decision.pi, sort_keys=True),
        "This-example instance LOMO shares s_m(x) = |ê(x)−ê(x_{-m})| / Σ|·|:",
        "  " + json.dumps(decision.instance_share, sort_keys=True),
        "Blend r = λπ + (1−λ)s used for routing:",
        "  " + json.dumps(decision.blend, sort_keys=True),
        "",
        f"Read order (inspect in this order): {decision.read_order}",
        f"Tools you MAY call: {decision.tools_enabled}",
        f"Tools you must NOT call: {decision.tools_disabled}",
        f"Pack mode: {decision.pack_mode}",
        f"User context concatenates ONLY these blocks, in this order: {decision.packed_blocks}",
        "Omitted blocks are not in the window (no empty headers).",
        f"Abstain / pack-all: {str(decision.abstain).lower()}",
        f"If abstain, ask the user for: {decision.ask_missing}",
        f"critic_must_cite: {must}",
        "",
        "Rules:",
        "1. Read packed_blocks in order. Do not retrieve omitted channels.",
        "2. First reasoning sentence names the packed expert and cites π (prior) and s.",
        "3. If critic_must_cite is set, evidence must come from that channel.",
        "4. If abstain is true, do not pack a subset; request ask_missing and stop.",
        "5. After the answer, emit one JSON line: "
        '{"cited": ["<channel>"], "confidence": 0-1}.',
        "",
        f"Router rationale: {decision.rationale}",
    ]
    return "\n".join(lines)


def openai_tool_schema(decision: RoutingDecision, catalog: Dict[str, ToolChannel]) -> List[dict]:
    """Filter the function-calling schema to enabled tools only."""
    by_id = {ch.tool_id: ch for ch in catalog.values()}
    out = []
    for tid in decision.tools_enabled:
        ch = by_id[tid]
        r = decision.blend.get(ch.name, 0.0)
        out.append(
            {
                "type": "function",
                "function": {
                    "name": ch.tool_id,
                    "description": f"{ch.description} Gap blend r_{ch.name}={r:.2f}. {ch.action}",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "note": {
                                "type": "string",
                                "description": "Optional focus inside this channel.",
                            }
                        },
                    },
                },
            }
        )
    return out


def critic_check(decision: RoutingDecision, cited: Sequence[str]) -> CriticVerdict:
    cited_list = [str(c).strip() for c in cited if str(c).strip()]
    must = decision.critic_must_cite
    if decision.abstain:
        ok = True
        mismatch = False
        reask = ""
    elif must is None:
        ok = True
        mismatch = False
        reask = ""
    else:
        mismatch = must not in cited_list
        ok = not mismatch
        reask = (
            ""
            if ok
            else (
                f"Your cited channels {cited_list} do not include '{must}', "
                f"which has the highest gap blend r={decision.blend.get(must)}. "
                f"Re-read {must} (tools {decision.tools_enabled}) and answer again. "
                f"Do not lead with {cited_list}."
            )
        )
    return CriticVerdict(
        ok=ok,
        cited=cited_list,
        must_cite=must,
        mismatch=mismatch,
        reask_prompt=reask,
    )


def _auc(y: np.ndarray, s: np.ndarray) -> float:
    y = _as_1d(y).astype(int)
    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, s))


def budgeted_inference_eval(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    *,
    gt: str | None = None,
    test_size: float = 0.35,
    seed: int = 0,
    n_estimators: int = 50,
    light: bool = False,
    lam: float = 0.4,
) -> dict:
    """Train-frozen experts; test-time routing under a one-block budget.

    Each block RF is an expert that only 'reads' B_m — the same constraint as
    packing one modality into an LLM context window or enabling one tool.
    """
    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    Xtr, Xte, Wtr, Wte = train_test_split(
        X, W, test_size=test_size, random_state=seed, stratify=W,
    )
    models = fit_frozen_gap(
        Xtr, Wtr, spec, seed=seed, n_estimators=n_estimators, light=light,
    )
    scores = predict_gap_scores(models, Xte)
    e_block = scores["e_block"]
    s = scores["instance_share"]
    r = blend_shares(models.pi, s, lam=lam if models.has_lomo else 1.0)
    names = models.names
    M = len(names)
    rng = np.random.default_rng(seed + 13)

    pi_i = int(np.argmax(models.pi))
    vimp_i = int(np.argmax(models.pi_vimp))
    inst_i = np.argmax(s, axis=1)
    blend_i = np.argmax(r, axis=1)
    rand_fixed_i = int(rng.integers(0, M))
    rand_row_i = rng.integers(0, M, size=len(Wte))

    def pick(idx) -> np.ndarray:
        if np.isscalar(idx) or (isinstance(idx, (int, np.integer))):
            return e_block[:, int(idx)]
        idx = np.asarray(idx, dtype=int)
        return e_block[np.arange(len(idx)), idx]

    policies = {
        "unconstrained": scores["e_full"],
        "pi_top1": pick(pi_i),
        "vimp_top1": pick(vimp_i),
        "instance_top1": pick(inst_i),
        "blend_top1": pick(blend_i),
        "random_fixed": pick(rand_fixed_i),
        "random_row": pick(rand_row_i),
        "pi_pool": e_block @ models.pi.reshape(-1),
        "uniform_pool": e_block.mean(axis=1),
    }
    if gt is not None and gt in names:
        policies["oracle_gt"] = pick(names.index(gt))

    aucs = {k: round(_auc(Wte, v), 4) for k, v in policies.items()}

    def hit(idx) -> float:
        if gt is None or gt not in names:
            return float("nan")
        g = names.index(gt)
        if np.isscalar(idx) or isinstance(idx, (int, np.integer)):
            return 1.0 if int(idx) == g else 0.0
        return float(np.mean(np.asarray(idx) == g))

    hits = {
        "pi_top1": round(hit(pi_i), 4),
        "vimp_top1": round(hit(vimp_i), 4),
        "instance_top1": round(hit(inst_i), 4),
        "blend_top1": round(hit(blend_i), 4),
        "random_fixed": round(hit(rand_fixed_i), 4),
        "random_row": round(hit(rand_row_i), 4),
    }
    cfg = RouterConfig()
    pop_E, pop_mode = population_pack(names, models.pi, config=cfg)
    order_pi = [names[i] for i in np.argsort(-models.pi)]
    pack_sets = {
        "same_expert_pi": [names[pi_i]],
        "same_expert_vimp": [names[vimp_i]],
        "same_expert_random": [names[rand_fixed_i]],
        "concat_top2_pi": order_pi[: min(2, M)],
        "pack_all": list(names),
        "abstain_policy": pop_E,
    }
    if gt is not None and gt in names:
        pack_sets["same_expert_oracle"] = [gt]
    pack_auc: Dict[str, float] = {}
    pack_width: Dict[str, int] = {}
    for key, E in pack_sets.items():
        a, w = holdout_packed_auc(
            Xtr, Wtr, Xte, Wte, spec, E, seed=seed + 21, n_estimators=n_estimators,
        )
        pack_auc[key] = round(float(a), 4)
        pack_width[key] = int(w)

    return {
        "n_train": int(len(Wtr)),
        "n_test": int(len(Wte)),
        "pi": models.pack_pi(),
        "pi_vimp": {n: round(float(v), 4) for n, v in zip(names, models.pi_vimp)},
        "pi_entropy": round(simplex_entropy(models.pi), 4),
        "same_expert": names[pi_i],
        "population_pack": {"mode": pop_mode, "blocks": pop_E},
        "selected_block": {
            "pi_top1": names[pi_i],
            "vimp_top1": names[vimp_i],
            "random_fixed": names[rand_fixed_i],
        },
        "auc": aucs,
        "hit_gt": hits,
        "pack_auc": pack_auc,
        "pack_width": pack_width,
        "gt": gt,
        "has_lomo": models.has_lomo,
        "lam": lam if models.has_lomo else 1.0,
        "note": (
            "auc.*_top1 is the frozen specialist ê_m (same expert reads only B_m). "
            "pack_auc is one RF on concatenated enabled columns (context packing). "
            "pi_pool mixes specialist scores — that is not packing."
        ),
    }


def worked_example_packet(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    *,
    seed: int = 0,
    n_estimators: int = 50,
    query_preview: str = "[features omitted; sidecar carries Z and s]",
) -> dict:
    """One end-to-end packet: frozen scores → route → prompt → critic demo."""
    models = fit_frozen_gap(X, W, spec, seed=seed, n_estimators=n_estimators, light=False)
    router = GapGuidedRouter(models)
    scores = predict_gap_scores(models, X)
    r = blend_shares(models.pi, scores["instance_share"], lam=router.config.lam)
    i = int(np.argmax(r.max(axis=1)))
    decision = router.route_row(scores, i)
    prompt = render_system_prompt(decision)
    tools = openai_tool_schema(decision, router.catalog)
    packed_user = pack_user_context(
        {n: f"<{n} payload>" for n in spec.names},
        decision,
    )
    wrong = critic_check(decision, ["text"] if decision.critic_must_cite != "text" else ["valence"])
    right = critic_check(
        decision,
        [decision.critic_must_cite] if decision.critic_must_cite else [],
    )
    return {
        "row_index": i,
        "decision": decision.as_dict(),
        "system_prompt": prompt,
        "user_message": packed_user if packed_user else query_preview,
        "openai_tools": tools,
        "critic_on_wrong_cite": asdict(wrong),
        "critic_on_correct_cite": asdict(right),
        "pi": models.pack_pi(),
    }
