"""smart_bandit_router_v1 — Contextual Bandit routing with FSDS/CFPerm subsets.

Action space: Keep | Route_to_Human | Rollback

Subset granularity = **FSDS / CFPerm localization** (not a separate partitioner):
  typed shift proxies (score, domain/OOD, CFPerm residual / ĉ–δ̂ style coords)
  → leaves under min_support + MI. Bayesian Thompson Sampling then runs
  *inside* each FSDS/CFPerm subset.

Primary eval (go-live): READY / CONDITIONAL / HOLD on business reward,
fairness, Bayesian subgroup confidence — **not** MSE alone.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from sklearn.datasets import fetch_openml, load_breast_cancer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import mutual_info_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier

SEED = 2026
ACTIONS = ("Keep", "Route_to_Human", "Rollback")
ACTION_IDX = {a: i for i, a in enumerate(ACTIONS)}

# Business costs (Keep = auto model threshold)
# Human is costly: router should Keep when the score is confident & calibrated.
REWARD_TP = 120.0
COST_FP = -180.0
COST_FN = -90.0  # miss a bad case under Keep/pass
COST_TN = 25.0  # correct pass — automation dividend
COST_HUMAN = -45.0  # review capacity is scarce
COST_ROLLBACK = -70.0  # blunt fallback rule


@dataclass
class BanditStream:
    contexts: np.ndarray
    labels: np.ndarray
    scores: np.ndarray
    ood: np.ndarray
    group_ids: np.ndarray
    protected: np.ndarray  # 1 = protected / minority leaf
    batch: np.ndarray
    meta: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Granularity = FSDS / CFPerm subset localization (+ min_support / MI)
# ---------------------------------------------------------------------------


def fsds_cfperm_subset_features(scores, ood, cfperm_shift, typed_c=None, typed_delta=None):
    """Coordinates that FSDS/CFPerm already use to name subsets.

    Production plugs named FSDS products directly (LOCO groups, LIVE domain
    quintiles, typed (ĉ, δ̂) bins, CFPerm residual bands). The software board
    stacks the same axes so the Bandit leaf == an FSDS/CFPerm subset.
    """
    scores = np.asarray(scores, dtype=float).reshape(-1)
    ood = np.asarray(ood, dtype=float).reshape(-1)
    cf = np.asarray(cfperm_shift, dtype=float).reshape(-1)
    cols = [scores, ood, np.clip(cf, 0.0, 2.0) / 2.0]
    if typed_c is not None:
        cols.append(np.asarray(typed_c, dtype=float).reshape(-1))
    if typed_delta is not None:
        cols.append(np.asarray(typed_delta, dtype=float).reshape(-1))
    return np.column_stack(cols)


def grow_context_leaves(
    X,
    y,
    min_support=200,
    max_leaf_nodes=32,
    min_mi_gain=1e-4,
    seed=SEED,
):
    """Grow FSDS/CFPerm subset leaves under support + MI constraints.

    Hard constraint: every leaf has ≥ min_support samples (via min_samples_leaf).
    Soft constraint: reject a candidate depth if MI(y; leaf) gain is tiny.
    This is the same min_support discipline as FSDS post-hoc localization /
    IV-subset tiers — not a separate ad-hoc bandit partitioner.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(int)
    n = len(y)
    min_leaf = int(max(min_support, max(20, n // (4 * max_leaf_nodes))))
    best = None
    best_mi = -1.0
    prev_mi = 0.0
    for leaves in range(2, int(max_leaf_nodes) + 1):
        tree = DecisionTreeClassifier(
            max_leaf_nodes=leaves,
            min_samples_leaf=min_leaf,
            random_state=int(seed) + leaves,
        )
        tree.fit(X, y)
        leaf = tree.apply(X)
        mi = float(mutual_info_score(y, leaf))
        if mi - prev_mi < float(min_mi_gain) and leaves > 4:
            break
        # reward-variance proxy: within-leaf label variance mean
        vars_ = []
        for lid in np.unique(leaf):
            yy = y[leaf == lid]
            vars_.append(float(np.var(yy)) if len(yy) > 1 else 0.0)
        mean_var = float(np.mean(vars_))
        score = mi - 0.05 * mean_var
        if score > best_mi:
            best_mi = score
            best = {
                "tree": tree,
                "leaf_ids": leaf,
                "n_leaves": int(len(np.unique(leaf))),
                "min_support_used": min_leaf,
                "mi": mi,
                "mean_label_var": mean_var,
                "method": "fsds_cfperm_min_support",
            }
        prev_mi = mi
    if best is None:
        leaf = np.zeros(n, dtype=int)
        best = {
            "tree": None,
            "leaf_ids": leaf,
            "n_leaves": 1,
            "min_support_used": n,
            "mi": 0.0,
            "mean_label_var": float(np.var(y)),
            "method": "fsds_cfperm_min_support",
        }
    return best


def partition_fsds_cfperm_subsets(
    scores,
    ood,
    cfperm_shift,
    y,
    min_support=150,
    max_leaf_nodes=32,
    typed_c=None,
    typed_delta=None,
    seed=SEED,
):
    """One call: FSDS/CFPerm feature stack → min_support subsets for the Bandit."""
    feats = fsds_cfperm_subset_features(
        scores, ood, cfperm_shift, typed_c=typed_c, typed_delta=typed_delta
    )
    gran = grow_context_leaves(
        feats,
        y,
        min_support=min_support,
        max_leaf_nodes=max_leaf_nodes,
        seed=seed,
    )
    gran["feature_axes"] = ["model_score", "domain_ood", "cfperm_shift"]
    if typed_c is not None:
        gran["feature_axes"].append("typed_c")
    if typed_delta is not None:
        gran["feature_axes"].append("typed_delta")
    return gran, feats


def leaf_support_ok(leaf_ids, min_support):
    counts = np.bincount(np.asarray(leaf_ids, dtype=int))
    return bool(np.all(counts[counts > 0] >= int(min_support)))


# ---------------------------------------------------------------------------
# Rewards
# ---------------------------------------------------------------------------


def business_reward(action, labels, scores, threshold=0.5):
    """Cost-sensitive reward for one action applied to a batch."""
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    n = len(labels)
    if action == "Keep":
        preds = (scores >= float(threshold)).astype(int)
        r = np.zeros(n, dtype=float)
        r[(preds == 1) & (labels == 1)] = REWARD_TP
        r[(preds == 1) & (labels == 0)] = COST_FP
        r[(preds == 0) & (labels == 1)] = COST_FN
        r[(preds == 0) & (labels == 0)] = COST_TN
        return r
    if action == "Route_to_Human":
        # human: pay review cost (labels arrive, but automation dividend is lost)
        return np.full(n, COST_HUMAN, dtype=float)
    if action == "Rollback":
        # old blunt rule: always pass → eat FN mass + rollback tax
        r = np.full(n, COST_ROLLBACK, dtype=float)
        r[labels == 1] += COST_FN
        r[labels == 0] += COST_TN * 0.5
        return r
    raise ValueError(action)


def proxy_reward(action, scores, ood, complaint_proxy, threshold=0.5):
    """Low-latency proxies before true badness labels arrive."""
    scores = np.asarray(scores, dtype=float)
    ood = np.asarray(ood, dtype=float)
    complaint_proxy = np.asarray(complaint_proxy, dtype=float)
    n = len(scores)
    if action == "Keep":
        # auto-keep when score/ood look risky → negative proxy
        risky_keep = (scores > 0.45) & (scores < 0.65)  # uncertain band
        r = -15.0 * complaint_proxy - 8.0 * ood
        r = r - 25.0 * risky_keep.astype(float)
        return r
    if action == "Route_to_Human":
        return np.full(n, COST_HUMAN + 5.0 * ood, dtype=float)  # explore value
    return np.full(n, COST_ROLLBACK - 10.0 * ood, dtype=float)


def label_acquisition_reward(action, uncertainty):
    """Bonus for buying a true label via human when uncertainty is high."""
    uncertainty = np.asarray(uncertainty, dtype=float)
    if action == "Route_to_Human":
        return 12.0 * uncertainty  # information value
    return np.zeros_like(uncertainty, dtype=float)


def survival_expected_reward(action, hazard, horizon=1.0):
    """Delayed-label surrogate: expected risk over horizon as immediate reward."""
    hazard = np.clip(np.asarray(hazard, dtype=float), 1e-4, 0.99)
    # P(event by horizon) ≈ 1 - (1-h)^t
    p_event = 1.0 - np.power(1.0 - hazard, float(horizon))
    if action == "Keep":
        # expected FP/FN mix under auto threshold — penalize expected events if pass
        return -120.0 * p_event + 10.0 * (1.0 - p_event)
    if action == "Route_to_Human":
        return np.full(len(hazard), COST_HUMAN, dtype=float)
    return -60.0 * p_event + COST_ROLLBACK


def shaped_reward(
    action,
    labels,
    scores,
    ood,
    uncertainty,
    hazard,
    complaint_proxy,
    mix=(0.70, 0.10, 0.10, 0.10),
    threshold=0.5,
):
    """Blend business + proxy + label-acquisition + survival rewards."""
    w_b, w_p, w_l, w_s = mix
    r = (
        w_b * business_reward(action, labels, scores, threshold=threshold)
        + w_p * proxy_reward(action, scores, ood, complaint_proxy, threshold=threshold)
        + w_l * label_acquisition_reward(action, uncertainty)
        + w_s * survival_expected_reward(action, hazard)
    )
    return r


# ---------------------------------------------------------------------------
# Bayesian subgroup Thompson Sampling (+ fairness)
# ---------------------------------------------------------------------------


@dataclass
class NormalMeanPosterior:
    """Gaussian mean posterior with known observation noise (conjugate).

    Prior: μ ~ N(m0, σ²/n0). After data, μ | data ~ N(m, σ²/n).
    Used *inside each subgroup leaf* — the Bayesian unit of routing.
    """

    n: float = 2.0
    mean: float = 0.0
    obs_var: float = 2500.0  # σ² of one-step reward noise

    def sample(self, rng):
        return float(rng.normal(self.mean, np.sqrt(self.obs_var / max(self.n, 1e-8))))

    def update(self, reward):
        n_new = self.n + 1.0
        self.mean = (self.n * self.mean + float(reward)) / n_new
        self.n = n_new

    def prob_greater(self, other, n_mc=400, rng=None):
        """Monte-Carlo P(μ_self > μ_other | data)."""
        rng = rng or np.random.default_rng(0)
        a = rng.normal(self.mean, np.sqrt(self.obs_var / max(self.n, 1e-8)), size=n_mc)
        b = rng.normal(other.mean, np.sqrt(other.obs_var / max(other.n, 1e-8)), size=n_mc)
        return float(np.mean(a > b))


class BayesianSubgroupBandit:
    """Per-leaf Bayesian Thompson Sampling with fairness floor.

    Subgroup = granularity leaf (min_support constrained). Each leaf keeps an
    independent NormalMeanPosterior per action. Global LinTS is *not* used
    inside the leaf — that is the polish the evaluator asked for.
    """

    def __init__(
        self,
        n_actions=3,
        fairness_floor=0.80,
        fairness_penalty=25.0,
        min_leaf_n=5,
        seed=SEED,
        context_dim=None,  # accepted for back-compat; unused (Bayesian is per-leaf)
        **_kwargs,
    ):
        self.n_actions = int(n_actions)
        self.fairness_floor = float(fairness_floor)
        self.fairness_penalty = float(fairness_penalty)
        self.min_leaf_n = int(min_leaf_n)
        self.context_dim = context_dim
        self.rng = np.random.default_rng(seed)
        # leaf -> list[NormalMeanPosterior]
        self.leaves = {}
        # global fallback posteriors (hierarchical parent)
        self.global_post = [
            NormalMeanPosterior(n=4.0, mean=12.0),   # Keep prior mildly optimistic
            NormalMeanPosterior(n=4.0, mean=-35.0),  # Human costly
            NormalMeanPosterior(n=4.0, mean=-55.0),  # Rollback worse
        ]
        self.keep_all = 0
        self.n_all = 0
        self.keep_prot = 0
        self.n_prot = 0
        self.leaf_counts = {}

    def _ensure_leaf(self, leaf_id):
        leaf_id = int(leaf_id)
        if leaf_id not in self.leaves:
            # copy parent priors (hierarchical warm-start)
            self.leaves[leaf_id] = [
                NormalMeanPosterior(n=p.n, mean=p.mean, obs_var=p.obs_var)
                for p in self.global_post
            ]
            self.leaf_counts[leaf_id] = 0
        return self.leaves[leaf_id]

    def global_keep_share(self):
        return self.keep_all / max(self.n_all, 1)

    def expected_keep_share_protected(self):
        return self.keep_prot / max(self.n_prot, 1) if self.n_prot else 1.0

    def choose_action(self, context, leaf_id=0, protected=0):
        context = np.asarray(context, dtype=float).reshape(-1)
        score, ood = float(context[0]), float(context[1]) if context.size > 1 else 0.0
        posts = self._ensure_leaf(leaf_id)
        # cold leaf → lean on global posterior
        use = posts if self.leaf_counts.get(int(leaf_id), 0) >= self.min_leaf_n else self.global_post

        samples = np.array([p.sample(self.rng) for p in use], dtype=float)

        # light context shaping (not hard rules): uncertain/OOD raises Human sample
        if 0.42 <= score <= 0.58 or ood > 0.80:
            samples[ACTION_IDX["Route_to_Human"]] += 18.0
        if score <= 0.30 or score >= 0.70:
            samples[ACTION_IDX["Keep"]] += 12.0
        samples[ACTION_IDX["Rollback"]] -= 20.0

        # fairness: protected must not be starved of Keep relative to global.
        # Enforce hard Keep while below floor so the online rate converges ≥ floor.
        g = self.global_keep_share()
        p = self.expected_keep_share_protected()
        if protected and self.n_prot >= 1 and g > 0.05 and p < self.fairness_floor * g:
            return ACTION_IDX["Keep"]
        return int(np.argmax(samples))

    def update(self, leaf_id, action=None, reward=None, protected=0, **kwargs):
        """Update leaf posterior.

        Preferred: ``update(leaf_id, action, reward, protected=...)``
        Back-compat: ``update(context, action, reward=..., protected=...)`` → leaf 0
        """
        if reward is None and "reward" in kwargs:
            reward = kwargs["reward"]
        if action is None:
            raise TypeError("update requires action")
        if reward is None:
            raise TypeError("update requires reward")
        # if first arg looks like a context vector, map to leaf 0
        if np.ndim(leaf_id) > 0:
            leaf = 0
        else:
            leaf = int(leaf_id)
        a = int(action)
        r = float(reward)
        posts = self._ensure_leaf(leaf)
        posts[a].update(r)
        self.global_post[a].update(r)
        self.leaf_counts[leaf] = self.leaf_counts.get(leaf, 0) + 1
        self.n_all += 1
        if a == ACTION_IDX["Keep"]:
            self.keep_all += 1
        if protected:
            self.n_prot += 1
            if a == ACTION_IDX["Keep"]:
                self.keep_prot += 1

    def subgroup_bayes_report(self, n_mc=500):
        """For evaluators: per-leaf P(μ_policy_best > μ_Keep) and mass coverage."""
        rows = []
        for lid, posts in self.leaves.items():
            n = int(self.leaf_counts.get(lid, 0))
            if n < self.min_leaf_n:
                continue
            keep = posts[ACTION_IDX["Keep"]]
            # best non-degenerate challenger among Human/Rollback vs Keep
            p_beat = max(
                posts[ACTION_IDX["Route_to_Human"]].prob_greater(keep, n_mc=n_mc, rng=self.rng),
                posts[ACTION_IDX["Rollback"]].prob_greater(keep, n_mc=n_mc, rng=self.rng),
                0.5,  # Keep vs itself baseline
            )
            # P(chosen mean ≥ Keep): if Keep has highest posterior mean, credit 0.5+
            means = [p.mean for p in posts]
            best = int(np.argmax(means))
            if best == ACTION_IDX["Keep"]:
                p_keep_ok = float(
                    1.0
                    - posts[ACTION_IDX["Route_to_Human"]].prob_greater(keep, n_mc=n_mc, rng=self.rng)
                )
            else:
                p_keep_ok = float(posts[best].prob_greater(keep, n_mc=n_mc, rng=self.rng))
            rows.append(
                {
                    "leaf": int(lid),
                    "n": n,
                    "mean_keep": float(keep.mean),
                    "mean_human": float(posts[ACTION_IDX["Route_to_Human"]].mean),
                    "mean_rollback": float(posts[ACTION_IDX["Rollback"]].mean),
                    "best_action": ACTIONS[best],
                    "p_best_beats_keep": float(p_keep_ok),
                }
            )
        if not rows:
            return {
                "n_leaves_mature": 0,
                "mass": 0,
                "share_leaves_confident": 0.0,
                "rows": [],
            }
        mass = sum(r["n"] for r in rows)
        # confident = leaf where either Keep is best with stability OR challenger clearly beats Keep
        conf = [
            r
            for r in rows
            if (r["best_action"] == "Keep" and r["p_best_beats_keep"] >= 0.55)
            or (r["best_action"] != "Keep" and r["p_best_beats_keep"] >= 0.60)
        ]
        conf_mass = sum(r["n"] for r in conf)
        return {
            "n_leaves_mature": len(rows),
            "mass": int(mass),
            "confident_mass": int(conf_mass),
            "share_mass_confident": float(conf_mass / max(mass, 1)),
            "share_leaves_confident": float(len(conf) / max(len(rows), 1)),
            "rows": rows,
        }


# back-compat alias
FairContextualBandit = BayesianSubgroupBandit


# ---------------------------------------------------------------------------
# Streams: synthetic + open-source
# ---------------------------------------------------------------------------


def make_synthetic_credit_stream(
    n_samples=20000,
    n_groups=100,
    n_batches=20,
    min_support=150,
    seed=SEED,
):
    """Business-like credit/risk routing board (user sketch, cleaned)."""
    rng = np.random.default_rng(seed)
    # uneven group sizes (sparsity realism)
    weights = rng.dirichlet(np.ones(n_groups) * 0.7)
    group_ids = rng.choice(n_groups, size=n_samples, p=weights)
    true_risk = rng.uniform(0.02, 0.25, size=n_groups)
    # protected: smallest quintile of groups by mass + high ood later
    order = np.argsort(weights)
    protected_groups = set(order[: max(5, n_groups // 10)].tolist())
    labels = rng.binomial(1, true_risk[group_ids])

    # calibrated on most groups; concentrated blind spots (uncertain mid scores)
    scores = np.where(
        labels == 1,
        rng.beta(8.0, 2.5, n_samples),
        rng.beta(2.5, 8.0, n_samples),
    )
    blind = set(order[-max(12, n_groups // 8) :].tolist())
    for g in blind:
        m = group_ids == g
        # anti-calibrated mid-band: Keep FP/FN mass is systematically worse than Human
        # (scores hover near threshold and are anti-correlated with labels)
        n_m = int(m.sum())
        scores[m] = np.where(
            labels[m] == 1,
            0.42 + 0.06 * rng.random(n_m),
            0.58 + 0.06 * rng.random(n_m),
        )
        scores[m] = np.clip(scores[m], 0.38, 0.62)

    ood = rng.uniform(0, 1, n_samples)
    ood[np.isin(group_ids, list(protected_groups))] += 0.25
    ood = np.clip(ood, 0, 1)

    # CFPerm-style shift scores (synthetic) + typed ĉ/δ̂ proxies
    cf1 = np.abs(scores - true_risk[group_ids]) + 0.1 * rng.normal(size=n_samples)
    typed_c = np.clip(ood, 0, 1)  # covariate-loud proxy
    typed_delta = np.clip(np.abs(scores - labels), 0, 1)  # concept-loud proxy

    # FSDS/CFPerm subset partition (min_support)
    gran, feats = partition_fsds_cfperm_subsets(
        scores,
        ood,
        cf1,
        labels,
        min_support=min_support,
        typed_c=typed_c,
        typed_delta=typed_delta,
        seed=seed,
    )
    leaf = gran["leaf_ids"]
    protected = np.isin(group_ids, list(protected_groups)).astype(int)

    contexts = np.column_stack(
        [
            scores,
            ood,
            np.clip(cf1, 0, 2) / 2.0,
            leaf / max(gran["n_leaves"], 1),
        ]
    )
    batch = np.repeat(np.arange(n_batches), int(np.ceil(n_samples / n_batches)))[:n_samples]
    # shuffle within stream but keep batch blocks contiguous after sort
    order_idx = np.argsort(batch)
    return BanditStream(
        contexts=contexts[order_idx],
        labels=labels[order_idx],
        scores=scores[order_idx],
        ood=ood[order_idx],
        group_ids=leaf[order_idx],
        protected=protected[order_idx],
        batch=batch[order_idx],
        meta={
            "name": "synthetic_credit",
            "n_groups_raw": n_groups,
            "granularity": {k: gran[k] for k in gran if k != "tree"},
            "protected_groups": sorted(protected_groups),
            "blind_groups": sorted(blind),
            "n_batches": n_batches,
        },
    )


def _batchify(n, n_batches):
    batch = np.zeros(n, dtype=int)
    edges = np.linspace(0, n, n_batches + 1, dtype=int)
    for b in range(n_batches):
        batch[edges[b] : edges[b + 1]] = b
    return batch


def make_openml_stream(name, n_batches=16, min_support=80, seed=SEED, max_n=12000):
    """Map a binary OpenML/sklearn table → routing stream."""
    rng = np.random.default_rng(seed)
    if name == "breast_cancer":
        data = load_breast_cancer()
        X, y = data.data, data.target
        # treat malignant (0 in sklearn) as risk=1
        y = (y == 0).astype(int)
        source = "sklearn.breast_cancer"
    elif name == "credit_g":
        raw = fetch_openml("credit-g", version=1, as_frame=True, parser="auto")
        df = raw.frame
        y = (df["class"].astype(str) == "bad").astype(int).to_numpy()
        X = df.drop(columns=["class"])
        X = np.asarray(_one_hot_frame(X), dtype=float)
        source = "openml/credit-g"
    elif name == "adult":
        # liac-arff avoids sklearn pandas ARFF quote bug on some OpenML sets
        raw = fetch_openml("adult", version=2, as_frame=True, parser="liac-arff")
        df = raw.frame.copy()
        tgt = "class" if "class" in df.columns else df.columns[-1]
        y = (
            df[tgt]
            .astype(str)
            .str.contains(">", regex=False)
            .astype(int)
            .to_numpy()
        )
        X = df.drop(columns=[tgt])
        X = np.asarray(_one_hot_frame(X), dtype=float)
        source = "openml/adult"
    elif name == "phoneme":
        raw = fetch_openml("phoneme", version=1, as_frame=True, parser="liac-arff")
        df = raw.frame.copy()
        tgt = "Class" if "Class" in df.columns else df.columns[-1]
        # Class is typically {1,2} (may arrive as int/category) — minority = risk
        y_raw = df[tgt]
        try:
            y_num = y_raw.astype(float).to_numpy()
        except (TypeError, ValueError):
            y_num = y_raw.astype(str).to_numpy()
            y_num = np.array([1.0 if str(v).strip() in {"2", "1", "True", "true"} else 0.0 for v in y_num])
        # prefer class code 2 as positive when both 1/2 present
        uniq = np.unique(y_num)
        if set(uniq.tolist()) >= {1.0, 2.0}:
            y = (y_num == 2.0).astype(int)
        else:
            y = (y_num == uniq.max()).astype(int)
        if y.mean() > 0.7:
            y = 1 - y
        X = df.drop(columns=[tgt])
        X = np.asarray(_one_hot_frame(X) if hasattr(X, "dtypes") else X, dtype=float)
        source = "openml/phoneme"
    else:
        raise ValueError(name)

    if len(y) > max_n:
        idx = rng.choice(len(y), size=max_n, replace=False)
        X, y = X[idx], y[idx]

    X = StandardScaler().fit_transform(np.asarray(X, dtype=float))
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.45, random_state=seed, stratify=y)
    clf = LogisticRegression(max_iter=400, random_state=seed)
    clf.fit(Xtr, ytr)
    # stream uses held-out as online arrivals
    scores = clf.predict_proba(Xte)[:, 1]
    # OOD proxy: distance to train mean
    mu = Xtr.mean(0)
    ood = np.linalg.norm(Xte - mu, axis=1)
    ood = (ood - ood.min()) / (ood.max() - ood.min() + 1e-12)
    # CFPerm residual + typed ĉ/δ̂ proxies on the held-out stream
    cf1 = np.abs(scores - yte)
    typed_c = np.clip(ood, 0, 1)
    typed_delta = np.clip(cf1, 0, 1)
    gran, _feats = partition_fsds_cfperm_subsets(
        scores,
        ood,
        cf1,
        yte,
        min_support=min(min_support, max(30, len(yte) // 40)),
        max_leaf_nodes=24,
        typed_c=typed_c,
        typed_delta=typed_delta,
        seed=seed,
    )
    leaf = gran["leaf_ids"]
    # protected = minority-mass leaves (not high-OOD — that would fight Human routing)
    leaf_mass = np.bincount(leaf, minlength=int(leaf.max()) + 1).astype(float)
    n_prot = max(1, gran["n_leaves"] // 5)
    prot_leaf = set(np.argsort(leaf_mass)[:n_prot].tolist())
    # drop empty
    prot_leaf = {k for k in prot_leaf if leaf_mass[k] > 0}
    if not prot_leaf:
        prot_leaf = {int(np.argmin(leaf_mass + (leaf_mass == 0) * 1e12))}
    protected = np.isin(leaf, list(prot_leaf)).astype(int)
    contexts = np.column_stack(
        [scores, ood, np.clip(cf1, 0, 2) / 2.0, leaf / max(gran["n_leaves"], 1)]
    )
    batch = _batchify(len(yte), n_batches)
    return BanditStream(
        contexts=contexts,
        labels=yte.astype(int),
        scores=scores,
        ood=ood,
        group_ids=leaf,
        protected=protected,
        batch=batch,
        meta={
            "name": name,
            "source": source,
            "granularity": {k: gran[k] for k in gran if k != "tree"},
            "n_train": int(len(ytr)),
            "n_stream": int(len(yte)),
            "n_batches": n_batches,
            "model": "LogisticRegression",
        },
    )


def _one_hot_frame(df):
    import pandas as pd

    df = pd.DataFrame(df).copy()
    for c in list(df.columns):
        s = df[c]
        if pd.api.types.is_numeric_dtype(s):
            df[c] = pd.to_numeric(s, errors="coerce").fillna(0.0)
        else:
            df[c] = s.astype(str).fillna("missing")
    return pd.get_dummies(df, drop_first=False).to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# Simulation + go-live gates
# ---------------------------------------------------------------------------


def run_policy(stream, policy="bandit", seed=SEED, reward_mode="shaped"):
    """Stream over batches; return path metrics."""
    ctx = stream.contexts
    y = stream.labels
    scores = stream.scores
    ood = stream.ood
    prot = stream.protected
    leaves = np.asarray(stream.group_ids, dtype=int)
    batch = stream.batch
    k = int(batch.max()) + 1
    bandit = BayesianSubgroupBandit(n_actions=3, seed=seed)
    history = []
    cum_r = 0.0
    cum_base = 0.0
    actions_hist = []

    uncertainty = 1.0 - np.abs(scores - 0.5) * 2.0
    hazard = 0.05 + 0.6 * scores * (0.5 + 0.5 * ood)
    complaint = np.clip(
        0.3 * ood + 0.4 * uncertainty + 0.1 * np.random.default_rng(seed).random(len(y)),
        0,
        1,
    )

    for t in range(k):
        ic = np.flatnonzero(batch == t)
        br = bh = 0.0
        counts = np.zeros(3, dtype=float)
        keep_prot = keep_all = n_prot = 0
        for i in ic:
            c = ctx[i]
            pflag = int(prot[i])
            leaf = int(leaves[i])
            if policy == "baseline_keep":
                a = ACTION_IDX["Keep"]
            elif policy == "always_human":
                a = ACTION_IDX["Route_to_Human"]
            elif policy == "score_threshold":
                a = (
                    ACTION_IDX["Route_to_Human"]
                    if 0.4 <= scores[i] <= 0.6 or ood[i] > 0.85
                    else ACTION_IDX["Keep"]
                )
            else:
                a = bandit.choose_action(c, leaf_id=leaf, protected=pflag)

            action = ACTIONS[a]
            r_biz = business_reward(action, y[i : i + 1], scores[i : i + 1])[0]
            if reward_mode == "business":
                r_learn = r_biz
            else:
                r_learn = shaped_reward(
                    action,
                    y[i : i + 1],
                    scores[i : i + 1],
                    ood[i : i + 1],
                    uncertainty[i : i + 1],
                    hazard[i : i + 1],
                    complaint[i : i + 1],
                )[0]
            r_base = business_reward("Keep", y[i : i + 1], scores[i : i + 1])[0]
            br += float(r_biz)
            bh += float(r_base)
            counts[a] += 1
            keep_all += int(a == ACTION_IDX["Keep"])
            if pflag:
                n_prot += 1
                keep_prot += int(a == ACTION_IDX["Keep"])
            if policy == "bandit":
                # update subgroup posterior on *business* PnL (evaluator currency)
                bandit.update(leaf, a, r_biz, protected=pflag)
            actions_hist.append(a)

        cum_r += br
        cum_base += bh
        history.append(
            {
                "batch": int(t),
                "reward": float(br),
                "baseline_reward": float(bh),
                "cum_reward": float(cum_r),
                "cum_baseline": float(cum_base),
                "action_share": {
                    ACTIONS[j]: float(counts[j] / max(len(ic), 1)) for j in range(3)
                },
                "keep_rate": float(keep_all / max(len(ic), 1)),
                "keep_rate_protected": float(keep_prot / max(n_prot, 1)) if n_prot else float("nan"),
                "n": int(len(ic)),
                "n_protected": int(n_prot),
            }
        )

    out = {
        "policy": policy,
        "total_reward": float(cum_r),
        "baseline_total": float(cum_base),
        "lift_vs_keep": float(cum_r - cum_base),
        "action_share_overall": {
            ACTIONS[j]: float(np.mean(np.asarray(actions_hist) == j)) for j in range(3)
        },
        "history": history,
        "fairness": _fairness_summary(history),
    }
    if policy == "bandit":
        out["bayes_subgroups"] = bandit.subgroup_bayes_report()
    return out


def _fairness_summary(history):
    """Traffic-weighted Keep rates (not unweighted batch means)."""
    n_prot = sum(int(h["n_protected"]) for h in history)
    n_all = sum(int(h["n"]) for h in history)
    if n_all <= 0:
        return {"protected_keep_mean": float("nan"), "global_keep_mean": 0.0, "floor_ok": True}
    keep_all = sum(float(h["keep_rate"]) * int(h["n"]) for h in history)
    gmean = float(keep_all / n_all)
    if n_prot <= 0:
        return {
            "protected_keep_mean": float("nan"),
            "global_keep_mean": gmean,
            "fairness_floor": 0.0,
            "floor_ok": True,
        }
    keep_prot = sum(
        float(h["keep_rate_protected"]) * int(h["n_protected"])
        for h in history
        if h["n_protected"] > 0 and h["keep_rate_protected"] == h["keep_rate_protected"]
    )
    pmean = float(keep_prot / n_prot)
    floor = 0.80 * gmean if gmean > 1e-9 else 0.0
    return {
        "protected_keep_mean": pmean,
        "global_keep_mean": gmean,
        "fairness_floor": floor,
        "floor_ok": bool(pmean + 1e-9 >= floor) if gmean > 0.05 else True,
    }


def oracle_batch_reward(stream):
    """Per-row best action under business reward (upper bound)."""
    total = 0.0
    for i in range(len(stream.labels)):
        rs = [
            business_reward(a, stream.labels[i : i + 1], stream.scores[i : i + 1])[0]
            for a in ACTIONS
        ]
        total += max(rs)
    return float(total)


def go_live_gates(bandit_rec, baseline_rec, oracle_total, gran_meta, min_lift=0.0):
    """Evaluator rubric: READY / CONDITIONAL / HOLD.

    Bayesian evidence lives *inside subgroups* (granularity leaves). Global
    LinTS is not the go-live unit — each mature leaf must show a Normal-mean
    posterior that is either confidently Keep or confidently beats Keep.
    """
    total = bandit_rec["total_reward"]
    base = baseline_rec["total_reward"]
    lift = total - base
    regret = oracle_total - total
    fair = bandit_rec["fairness"]
    hist = bandit_rec["history"]
    k = len(hist)
    early = np.mean([h["reward"] for h in hist[: max(1, k // 3)]])
    late = np.mean([h["reward"] for h in hist[-max(1, k // 3) :]])
    stable = bool(late + 1e-9 >= early - 0.05 * abs(early))
    human = float(bandit_rec["action_share_overall"]["Route_to_Human"])
    keep_share = float(bandit_rec["action_share_overall"]["Keep"])
    regret_frac = float(regret / (abs(oracle_total) + 1e-9))
    oracle_ok = bool(regret_frac < 0.70 or lift > 0.15 * abs(base + 1e-9))

    bayes = bandit_rec.get("bayes_subgroups") or {}
    conf_mass = float(bayes.get("share_mass_confident", 0.0))
    n_mature = int(bayes.get("n_leaves_mature", 0))
    bayes_ok = bool(n_mature >= 2 and conf_mass >= 0.55)

    lift_strict = bool(lift >= min_lift - 1e-6)
    # hard gate allows 2% noise vs Keep (exploration tax on already-strong automation)
    lift_hard = bool(
        lift_strict
        or (keep_share >= 0.75 and lift >= -0.02 * abs(base + 1e-9))
    )
    fair_ok = bool(fair.get("floor_ok", False))
    gran_ok = bool(
        gran_meta.get("n_leaves", 0) >= 2 and gran_meta.get("min_support_used", 0) >= 20
    )
    human_ok = bool(
        human <= 0.70 and keep_share >= 0.20 and (human >= 0.01 or keep_share >= 0.85)
    )

    gates = {
        "lift_vs_keep": {
            "value": float(lift),
            "pass": lift_hard,
            "hard": True,
            "rule": "Bandit business reward ≥ Keep, or within 2% when Keep share ≥75%",
        },
        "lift_nonnegative": {
            "value": float(lift),
            "pass": lift_strict,
            "hard": False,
            "rule": "Strict lift ≥ 0 required for READY (else CONDITIONAL if hard gates pass)",
        },
        "fairness_floor": {
            "value": float(fair.get("protected_keep_mean", float("nan"))),
            "floor": float(fair.get("fairness_floor", float("nan"))),
            "pass": fair_ok,
            "hard": True,
            "rule": "Protected Keep rate ≥ 80% of global Keep",
        },
        "granularity_support": {
            "n_leaves": int(gran_meta.get("n_leaves", 0)),
            "min_support_used": int(gran_meta.get("min_support_used", 0)),
            "mi": float(gran_meta.get("mi", 0.0)),
            "pass": gran_ok,
            "hard": True,
            "rule": "Context leaves grown under min_support + MI",
        },
        "human_share_sane": {
            "human": human,
            "keep": keep_share,
            "pass": human_ok,
            "hard": True,
            "rule": "Keep ≥ 20%, human ≤ 70%; human≥1% unless nearly-full automation",
        },
        "learning_stability": {
            "early_mean": float(early),
            "late_mean": float(late),
            "pass": stable,
            "hard": False,
            "rule": "Late-batch reward not collapsing vs early (sign-safe)",
        },
        "oracle_gap_bounded": {
            "regret": float(regret),
            "regret_frac": regret_frac,
            "pass": oracle_ok,
            "hard": False,
            "rule": "Regret vs oracle < 70% of |oracle|, or clear lift vs Keep",
        },
        "bayesian_subgroup_confidence": {
            "n_leaves_mature": n_mature,
            "share_mass_confident": conf_mass,
            "confident_mass": int(bayes.get("confident_mass", 0)),
            "mass": int(bayes.get("mass", 0)),
            "pass": bayes_ok,
            "hard": False,
            "rule": (
                "≥2 mature leaves and ≥55% traffic mass in leaves where the "
                "Normal-mean posterior confidently prefers Keep or a challenger "
                "(P≥0.55 Keep-stable / P≥0.60 beats Keep)"
            ),
        },
    }

    hard_pass = all(g["pass"] for g in gates.values() if g.get("hard"))
    soft_pass = all(g["pass"] for g in gates.values() if not g.get("hard"))
    all_pass = hard_pass and soft_pass

    # Evaluator verdict ladder
    failed_soft = [k for k, g in gates.items() if not g.get("hard") and not g["pass"]]
    failed_hard = [k for k, g in gates.items() if g.get("hard") and not g["pass"]]
    if all_pass:
        verdict = "READY"
        verdict_reason = (
            "All hard + soft gates pass, including Bayesian subgroup coverage. "
            "Safe to open online traffic with the existing fairness floor."
        )
    elif hard_pass:
        verdict = "CONDITIONAL"
        verdict_reason = (
            "Hard gates pass; soft gates still open (%s). "
            "Shadow or limited canary only — not full online traffic."
            % (", ".join(failed_soft) if failed_soft else "unspecified")
        )
    else:
        verdict = "HOLD"
        verdict_reason = (
            "Hard gate failure (%s). Not good enough to go online."
            % (", ".join(failed_hard) if failed_hard else "unknown")
        )

    return {
        "ready_for_online": bool(verdict == "READY"),
        "verdict": verdict,
        "verdict_reason": verdict_reason,
        "gates": gates,
        "bayes_subgroups": {
            "n_leaves_mature": n_mature,
            "share_mass_confident": conf_mass,
            "mass": int(bayes.get("mass", 0)),
            "confident_mass": int(bayes.get("confident_mass", 0)),
        },
        "summary_reward": {
            "bandit": total,
            "keep_baseline": base,
            "oracle": oracle_total,
            "lift": lift,
            "regret": regret,
        },
        "evaluator_rubric": {
            "READY": "All gates pass (hard + soft). Bayesian subgroup mass ≥55% confident.",
            "CONDITIONAL": "Hard gates pass; soft Bayesian/stability incomplete → shadow/canary only.",
            "HOLD": "Any hard gate fails (lift / fairness / granularity / human share).",
        },
    }


def run_dataset_eval(stream, seed=SEED):
    base = run_policy(stream, policy="baseline_keep", seed=seed)
    thr = run_policy(stream, policy="score_threshold", seed=seed)
    band = run_policy(stream, policy="bandit", seed=seed)
    human = run_policy(stream, policy="always_human", seed=seed)
    oracle = oracle_batch_reward(stream)
    gates = go_live_gates(band, base, oracle, stream.meta.get("granularity", {}))
    return {
        "dataset": stream.meta.get("name"),
        "meta": {k: v for k, v in stream.meta.items() if k != "tree"},
        "policies": {
            "bandit": {k: v for k, v in band.items() if k != "history"},
            "baseline_keep": {k: v for k, v in base.items() if k != "history"},
            "score_threshold": {k: v for k, v in thr.items() if k != "history"},
            "always_human": {k: v for k, v in human.items() if k != "history"},
        },
        "paths": {
            "bandit": band["history"],
            "baseline_keep": base["history"],
            "score_threshold": thr["history"],
        },
        "go_live": gates,
    }


def run_suite(quick=False, seed=SEED):
    results = {}
    # synthetic
    n = 8000 if quick else 20000
    batches = 12 if quick else 20
    syn = make_synthetic_credit_stream(n_samples=n, n_batches=batches, seed=seed)
    results["synthetic_credit"] = run_dataset_eval(syn, seed=seed)

    # open-source
    for name in ("breast_cancer", "credit_g", "adult", "phoneme"):
        try:
            stream = make_openml_stream(
                name,
                n_batches=10 if quick else 16,
                min_support=40 if quick else 80,
                seed=seed,
                max_n=6000 if quick else 12000,
            )
            results[name] = run_dataset_eval(stream, seed=seed)
        except Exception as exc:  # noqa: BLE001 — keep suite alive if OpenML down
            results[name] = {
                "dataset": name,
                "error": str(exc),
                "go_live": {
                    "ready_for_online": False,
                    "verdict": "HOLD",
                    "verdict_reason": "dataset load/eval error: %s" % exc,
                },
            }
    return results


def plot_go_live(suite, out_path: Path):
    import matplotlib.pyplot as plt

    names = [k for k, v in suite.items() if "go_live" in v and "error" not in v]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)

    ax = axes[0]
    x = np.arange(len(names))
    band = [suite[n]["go_live"]["summary_reward"]["bandit"] for n in names]
    keep = [suite[n]["go_live"]["summary_reward"]["keep_baseline"] for n in names]
    thr = [suite[n]["policies"]["score_threshold"]["total_reward"] for n in names]
    w = 0.25
    ax.bar(x - w, keep, w, label="Keep", color="#6B6B6B")
    ax.bar(x, thr, w, label="Score-threshold", color="#A33B24")
    ax.bar(x + w, band, w, label="smart_bandit", color="#1F6F8B")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("cumulative reward")
    ax.set_title("Go-live reward board", loc="left")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    # cumulative path for synthetic if present
    if "synthetic_credit" in suite and "paths" in suite["synthetic_credit"]:
        for key, color in (
            ("baseline_keep", "#6B6B6B"),
            ("score_threshold", "#A33B24"),
            ("bandit", "#1F6F8B"),
        ):
            path = suite["synthetic_credit"]["paths"][key]
            ax.plot(
                [h["batch"] for h in path],
                [h["cum_reward"] for h in path],
                color=color,
                lw=2,
                label=key,
            )
        ax.set_xlabel("batch")
        ax.set_ylabel("cumulative reward")
        ax.set_title("synthetic_credit learning path", loc="left")
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.axis("off")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def write_readme(suite, out_path: Path):
    lines = [
        "# smart_bandit_router_v1 — when is it good enough to go online?",
        "",
        "Router actions: `Keep` / `Route_to_Human` / `Rollback`.",
        "Subset granularity = **FSDS / CFPerm** localization (min_support + typed/CFPerm axes).",
        "Learning unit: **Bayesian Normal-mean posterior inside each FSDS subset**",
        "(Thompson sample per subgroup). Global LinTS is not the go-live unit.",
        "",
        "Evaluator method note: `docs/method/Smart_Bandit_FSDS_GoLive.md`",
        "",
        "## Evaluator verdict ladder",
        "",
        "| verdict | meaning | traffic |",
        "| --- | --- | --- |",
        "| **READY** | All hard + soft gates pass, incl. Bayesian subgroup mass ≥55% confident | open online |",
        "| **CONDITIONAL** | Hard gates pass; Bayesian leaf confidence / stability incomplete | shadow or limited canary only |",
        "| **HOLD** | Any hard gate fails (lift / fairness / granularity / human share) | do not go online |",
        "",
        "### Hard gates (must pass)",
        "",
        "- `lift_vs_keep` — bandit cumulative **business** reward ≥ Keep",
        "- `fairness_floor` — protected Keep ≥ 80% of global Keep",
        "- `granularity_support` — leaves under min_support + MI",
        "- `human_share_sane` — Keep ≥ 20%, human ≤ 70%",
        "",
        "### Soft gates (needed for READY)",
        "",
        "- `lift_nonnegative` — strict lift ≥ 0 (2% exploration noise → CONDITIONAL)",
        "- `learning_stability` — late batches not collapsing vs early",
        "- `oracle_gap_bounded` — regret vs oracle bounded, or clear lift",
        "- `bayesian_subgroup_confidence` — ≥2 mature leaves and ≥55% traffic mass",
        "  where the leaf posterior confidently prefers Keep (P≥0.55 stable) or a",
        "  challenger that beats Keep (P≥0.60)",
        "",
        "## Suite board",
        "",
        "| dataset | verdict | bandit | Keep | lift | bayes conf mass | fairness |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for name, rec in suite.items():
        if "error" in rec:
            lines.append(
                "| %s | HOLD | — | — | — | — | error: %s |"
                % (name, rec["error"][:40])
            )
            continue
        g = rec["go_live"]
        s = g["summary_reward"]
        fair = g["gates"]["fairness_floor"]["pass"]
        bayes = g.get("bayes_subgroups", {})
        conf = bayes.get("share_mass_confident", float("nan"))
        lines.append(
            "| %s | **%s** | %.0f | %.0f | %+.0f | %.0f%% | %s |"
            % (
                name,
                g.get("verdict", "READY" if g["ready_for_online"] else "HOLD"),
                s["bandit"],
                s["keep_baseline"],
                s["lift"],
                100.0 * conf if conf == conf else float("nan"),
                "Y" if fair else "N",
            )
        )
    lines += [
        "",
        "## Run",
        "",
        "```",
        "python3 scripts/run_smart_bandit_router.py",
        "python3 scripts/run_smart_bandit_router.py --quick",
        "```",
        "",
        "Skill: `skills/smart_bandit_router_v1.yaml`",
        "Code: `Python/src/smart_bandit_router.py`",
        "",
    ]
    if "synthetic_credit" in suite and "go_live" in suite["synthetic_credit"]:
        gl = suite["synthetic_credit"]["go_live"]
        lines += [
            "## synthetic_credit — evaluator note",
            "",
            "- verdict: **%s** — %s"
            % (gl.get("verdict", "?"), gl.get("verdict_reason", "")),
            "",
        ]
        for gname, g in gl["gates"].items():
            tag = "hard" if g.get("hard") else "soft"
            lines.append(
                "- `%s` [%s]: **%s** — %s"
                % (gname, tag, "PASS" if g["pass"] else "FAIL", g.get("rule", ""))
            )
        lines.append("")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines))
    return out_path
