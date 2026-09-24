"""Round-1 probe: OnlineRFPerm-style evaluator monitoring on five classic agent tasks.

Tasks (programmatic, auto-scored; no external LLM):
  1. Game of 24          — Tree-of-Thoughts (Yao et al., 2023)
  2. Blocksworld         — PlanBench (Valmeekam et al., 2023)
  3. MiniGrid DoorKey    — BabyAI / MiniGrid (Chevalier-Boisvert et al.)
  4. Two-hop QA          — HotpotQA structure (Yang et al., 2018)
  5. Attribute shopping  — WebShop structure (Yao et al., 2022)

The search policy is a small beam. A noisy state evaluator is calibrated on a
reference prefix and then suffers concept drift (it starts trusting a spurious
cue). The monitor follows OnlineRFPerm's online procedure: batch loss of the
evaluator, empirical p-values against the reference stream, then LORD with a
summable spending sequence. On rejection the policy falls back to a stationary
heuristic that was available the whole time. An oracle beam is a ceiling only.

Localization follows the FSDS covariate-shift half: which step feature's
mean differs between the reference and drifted batches, judged by a
permutation test (CFPerm's permute-the-batch-label pattern, per feature).
"""

from __future__ import annotations

import json
import math
import random
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "results_round7.json"
# Weight on the spurious cue after drift. Below 1 so the calibrated score remains
# in the judge and post-drift success is mixed rather than identically zero.
# Overwritten per task in main. Weight on the spurious cue after the drift point.
CURRENT_MIX = 0.55
# Per-episode spread around CURRENT_MIX. DoorKey is the same start every episode,
# so a single mix is all-or-nothing; a small spread puts some episodes on each side.
MIX_JITTER = 0.0


def drifted_judge(calibrated: float, spurious: float, t: int, drift_at: int) -> float:
    if t < drift_at:
        return calibrated
    w = CURRENT_MIX
    if MIX_JITTER:
        w = min(1.0, max(0.0, CURRENT_MIX + ((t % 5) - 2) * MIX_JITTER))
    return (1.0 - w) * calibrated + w * spurious
SEED = 2026


# ---------------------------------------------------------------------------
# Online monitor: empirical p-values + LORD (summable gamma)
# ---------------------------------------------------------------------------

def empirical_pvals(losses: np.ndarray, burnin: int) -> np.ndarray:
    """Match OnlineRFPerm::empirical_pval: large loss -> small p-value."""
    losses = np.asarray(losses, dtype=float)
    n = len(losses)
    pvals = np.ones(n)
    for i in range(burnin, n):
        past = losses[burnin:i]
        if len(past) == 0:
            pvals[i] = 1.0
        else:
            pvals[i] = (1.0 + np.sum(past >= losses[i])) / (len(past) + 1.0)
    return pvals


def gamma_seq(n: int) -> np.ndarray:
    t = np.arange(1, n + 1, dtype=float)
    g = 1.0 / (t * np.log(np.e * t) ** 2)
    return g / g.sum()


def lord_reject(pvals: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """LORD as in OnlinePermOOB_core.R, with summable gamma so wealth is a test."""
    n = len(pvals)
    gamma = gamma_seq(n)
    alpha_t = np.zeros(n)
    reject = np.zeros(n, dtype=bool)
    rej_times: list[int] = []
    wealth = alpha
    for t in range(n):
        if not np.isfinite(pvals[t]) or wealth <= 1e-12:
            continue
        a = gamma[t] * alpha
        for tau in rej_times:
            lag = t - tau
            if 0 < lag <= n:
                a += gamma[lag - 1] * alpha
        alpha_t[t] = min(a, wealth)
        reject[t] = pvals[t] <= alpha_t[t]
        wealth = wealth - alpha_t[t] + (alpha if reject[t] else 0.0)
        if reject[t]:
            rej_times.append(t)
    return reject


def first_reject_after(reject: np.ndarray, start: int) -> int | None:
    idx = np.flatnonzero(reject[start:])
    if len(idx) == 0:
        return None
    return int(start + idx[0])


class Gate:
    """Frozen-reference quantile rule, plus a LORD clock that often has no power.

    Permutation p-values on a short stream cannot get small enough for a summable
    spending sequence, so LORD is recorded but does not steer the search.
    The steering rule matches the batch test: alarm after k consecutive losses
    above the reference quantile.
    """

    def __init__(self, burnin: int = 12, k: int = 2, quantile: float = 0.95, trial: int = 4):
        self.burnin = burnin
        self.k = k
        self.quantile = quantile
        self.losses: list[float] = []
        self.switch_at: int | None = None
        self.lord_at: int | None = None
        self.streak = 0
        self.trial = trial
        self.reverted = False
        self.kept = False
        self.ref_heur: list[float] = []

    def fallback_wins(self, t: int) -> bool:
        """Reference heuristic loss versus the judge's recent loss. No post-drift trial."""
        if len(self.ref_heur) < self.burnin or t < self.burnin:
            return False
        recent = self.losses[max(self.burnin, t - 8) : t]
        if len(recent) < 4:
            return False
        return float(np.mean(self.ref_heur)) + 1e-9 < float(np.mean(recent))

    def observe(self, t: int, loss: float) -> None:
        self.losses.append(loss)
        if t >= self.burnin and self.switch_at is None:
            thr = float(np.quantile(self.losses[: self.burnin], self.quantile))
            if loss > thr + 1e-9:
                self.streak += 1
            else:
                self.streak = 0
            if self.streak >= self.k:
                self.switch_at = t + 1
        if (
            self.trial > 0
            and self.switch_at is not None
            and not self.reverted
            and not self.kept
            and len(self.losses) >= self.switch_at + self.trial
        ):
            tried = self.losses[self.switch_at : self.switch_at + self.trial]
            before = self.losses[self.switch_at - self.trial : self.switch_at]
            if before and float(np.mean(tried)) > float(np.mean(before)) + 1e-9:
                self.reverted = True
            else:
                self.kept = True
        if t >= self.burnin and self.lord_at is None:
            pvals = empirical_pvals(np.asarray(self.losses), burnin=5)
            if lord_reject(pvals)[t]:
                self.lord_at = t

    def use_stable(self, t: int, policy: str) -> bool:
        if policy == "stable":
            return True
        if policy == "gated" and self.switch_at is not None and t >= self.switch_at:
            return True
        if policy == "confirm":
            ok = self.fallback_wins(t)
            if ok and self.switch_at is None:
                self.switch_at = t
            return ok
        return False


# ---------------------------------------------------------------------------
# Shared beam search
# ---------------------------------------------------------------------------

@dataclass
class StepView:
    value: float
    judge: float
    stable: float
    spurious: float
    progress: float


def beam_search(root, expand, score_fn, beam: int, depth: int, rng: random.Random):
    """score_fn(state) -> StepView. Returns (final_or_None, chosen-step view, expansions)."""
    frontier = [root]
    chosen = score_fn(root)
    expansions = 0
    for step in range(depth):
        cand = []
        for state in frontier:
            for nxt, done in expand(state):
                expansions += 1
                view = score_fn(nxt)
                cand.append((view.value + rng.uniform(-1e-9, 1e-9), nxt, done, view))
        if not cand:
            break
        cand.sort(key=lambda z: z[0], reverse=True)
        if step == 0:
            chosen = cand[0][3]
        for _, nxt, done, _ in cand[:beam]:
            if done:
                return nxt, chosen, expansions
        frontier = [nxt for _, nxt, _, _ in cand[:beam]]
    return None, chosen, expansions


# ---------------------------------------------------------------------------
# Task 1: Game of 24
# ---------------------------------------------------------------------------

_G24_MEMO: dict[tuple, bool] = {}


def _g24_solvable(nums: tuple[float, ...], target: float = 24.0) -> bool:
    key = tuple(sorted(round(x, 6) for x in nums))
    if key in _G24_MEMO:
        return _G24_MEMO[key]
    if len(key) == 1:
        ok = abs(key[0] - target) < 1e-4
        _G24_MEMO[key] = ok
        return ok
    ok = False
    xs = list(key)
    for i in range(len(xs)):
        for j in range(len(xs)):
            if i == j:
                continue
            a, b = xs[i], xs[j]
            rest = [xs[k] for k in range(len(xs)) if k != i and k != j]
            cands = [a + b, a - b, a * b]
            if abs(b) > 1e-8:
                cands.append(a / b)
            # commute duplicates
            if i > j:
                cands = [a - b] + ([a / b] if abs(b) > 1e-8 else [])
            for v in cands:
                if _g24_solvable(tuple(rest + [v]), target):
                    ok = True
                    break
            if ok:
                break
        if ok:
            break
    _G24_MEMO[key] = ok
    return ok


def _g24_children(nums: tuple[float, ...]):
    xs = list(nums)
    seen = set()
    out = []
    for i in range(len(xs)):
        for j in range(len(xs)):
            if i == j:
                continue
            a, b = xs[i], xs[j]
            rest = [xs[k] for k in range(len(xs)) if k != i and k != j]
            pairs = [(a + b, "+"), (a - b, "-"), (a * b, "*")]
            if abs(b) > 1e-8:
                pairs.append((a / b, "/"))
            for v, _ in pairs:
                key = tuple(sorted(round(x, 6) for x in rest + [v]))
                if key in seen:
                    continue
                seen.add(key)
                done = len(key) == 1 and abs(key[0] - 24) < 1e-4
                out.append((key, done))
    return out


def make_game24(rng: random.Random, n: int):
    puzzles = []
    guard = 0
    while len(puzzles) < n and guard < n * 40:
        guard += 1
        nums = tuple(rng.randint(1, 9) for _ in range(4))
        if _g24_solvable(nums):
            puzzles.append(nums)
    return puzzles


def run_game24(puzzles, drift_at: int, policy: str, beam: int, rng: random.Random):
    losses = []
    successes = []
    feats = []

    def make_score(t, state, use_stable):
        solv = 1.0 if _g24_solvable(state) else 0.0
        # stationary heuristic: any intermediate close to 24, else large partial products
        close = max(math.exp(-abs(x - 24) / 8) for x in state)
        stable = close if len(state) == 1 else 0.45 * close + 0.15
        spurious = sum(state) / (13 * len(state))
        calibrated = 0.85 * solv + 0.15 * stable
        judge = drifted_judge(calibrated, spurious, t, drift_at)
        steer = solv if policy == "oracle" else (stable if use_stable else judge)
        return StepView(steer, judge, stable, spurious, solv)

    gate = Gate(k=2, trial=0 if policy == "confirm" else 4)
    for t, puzzle in enumerate(puzzles):
        if policy == "confirm" and t < gate.burnin:
            held, _, _ = beam_search(
                puzzle, _g24_children,
                lambda state, t=t: make_score(t, state, True),
                beam=beam, depth=3, rng=random.Random(10_000 + t),
            )
            gate.ref_heur.append(0.0 if held is not None else 1.0)
        use_stable = gate.use_stable(t, policy)

        final, chosen, _ = beam_search(
            puzzle, _g24_children, lambda state, t=t, use_stable=use_stable: make_score(t, state, use_stable), beam=beam, depth=3, rng=rng
        )
        y = 1.0 if final is not None else 0.0
        loss = 1.0 - y
        losses.append(loss)
        successes.append(y)
        feats.append(
            {
                "judge": chosen.judge,
                "stable": chosen.stable,
                "spurious": chosen.spurious,
                "progress": chosen.progress,
                "y": y,
            }
        )
        if policy in ("gated", "confirm"):
            gate.observe(t, loss)
    return {
        "success": successes,
        "loss": losses,
        "feats": feats,
        "switch_at": gate.switch_at,
        "reverted": gate.reverted,
        "lord_at": gate.lord_at,
    }


# ---------------------------------------------------------------------------
# Task 2: Blocksworld (4 blocks)
# ---------------------------------------------------------------------------

BLOCKS = ("A", "B", "C", "D")


def _bw_clear(stacks, block):
    for st in stacks:
        if block in st and st[-1] != block:
            return False
    return True


def _bw_holding_stack(stacks, block):
    for i, st in enumerate(stacks):
        if block in st:
            return i
    raise KeyError(block)


def bw_goal():
    return (("A", "B", "C", "D"),)


def bw_success(state) -> bool:
    return state == bw_goal()


def bw_neighbors(state):
    stacks = [list(st) for st in state]
    moves = []
    # move clear block onto another clear block, or to table if not already alone
    clears = [st[-1] for st in stacks if st]
    for b in clears:
        si = _bw_holding_stack(tuple(tuple(s) for s in stacks), b)
        # to table
        if len(stacks[si]) > 1:
            new_stacks = [list(st) for st in stacks]
            new_stacks[si].pop()
            new_stacks.append([b])
            moves.append(_bw_canon(new_stacks))
        for c in clears:
            if c == b:
                continue
            new_stacks = [list(st) for st in stacks]
            new_stacks[si].pop()
            ti = _bw_holding_stack(tuple(tuple(s) for s in new_stacks), c)
            new_stacks[ti].append(b)
            new_stacks = [st for st in new_stacks if st]
            moves.append(_bw_canon(new_stacks))
    # unique
    uniq = []
    seen = set()
    for m in moves:
        if m not in seen:
            seen.add(m)
            uniq.append(m)
    return uniq


def _bw_canon(stacks):
    cleaned = [tuple(st) for st in stacks if st]
    return tuple(sorted(cleaned, key=lambda st: st[0]))


def bw_distances():
    goal = bw_goal()
    dist = {goal: 0}
    q = deque([goal])
    while q:
        cur = q.popleft()
        for nxt in bw_neighbors(cur):
            if nxt not in dist:
                dist[nxt] = dist[cur] + 1
                q.append(nxt)
    return dist


BW_DIST = None


def bw_random_start(rng: random.Random, dist_map, lo=2, hi=4):
    cands = [s for s, d in dist_map.items() if lo <= d <= hi]
    return rng.choice(cands)


def bw_stable(state) -> float:
    """Longest correct prefix of the goal tower, in any stack."""
    goal = ("A", "B", "C", "D")
    best = 0
    for st in state:
        k = 0
        for a, b in zip(st, goal):
            if a != b:
                break
            k += 1
        best = max(best, k)
    return best / 4


def bw_spurious(state) -> float:
    # Graded anti-tower cue: more stacks, and D not buried. Not the goal prefix.
    buried = 0.0
    for st in state:
        if "D" in st:
            buried = st.index("D") / 3
    return 0.6 * (len(state) / 4) + 0.4 * (1.0 - buried)


def run_blocksworld(starts, drift_at, policy, beam, rng, dist_map):
    successes, losses, feats = [], [], []
    gate = Gate(k=2, trial=0 if policy == "confirm" else 4)

    def expand(state):
        out = []
        for nxt in bw_neighbors(state):
            out.append((nxt, bw_success(nxt)))
        return out

    def score_at(state, t, use_stable):
        d = dist_map.get(state, 8)
        progress = 1.0 - d / 8
        stable = bw_stable(state)
        spurious = bw_spurious(state)
        calibrated = 0.8 * progress + 0.2 * stable
        judge = drifted_judge(calibrated, spurious, t, drift_at)
        steer = progress if policy == "oracle" else (stable if use_stable else judge)
        return StepView(steer, judge, stable, spurious, progress)

    for t, start in enumerate(starts):
        if policy == "confirm" and t < gate.burnin:
            held, _, _ = beam_search(
                start, expand,
                lambda state, t=t: score_at(state, t, True),
                beam=beam, depth=4, rng=random.Random(10_000 + t),
            )
            gate.ref_heur.append(0.0 if held is not None else 1.0)
        use_stable = gate.use_stable(t, policy)
        final, chosen, _ = beam_search(
            start, expand, lambda state, t=t, use_stable=use_stable: score_at(state, t, use_stable),
            beam=beam, depth=4, rng=rng,
        )
        y = 1.0 if final is not None else 0.0
        loss = 1.0 - y
        successes.append(y)
        losses.append(loss)
        feats.append(
            {
                "judge": chosen.judge,
                "stable": chosen.stable,
                "spurious": chosen.spurious,
                "progress": chosen.progress,
                "y": y,
            }
        )
        if policy in ("gated", "confirm"):
            gate.observe(t, loss)
    return {
        "success": successes,
        "loss": losses,
        "feats": feats,
        "switch_at": gate.switch_at,
        "reverted": gate.reverted,
        "lord_at": gate.lord_at,
    }


# ---------------------------------------------------------------------------
# Task 3: MiniGrid DoorKey (5x5)
# ---------------------------------------------------------------------------

# 0 empty, 1 wall. Agent starts at (1,1) facing east (0).
# Key at (1,3), door at (3,2), goal at (3,3).
GRID = (
    (1, 1, 1, 1, 1),
    (1, 0, 0, 0, 1),
    (1, 0, 1, 0, 1),
    (1, 0, 0, 0, 1),
    (1, 1, 1, 1, 1),
)
KEY = (1, 3)
DOOR = (3, 2)
GOAL = (3, 3)
DIRS = ((1, 0), (0, 1), (-1, 0), (0, -1))  # E N W S


def _dk_passable(x, y, door_open):
    if not (0 <= x < 5 and 0 <= y < 5):
        return False
    if (x, y) == DOOR:
        return door_open
    return GRID[y][x] == 0


def dk_neighbors(state):
    x, y, d, has_key, door_open = state
    out = []
    # turn left / right
    out.append(((x, y, (d - 1) % 4, has_key, door_open), False))
    out.append(((x, y, (d + 1) % 4, has_key, door_open), False))
    dx, dy = DIRS[d]
    nx, ny = x + dx, y + dy
    if _dk_passable(nx, ny, door_open):
        done = (nx, ny) == GOAL
        out.append(((nx, ny, d, has_key, door_open), done))
    # pickup
    if (x, y) == KEY and not has_key:
        out.append(((x, y, d, 1, door_open), False))
    # toggle door if adjacent and holding key
    if has_key and not door_open and abs(x - DOOR[0]) + abs(y - DOOR[1]) == 1:
        out.append(((x, y, d, has_key, 1), False))
    return out


def dk_distances():
    # reverse BFS is awkward with irreversible key/door; forward from start is enough
    # Precompute dist-to-goal by reverse from goal over the reversible skeleton:
    # search all states forward once from a dummy and record dist via Dijkstra from goal backwards.
    goal_states = []
    for d in range(4):
        for has in (0, 1):
            goal_states.append((GOAL[0], GOAL[1], d, has, 1))
    dist = {s: 0 for s in goal_states}
    q = deque(goal_states)
    # build reverse edges by enumerating the small state space
    states = []
    for x in range(5):
        for y in range(5):
            if GRID[y][x] == 1 and (x, y) != DOOR:
                continue
            for d in range(4):
                for has in (0, 1):
                    for door in (0, 1):
                        if (x, y) == DOOR and door == 0:
                            continue
                        if (x, y) == KEY:
                            pass
                        states.append((x, y, d, has, door))
    rev = {s: [] for s in states}
    for s in states:
        for nxt, _ in dk_neighbors(s):
            if nxt in rev:
                rev[nxt].append(s)
    while q:
        cur = q.popleft()
        for prev in rev.get(cur, []):
            if prev not in dist:
                dist[prev] = dist[cur] + 1
                q.append(prev)
    return dist


def dk_stable(state) -> float:
    x, y, _d, has_key, door_open = state
    if not has_key:
        target, stage = KEY, 0
    elif not door_open:
        target, stage = DOOR, 1
    else:
        target, stage = GOAL, 2
    man = abs(x - target[0]) + abs(y - target[1])
    return stage / 3 + (1.0 / 3.0) * (1.0 - min(man, 8) / 8)


def dk_spurious(state) -> float:
    # Dense distractor: left column, facing west, higher row. Independent of key and door.
    x, y, d, _has_key, _door_open = state
    return (0.45 if x == 1 else 0.0) + (0.35 if d == 2 else 0.0) + 0.05 * y


def run_doorkey(n, drift_at, policy, beam, rng, dist_map):
    start = (1, 1, 0, 0, 0)
    successes, losses, feats = [], [], []
    gate = Gate(k=2, trial=0 if policy == "confirm" else 4)

    def expand(state):
        return dk_neighbors(state)

    def score_at(state, t, use_stable):
        d = dist_map.get(state, 30)
        progress = max(0.0, 1.0 - d / 20)
        stable = dk_stable(state)
        spurious = dk_spurious(state)
        calibrated = 0.85 * progress + 0.15 * stable
        judge = drifted_judge(calibrated, spurious, t, drift_at)
        steer = progress if policy == "oracle" else (stable if use_stable else judge)
        return StepView(steer, judge, stable, spurious, progress)

    for t in range(n):
        if policy == "confirm" and t < gate.burnin:
            held, _, _ = beam_search(
                start, expand, lambda state, t=t: score_at(state, t, True),
                beam=beam, depth=12, rng=random.Random(11_000 + t),
            )
            gate.ref_heur.append(0.0 if held is not None else 1.0)
        use_stable = gate.use_stable(t, policy)
        final, chosen, _ = beam_search(
            start, expand, lambda state, t=t, use_stable=use_stable: score_at(state, t, use_stable),
            beam=beam, depth=12, rng=rng,
        )
        y = 1.0 if final is not None else 0.0
        loss = 1.0 - y
        successes.append(y)
        losses.append(loss)
        feats.append(
            {
                "judge": chosen.judge,
                "stable": chosen.stable,
                "spurious": chosen.spurious,
                "progress": chosen.progress,
                "y": y,
            }
        )
        if policy in ("gated", "confirm"):
            gate.observe(t, loss)
    return {
        "success": successes,
        "loss": losses,
        "feats": feats,
        "switch_at": gate.switch_at,
        "reverted": gate.reverted,
        "lord_at": gate.lord_at,
    }


# ---------------------------------------------------------------------------
# Task 4: HotpotQA-style two-hop retrieval
# ---------------------------------------------------------------------------

# Each question: gold passage ids that must both be selected, in two picks.
# Passages have lexical overlap (stable) and a length distractor (spurious).

PASSAGES = [
    "paris is the capital of france and sits on the seine",
    "marie curie conducted radium research in paris",
    "the seine is a river in france",
    "london is the capital of england on the thames",
    "shakespeare wrote hamlet in london",
    "the thames is a river in england",
    "tokyo is the capital of japan",
    "murakami wrote novels in tokyo",
    "osaka is a city in japan",
    "berlin is the capital of germany",
    "einstein studied in berlin",
    "the spree is a river in berlin",
    "rome is the capital of italy",
    "fibonacci published in pisa near rome",
    "the tiber is a river in rome",
    "madrid is the capital of spain",
    "cervantes wrote don quixote in madrid",
    "the manzanares is a river in madrid",
    "cairo is the capital of egypt",
    "the nile flows through cairo",
]


def hop_questions():
    # (question tokens, gold set)
    pairs = [
        ("where did marie curie conduct research", {0, 1}),
        ("which river is associated with shakespeare city", {3, 4, 5}),
        ("who wrote novels in the capital of japan", {6, 7}),
        ("who studied in the capital of germany", {9, 10}),
        ("which river flows through the capital of egypt", {18, 19}),
        ("who published near the capital of italy", {12, 13}),
        ("who wrote don quixote and in which capital", {15, 16}),
        ("which river runs through paris", {0, 2}),
    ]
    return pairs


def hop_overlap(question: str, pid: int) -> float:
    qt = set(question.split())
    pt = set(PASSAGES[pid].split())
    return len(qt & pt) / max(1, len(qt))


def run_hotpot(n, drift_at, policy, beam, rng):
    bank = hop_questions()
    successes, losses, feats = [], [], []
    gate = Gate(k=2, trial=0 if policy == "confirm" else 4)

    for t in range(n):
        question, gold = bank[t % len(bank)]
        # state is frozenset of picked ids; need |picked ∩ gold| == 2 and len==2
        # gold sets of size 2 are strict; size 3 means any 2 of the chain count as success
        need = 2

        def done_fn(picked):
            return len(picked) >= 2 and len(picked & gold) >= need

        def expand(picked):
            if len(picked) >= 2:
                return []
            out = []
            for pid in range(len(PASSAGES)):
                if pid in picked:
                    continue
                nxt = frozenset(set(picked) | {pid})
                out.append((nxt, done_fn(nxt)))
            return out

        def score_at(picked, t, use_stable):
            if not picked:
                progress = 0.0
                stable = 0.0
                spurious = 0.0
            else:
                progress = len(set(picked) & gold) / need
                stable = float(np.mean([hop_overlap(question, pid) for pid in picked]))
                spurious = float(np.mean([len(PASSAGES[pid].split()) / 12 for pid in picked]))
            calibrated = 0.75 * progress + 0.25 * stable
            judge = drifted_judge(calibrated, spurious, t, drift_at)
            steer = progress if policy == "oracle" else (stable if use_stable else judge)
            return StepView(steer, judge, stable, spurious, progress)

        if policy == "confirm" and t < gate.burnin:
            held, _, _ = beam_search(
                frozenset(), expand, lambda picked, t=t: score_at(picked, t, True),
                beam=beam, depth=2, rng=random.Random(12_000 + t),
            )
            gate.ref_heur.append(0.0 if held is not None else 1.0)
        use_stable = gate.use_stable(t, policy)
        final, chosen, _ = beam_search(
            frozenset(), expand,
            lambda picked, t=t, use_stable=use_stable: score_at(picked, t, use_stable),
            beam=beam, depth=2, rng=rng,
        )
        y = 1.0 if final is not None else 0.0
        loss = 1.0 - y
        successes.append(y)
        losses.append(loss)
        feats.append(
            {
                "judge": chosen.judge,
                "stable": chosen.stable,
                "spurious": chosen.spurious,
                "progress": chosen.progress,
                "y": y,
            }
        )
        if policy in ("gated", "confirm"):
            gate.observe(t, loss)
    return {
        "success": successes,
        "loss": losses,
        "feats": feats,
        "switch_at": gate.switch_at,
        "reverted": gate.reverted,
        "lord_at": gate.lord_at,
    }


# ---------------------------------------------------------------------------
# Task 5: WebShop-style attribute purchase
# ---------------------------------------------------------------------------

# item: (color, price_bucket, category) and review_count distractor
COLORS = ("red", "blue", "black")
PRICES = ("cheap", "mid", "dear")
CATS = ("shoe", "hat", "bag")


def webshop_catalog(rng: random.Random):
    items = []
    for c in COLORS:
        for p in PRICES:
            for cat in CATS:
                reviews = rng.randint(1, 100)
                items.append({"color": c, "price": p, "cat": cat, "reviews": reviews})
    return items


def run_webshop(n, drift_at, policy, beam, rng):
    successes, losses, feats = [], [], []
    gate = Gate(k=2, trial=0 if policy == "confirm" else 4)
    catalog = webshop_catalog(rng)

    for t in range(n):
        # query wants an exact triple; one item matches
        want = {
            "color": COLORS[t % 3],
            "price": PRICES[(t // 3) % 3],
            "cat": CATS[(t // 9) % 3],
        }

        def match_count(item):
            return sum(item[k] == want[k] for k in ("color", "price", "cat"))

        def expand(state):
            # state is picked index or None. one-step buy.
            if state is not None:
                return []
            out = []
            for i, item in enumerate(catalog):
                out.append((i, match_count(item) == 3))
            return out

        def score_at(state, t, use_stable):
            if state is None:
                return StepView(0.0, 0.0, 0.0, 0.0, 0.0)
            item = catalog[state]
            progress = match_count(item) / 3
            stable = 1.0 if item["color"] == want["color"] else 0.0
            spurious = item["reviews"] / 100
            calibrated = 0.8 * progress + 0.2 * stable
            judge = drifted_judge(calibrated, spurious, t, drift_at)
            steer = progress if policy == "oracle" else (stable if use_stable else judge)
            return StepView(steer, judge, stable, spurious, progress)

        if policy == "confirm" and t < gate.burnin:
            held, _, _ = beam_search(
                None, expand, lambda state, t=t: score_at(state, t, True),
                beam=max(beam, 3), depth=1, rng=random.Random(13_000 + t),
            )
            gate.ref_heur.append(0.0 if held is not None else 1.0)
        use_stable = gate.use_stable(t, policy)
        final, chosen, _ = beam_search(
            None, expand,
            lambda state, t=t, use_stable=use_stable: score_at(state, t, use_stable),
            beam=max(beam, 3), depth=1, rng=rng,
        )
        y = 1.0 if final is not None else 0.0
        loss = 1.0 - y
        successes.append(y)
        losses.append(loss)
        feats.append(
            {
                "judge": chosen.judge,
                "stable": chosen.stable,
                "spurious": chosen.spurious,
                "progress": chosen.progress,
                "y": y,
            }
        )
        if policy in ("gated", "confirm"):
            gate.observe(t, loss)
    return {
        "success": successes,
        "loss": losses,
        "feats": feats,
        "switch_at": gate.switch_at,
        "reverted": gate.reverted,
        "lord_at": gate.lord_at,
    }


# ---------------------------------------------------------------------------
# Localization: permute batch labels, per-feature mean shift
# ---------------------------------------------------------------------------

def localize(feats, drift_at, n_perm=400, seed=0):
    rng = np.random.default_rng(seed)
    keys = ["judge", "stable", "spurious", "progress"]
    pre = feats[:drift_at]
    post = feats[drift_at:]
    if len(pre) < 5 or len(post) < 5:
        return []
    rows = []
    for k in keys:
        a = np.array([f[k] for f in pre], dtype=float)
        b = np.array([f[k] for f in post], dtype=float)
        obs = abs(a.mean() - b.mean())
        pooled = np.concatenate([a, b])
        n_a = len(a)
        hits = 0
        for _ in range(n_perm):
            rng.shuffle(pooled)
            diff = abs(pooled[:n_a].mean() - pooled[n_a:].mean())
            hits += diff >= obs - 1e-15
        p = (1 + hits) / (n_perm + 1)
        rows.append({"feature": k, "abs_mean_shift": float(obs), "p": float(p)})
    rows.sort(key=lambda r: r["p"])
    return rows


def _sse(design: np.ndarray, y: np.ndarray) -> float:
    beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    return float(resid @ resid)


def localize_fsds(feats, drift_at, n_perm=300, seed=0):
    """Univariate concept-drift LOCO.

    For each step feature x, compare y ~ 1 + batch with y ~ 1 + batch + x + batch:x.
    The SSE drop is how much that feature's change across batches predicts success.
    A permutation of x gives a p-value. This is the lightweight FSDS question:
    which feature carries the success drop, not which feature's mean moved.
    """
    keys = ["judge", "stable", "spurious", "progress"]
    n = len(feats)
    if drift_at < 5 or n - drift_at < 5:
        return []
    y = np.array([f["y"] for f in feats], dtype=float)
    batch = np.zeros(n)
    batch[drift_at:] = 1.0
    base = np.column_stack([np.ones(n), batch])
    sse_base = _sse(base, y)
    rng = np.random.default_rng(seed)
    rows = []
    for k in keys:
        x = np.array([f[k] for f in feats], dtype=float)
        sse = _sse(np.column_stack([base, x, batch * x]), y)
        drop = sse_base - sse
        hits = 0
        for _ in range(n_perm):
            xp = rng.permutation(x)
            sse_p = _sse(np.column_stack([base, xp, batch * xp]), y)
            hits += (sse_base - sse_p) >= drop - 1e-12
        rows.append(
            {
                "feature": k,
                "sse_drop": float(drop),
                "p": float((1 + hits) / (n_perm + 1)),
            }
        )
    rows.sort(key=lambda r: (-r["sse_drop"], r["p"]))
    return rows


def rate(xs, a, b):
    sl = xs[a:b]
    if not sl:
        return None
    return float(np.mean(sl))


def null_fdr(n=300, alpha=0.05, reps=20):
    """Stationary losses: how often LORD rejects. Descriptive, not a proof."""
    rng = np.random.default_rng(0)
    counts = []
    for _ in range(reps):
        losses = rng.random(n)
        pvals = empirical_pvals(losses, burnin=20)
        rej = lord_reject(pvals, alpha=alpha)
        counts.append(float(np.mean(rej[20:])))
    return {"mean_reject_rate": float(np.mean(counts)), "alpha": alpha, "reps": reps, "n": n}


def summarize(name, out, drift_at, n):
    suc = out["success"]
    pre = rate(suc, 0, drift_at)
    post = rate(suc, drift_at, n)
    fa = None
    if out["switch_at"] is not None:
        fa = out["switch_at"] <= drift_at
    seed = abs(hash(name)) % 10_000
    loc = localize(out["feats"], drift_at, seed=seed)
    fsds = localize_fsds(out["feats"], drift_at, seed=seed)
    top = loc[0]["feature"] if loc else None
    fsds_top = fsds[0]["feature"] if fsds else None
    return {
        "task": name,
        "pre_success": pre,
        "post_success": post,
        "switch_at": out["switch_at"],
        "false_alarm_before_drift": fa,
        "detection_delay": None
        if out["switch_at"] is None
        else out["switch_at"] - drift_at,
        "lord_at": out.get("lord_at"),
        "top_shifted_feature": top,
        "localization": loc,
        "fsds_top_feature": fsds_top,
        "fsds_localization": fsds,
    }


def main():
    rng = random.Random(SEED)
    n = 40
    drift_at = 16
    policies = ["noisy", "stable", "gated", "confirm", "oracle"]
    task_mix = {
        "game24": 0.28,
        "blocksworld": 0.30,
        "minigrid_doorkey": 0.169,
        "hotpot_twohop": 0.80,
        "webshop_attr": 0.65,
    }
    task_jitter = {
        "game24": 0.0,
        "blocksworld": 0.0,
        "minigrid_doorkey": 0.003,
        "hotpot_twohop": 0.0,
        "webshop_attr": 0.0,
    }
    report = {
        "null_monitor": null_fdr(),
        "drift_mix": task_mix,
        "mix_jitter": task_jitter,
        "tasks": {},
    }

    print("precomputing blocksworld distances...")
    dist_bw = bw_distances()
    print("blocksworld states", len(dist_bw))
    print("precomputing doorkey distances...")
    dist_dk = dk_distances()
    print("doorkey states", len(dist_dk))

    g24 = make_game24(rng, n)
    print("game24 puzzles", len(g24))
    bw_starts = [bw_random_start(rng, dist_bw) for _ in range(n)]

    runners = {
        "game24": lambda policy: run_game24(g24, drift_at, policy, beam=4, rng=random.Random(SEED)),
        "blocksworld": lambda policy: run_blocksworld(
            bw_starts, drift_at, policy, beam=4, rng=random.Random(SEED), dist_map=dist_bw
        ),
        "minigrid_doorkey": lambda policy: run_doorkey(
            n, drift_at, policy, beam=4, rng=random.Random(SEED + 1), dist_map=dist_dk
        ),
        "hotpot_twohop": lambda policy: run_hotpot(
            n, drift_at, policy, beam=3, rng=random.Random(SEED + 2)
        ),
        "webshop_attr": lambda policy: run_webshop(
            n, drift_at, policy, beam=3, rng=random.Random(SEED + 3)
        ),
    }

    global CURRENT_MIX, MIX_JITTER
    for task, fn in runners.items():
        CURRENT_MIX = task_mix[task]
        MIX_JITTER = task_jitter[task]
        report["tasks"][task] = {"mix": CURRENT_MIX, "jitter": MIX_JITTER}
        for policy in policies:
            print(f"running {task} / {policy}")
            out = fn(policy)
            if policy in ("gated", "confirm"):
                summary = summarize(task, out, drift_at, len(out["success"]))
                summary["reverted"] = out.get("reverted")
            elif policy == "noisy":
                seed = abs(hash(task + "-noisy")) % 10_000
                fsds = localize_fsds(out["feats"], drift_at, seed=seed)
                summary = {
                    "pre_success": rate(out["success"], 0, drift_at),
                    "post_success": rate(out["success"], drift_at, len(out["success"])),
                    "fsds_top_feature": fsds[0]["feature"] if fsds else None,
                    "fsds_localization": fsds,
                }
            else:
                summary = {
                    "pre_success": rate(out["success"], 0, drift_at),
                    "post_success": rate(out["success"], drift_at, len(out["success"])),
                }
            report["tasks"][task][policy] = summary
            print(
                f"  pre={summary['pre_success']:.2f} post={summary['post_success']:.2f}"
                + (
                    f" switch={summary.get('switch_at')}"
                    if policy in ("gated", "confirm")
                    else ""
                )
            )

    OUT.write_text(json.dumps(report, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
