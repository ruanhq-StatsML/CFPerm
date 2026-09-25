"""Layered ToT features, then the concept-drift and mean-shift tests.

L4 is a fixed encoding of the committed state, not a 768-d text encoder.
Residual of the current episode is not a column: that Y is the target.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import tot_agent_probe as probe

OUT = Path(__file__).resolve().parent / "results_fsds_tot_layers.json"
SEED = 2026
N = 40
DRIFT_AT = 16
REF = 12

LAYERS = {
    "L1_graph": ["depth", "n_children", "n_siblings", "path_len", "branch"],
    "L2_history": ["prev_score", "score_diff", "roll_mean", "roll_std", "prev_y"],
    "L3_stat": ["score_std", "score_z", "score_entropy"],
    "L4_embed": ["e0", "e1", "e2", "e3"],
    "L5_interact": ["child_per_depth", "depth_x_gap", "score_x_len"],
    "L6_context": ["beam", "max_depth", "step_frac"],
}


def _entropy(xs: list[float]) -> float:
    a = np.abs(np.asarray(xs, float))
    s = float(a.sum())
    if s <= 1e-12:
        return 0.0
    p = a / s
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def search_log(root, expand, score_fn, beam: int, depth: int, rng):
    frontier = [root]
    trace = []
    expansions = 0
    for step in range(depth):
        cand = []
        for state in frontier:
            kids = expand(state)
            n_sib = max(len(kids) - 1, 0)
            for nxt, done in kids:
                expansions += 1
                view = score_fn(nxt)
                cand.append((view.value + rng.uniform(-1e-9, 1e-9), nxt, done, view, len(kids), n_sib, step))
        if not cand:
            break
        cand.sort(key=lambda z: z[0], reverse=True)
        top = cand[0]
        second = cand[1][0] if len(cand) > 1 else top[0]
        gap = float(top[0] - second)
        hit = None
        for item in cand[:beam]:
            if item[2]:
                hit = item
                break
        chosen = hit if hit is not None else top
        _, nxt, done, view, n_ch, n_sib, step_i = chosen
        trace.append({
            "step": step_i,
            "depth": step_i + 1,
            "n_children": n_ch,
            "n_siblings": n_sib,
            "gap": gap,
            "score": float(view.value),
            "state": nxt,
            "beam": beam,
            "max_depth": depth,
        })
        if done:
            return nxt, trace, expansions
        frontier = [item[1] for item in cand[:beam]]
    return None, trace, expansions


def episode_row(trace, y: float, prev_y: float, prev_score: float) -> dict:
    scores = [r["score"] for r in trace] or [0.0]
    depths = [r["depth"] for r in trace] or [0]
    children = [r["n_children"] for r in trace] or [0]
    last = trace[-1]
    state = last["state"]
    if isinstance(state, tuple):
        vals = sorted(float(v) for v in state)[:4]
        vals = vals + [0.0] * (4 - len(vals))
        emb = [v / 24.0 for v in vals]
    else:
        emb = [0.0, 0.0, 0.0, 0.0]
        if isinstance(state, int):
            emb[int(state) % 4] = 1.0
    sd = float(np.std(scores))
    mu = float(np.mean(scores))
    z = (scores[-1] - mu) / sd if sd > 1e-8 else 0.0
    depth_m = float(np.mean(depths))
    child_m = float(np.mean(children))
    gap_m = float(np.mean([r["gap"] for r in trace])) if trace else 0.0
    path_len = float(len(trace))
    return {
        "depth": depth_m,
        "n_children": child_m,
        "n_siblings": float(np.mean([r["n_siblings"] for r in trace])) if trace else 0.0,
        "path_len": path_len,
        "branch": child_m,
        "prev_score": float(prev_score),
        "score_diff": float(scores[-1] - scores[0]) if scores else 0.0,
        "roll_mean": mu,
        "roll_std": sd,
        "prev_y": float(prev_y),
        "score_std": sd,
        "score_z": float(z),
        "score_entropy": _entropy(scores),
        "e0": emb[0], "e1": emb[1], "e2": emb[2], "e3": emb[3],
        "child_per_depth": child_m / max(depth_m, 1.0),
        "depth_x_gap": depth_m * gap_m,
        "score_x_len": mu * path_len,
        "beam": float(last["beam"]),
        "max_depth": float(last["max_depth"]),
        "step_frac": float(last["step"]) / max(float(last["max_depth"]), 1.0),
        "y": float(y),
    }


def game24_rows():
    probe.CURRENT_MIX = 0.28
    probe.MIX_JITTER = 0.0
    puzzles = probe.make_game24(probe.random.Random(SEED), N)

    def score_at(state, t):
        solv = 1.0 if probe._g24_solvable(state) else 0.0
        close = max(probe.math.exp(-abs(x - 24) / 8) for x in state)
        stable = close if len(state) == 1 else 0.45 * close + 0.15
        spurious = sum(state) / (13 * len(state))
        calibrated = 0.85 * solv + 0.15 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
        return probe.StepView(judge, judge, stable, spurious, solv)

    rows = []
    prev_y, prev_score = 0.0, 0.0
    for t, puzzle in enumerate(puzzles):
        final, trace, _ = search_log(
            puzzle, probe._g24_children, lambda s, t=t: score_at(s, t),
            beam=4, depth=3, rng=probe.random.Random(SEED + 3000 + t),
        )
        y = 1.0 if final is not None else 0.0
        row = episode_row(trace, y, prev_y, prev_score)
        row["t"] = t
        rows.append(row)
        prev_y, prev_score = y, row["roll_mean"]
    return rows


def offer_rows():
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    offers = (("habit", 0.50, 0.15, 0.90), ("growth", 0.67, 0.12, 0.55), ("click", 0.20, 0.98, 0.10))

    def activates(i, t):
        name = offers[i][0]
        if name == "habit":
            return t % 2 == 0
        if name == "growth":
            return t % 3 != 0
        return t % 5 == 0

    def expand(_s, t):
        return [(i, activates(i, t)) for i in range(3)]

    def score_at(i, t):
        _, stable, spurious, progress = offers[i]
        calibrated = 0.8 * progress + 0.2 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
        return probe.StepView(judge, judge, stable, spurious, progress)

    rows = []
    prev_y, prev_score = 0.0, 0.0
    for t in range(N):
        final, trace, _ = search_log(
            None, lambda s, t=t: expand(s, t), lambda s, t=t: score_at(s, t),
            beam=1, depth=1, rng=probe.random.Random(SEED + 7000 + t),
        )
        y = 1.0 if final is not None else 0.0
        row = episode_row(trace, y, prev_y, prev_score)
        row["t"] = t
        rows.append(row)
        prev_y, prev_score = y, row["roll_mean"]
    return rows


def _loco(feats, keys):
    y = np.array([f["y"] for f in feats], float)
    batch = np.zeros(len(feats))
    batch[REF:] = 1.0
    base = np.column_stack([np.ones(len(feats)), batch])
    sse_base = probe._sse(base, y)
    rng = np.random.default_rng(SEED)
    rows = []
    for k in keys:
        x = np.array([f[k] for f in feats], float)
        drop = sse_base - probe._sse(np.column_stack([base, x, batch * x]), y)
        hits = 0
        for _ in range(80):
            xp = rng.permutation(x)
            hits += (sse_base - probe._sse(np.column_stack([base, xp, batch * xp]), y)) >= drop - 1e-12
        rows.append({"feature": k, "sse_drop": float(drop), "p": float((1 + hits) / 81)})
    rows.sort(key=lambda r: (-r["sse_drop"], r["p"]))
    return rows


def _shift(feats, keys):
    rng = np.random.default_rng(SEED)
    pre, post = feats[:DRIFT_AT], feats[DRIFT_AT:]
    rows = []
    for k in keys:
        a = np.array([f[k] for f in pre], float)
        b = np.array([f[k] for f in post], float)
        obs = abs(float(a.mean() - b.mean()))
        pooled = np.concatenate([a, b])
        hits = 0
        for _ in range(80):
            rng.shuffle(pooled)
            diff = abs(float(pooled[: len(a)].mean() - pooled[len(a) :].mean()))
            hits += diff >= obs - 1e-15
        rows.append({"feature": k, "abs_mean_shift": obs, "p": float((1 + hits) / 81)})
    return rows


def layer_report(rows: list[dict]) -> dict:
    keys = [k for ks in LAYERS.values() for k in ks]
    feats = [{k: r[k] for k in keys} | {"y": r["y"]} for r in rows]
    concept = _loco(feats, keys)
    covariate = _shift(feats, keys)
    by_c = {r["feature"]: r for r in concept}
    by_v = {r["feature"]: r for r in covariate}
    layers = []
    for name, cols in LAYERS.items():
        drops = [max(by_c[c]["sse_drop"], 0.0) for c in cols]
        shifts = [by_v[c]["abs_mean_shift"] for c in cols]
        best = max(cols, key=lambda c: by_c[c]["sse_drop"])
        layers.append({
            "layer": name,
            "sse_sum": float(sum(drops)),
            "sse_best": float(by_c[best]["sse_drop"]),
            "best_feature": best,
            "best_p": float(by_c[best]["p"]),
            "mean_shift_max": float(max(shifts)),
            "shift_feature": max(cols, key=lambda c: by_v[c]["abs_mean_shift"]),
        })
    layers.sort(key=lambda r: -r["sse_sum"])
    total = sum(r["sse_sum"] for r in layers) or 1.0
    for r in layers:
        r["sse_share"] = r["sse_sum"] / total
    y = np.array([r["y"] for r in rows])
    return {
        "pre_Y": float(y[:DRIFT_AT].mean()),
        "post_Y": float(y[DRIFT_AT:].mean()),
        "layers": layers,
        "top_features": concept[:8],
    }


def main():
    report = {
        "game24": layer_report(game24_rows()),
        "offer": layer_report(offer_rows()),
        "note": "L4 encodes the committed state in four numbers. It is not a 768-d text embedding. Current-episode residual is excluded.",
    }
    OUT.write_text(json.dumps(report, indent=2))
    for task, block in report.items():
        if task == "note":
            continue
        print(task, "Y", round(block["pre_Y"], 3), "->", round(block["post_Y"], 3))
        for layer in block["layers"]:
            print(
                f"  {layer['layer']:12} share={layer['sse_share']:.3f} best={layer['best_feature']} "
                f"sse={layer['sse_best']:.3f} p={layer['best_p']:.3f} shift={layer['mean_shift_max']:.3f}"
            )
    print("wrote", OUT)


if __name__ == "__main__":
    main()
