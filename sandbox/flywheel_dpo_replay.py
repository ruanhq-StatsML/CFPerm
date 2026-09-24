"""DPO + label pollution + decision-path JSON replay for the forecast flywheel.

Pipeline (as requested)
-----------------------
1. Each flywheel step emits a JSON with a ``decision_path``.
2. Pollute labels on a random region:
     ``y[indices_subset] = y[np.random.permutation(indices_subset)]``
   and/or randomly sample certain periods for contamination.
3. Shuffle the decision_path for replay (same observations, different action order).
4. Feed the flywheel with a light **DPO** preference update:
     chosen = low-MAE clean path, rejected = polluted / shuffled-worse path
     L = -log σ(β (r_chosen − r_rejected))

Stays in sandbox — not AGOD excess/PO.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from sandbox.forecast_bakeoff import Pack, make_supervised
from sandbox.forecast_flywheel import run_flywheel


ACTIONS = ("idle", "log_surprise", "retrain", "fallback_naive")


def sample_period_indices(
    n: int,
    *,
    n_periods: int = 3,
    period_len: int = 40,
    seed: int = 0,
    prefer_range: tuple | None = None,
) -> np.ndarray:
    """Randomly sample contiguous non-overlapping periods; return flat indices.

    If ``prefer_range=(lo, hi)`` is set, periods are drawn inside that window
    (so pollution hits the flywheel eval region, not only warm-up).
    """
    rng = np.random.default_rng(seed)
    if n < period_len + 1:
        return np.arange(n)
    if prefer_range is not None:
        lo, hi = int(prefer_range[0]), int(prefer_range[1])
        lo = max(0, lo)
        hi = min(n, hi)
        if hi - lo < period_len:
            lo, hi = 0, n
        max_start = max(hi - period_len, lo)
        grid = list(range(lo, max_start + 1, max(period_len // 2, 1)))
    else:
        max_start = n - period_len
        grid = list(range(0, max_start + 1, period_len))
    if not grid:
        grid = [0]
    n_take = min(n_periods, len(grid))
    starts = rng.choice(grid, size=n_take, replace=False)
    idx = []
    for s in starts:
        idx.extend(range(int(s), int(min(s + period_len, n))))
    return np.unique(np.asarray(idx, dtype=int))


def pollute_labels_region(
    y: np.ndarray,
    indices_subset: np.ndarray,
    *,
    seed: int = 0,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Shuffle labels inside ``indices_subset`` (region pollution).

    ``y_out[indices_subset] = y[np.random.permutation(indices_subset)]``
    """
    y = np.asarray(y, dtype=float).copy()
    idx = np.asarray(indices_subset, dtype=int)
    idx = idx[(idx >= 0) & (idx < len(y))]
    meta = {"n_polluted": 0, "indices_head": []}
    if idx.size < 2:
        return y, meta
    rng = np.random.default_rng(seed)
    perm = rng.permutation(idx)
    y[idx] = y[perm]
    meta = {
        "n_polluted": int(idx.size),
        "indices_head": idx[:12].tolist(),
        "method": "np.random.permutation(indices_subset) label shuffle",
    }
    return y, meta


def tick_to_step_json(
    tick: Any,
    *,
    step: int,
    decision_path: Sequence[str],
    polluted: bool = False,
) -> Dict[str, Any]:
    """One step JSON for logging / replay."""
    return {
        "step": int(step),
        "t": int(tick.t),
        "y_true": float(tick.y_true),
        "y_hat": float(tick.y_hat),
        "residual": float(tick.residual),
        "surprise": float(tick.surprise),
        "action": str(tick.action),
        "model": str(tick.model),
        "decision_path": list(decision_path),
        "polluted": bool(polluted),
    }


def build_decision_path(action: str, model: str, surprise: float) -> List[str]:
    """Canonical path for one step (nodes an agent traversed)."""
    path = ["observe", f"forecast:{model}"]
    if surprise >= 2.5:
        path.append("flag_surprise")
    path.append(f"decide:{action}")
    if action == "retrain":
        path.append("refit")
    if action == "fallback_naive":
        path.append("switch:naive_last")
    path.append("log")
    return path


def emit_step_jsonl(
    flywheel_run: Dict[str, Any],
    ticks: Sequence[Any],
    *,
    polluted: bool = False,
) -> List[Dict[str, Any]]:
    steps = []
    for i, tick in enumerate(ticks):
        path = build_decision_path(tick.action, tick.model, tick.surprise)
        steps.append(
            tick_to_step_json(tick, step=i, decision_path=path, polluted=polluted)
        )
    return steps


def shuffle_decision_path(
    path: Sequence[str],
    *,
    seed: int,
    keep_ends: bool = True,
) -> List[str]:
    """Shuffle middle nodes of a decision path for replay."""
    p = list(path)
    if len(p) <= 2:
        return p
    rng = np.random.default_rng(seed)
    if keep_ends:
        mid = p[1:-1]
        rng.shuffle(mid)
        return [p[0]] + mid + [p[-1]]
    rng.shuffle(p)
    return p


def replay_with_shuffled_paths(
    steps: Sequence[Dict[str, Any]],
    *,
    seed: int = 0,
) -> List[Dict[str, Any]]:
    """Replay: copy step JSONs but shuffle each decision_path."""
    out = []
    for i, s in enumerate(steps):
        sp = dict(s)
        sp["decision_path_orig"] = list(s.get("decision_path") or [])
        sp["decision_path"] = shuffle_decision_path(
            s.get("decision_path") or [], seed=seed + i
        )
        sp["replay"] = True
        out.append(sp)
    return out


def path_reward(steps: Sequence[Dict[str, Any]]) -> float:
    """Scalar path score for DPO: prefer low |residual| and fewer surprises."""
    if not steps:
        return 0.0
    mae = float(np.mean([abs(float(s["residual"])) for s in steps]))
    sur = float(np.mean([float(s["surprise"]) for s in steps]))
    # higher reward is better
    return float(-mae - 0.1 * sur)


def dpo_preference_loss(
    r_chosen: float,
    r_rejected: float,
    *,
    beta: float = 1.0,
) -> float:
    """Bradley–Terry / DPO-style loss: -log σ(β (r_c − r_r))."""
    z = beta * (r_chosen - r_rejected)
    # stable softplus form: -log σ(z) = softplus(-z)
    return float(np.logaddexp(0.0, -z))


def dpo_update_from_pairs(
    pairs: Sequence[Tuple[float, float]],
    *,
    beta: float = 1.0,
) -> Dict[str, float]:
    """Aggregate DPO loss over preference pairs (chosen, rejected) rewards."""
    if not pairs:
        return {"n_pairs": 0.0, "mean_loss": float("nan"), "mean_margin": float("nan")}
    losses = [dpo_preference_loss(c, r, beta=beta) for c, r in pairs]
    margins = [c - r for c, r in pairs]
    return {
        "n_pairs": float(len(pairs)),
        "mean_loss": float(np.mean(losses)),
        "mean_margin": float(np.mean(margins)),
        "frac_chosen_better": float(np.mean([m > 0 for m in margins])),
    }


def run_pollute_replay_dpo(
    pack: Pack,
    *,
    model_name: str = "hgb",
    seed: int = 0,
    n_periods: int = 3,
    period_len: int = 50,
    warm: int = 180,
    max_steps: int = 200,
    out_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Full demo: clean flywheel → polluted region → path shuffle replay → DPO."""
    Z, y = make_supervised(pack.X, pack.y, n_lags=5)
    # Pollute periods that overlap the flywheel eval window (after warm).
    # Supervised index i ≈ original index i+n_lags; target y[i] is original[i+n_lags].
    n_lags = 5
    eval_lo = warm + n_lags
    eval_hi = warm + max_steps + n_lags
    idx = sample_period_indices(
        len(pack.y),
        n_periods=n_periods,
        period_len=period_len,
        seed=seed,
        prefer_range=(eval_lo, min(len(pack.y), eval_hi)),
    )
    y_dirty, poll_meta = pollute_labels_region(pack.y, idx, seed=seed + 1)
    pack_clean = pack
    pack_dirty = Pack(name=pack.name + "_polluted", X=pack.X, y=y_dirty)

    clean = run_flywheel(
        pack_clean,
        model_name=model_name,
        warm=warm,
        max_steps=max_steps,
        seed=seed,
    )
    dirty = run_flywheel(
        pack_dirty,
        model_name=model_name,
        warm=warm,
        max_steps=max_steps,
        seed=seed,
    )
    # recover ticks from a second internal run — run_flywheel doesn't return ticks.
    # Re-run thin collectors:
    clean_steps = _collect_steps(pack_clean, model_name, warm, max_steps, seed, polluted=False)
    dirty_steps = _collect_steps(pack_dirty, model_name, warm, max_steps, seed, polluted=True)
    replay_steps = replay_with_shuffled_paths(clean_steps, seed=seed + 7)

    r_clean = path_reward(clean_steps)
    r_dirty = path_reward(dirty_steps)
    r_replay = path_reward(replay_steps)
    # DPO pairs: prefer clean over polluted; prefer clean over shuffled-path replay
    # (replay keeps residuals — reward same as clean — so use path-structure penalty)
    r_replay_penalized = r_replay - 0.05 * _path_shuffle_distance(clean_steps, replay_steps)
    pairs = [(r_clean, r_dirty), (r_clean, r_replay_penalized)]
    dpo = dpo_update_from_pairs(pairs, beta=1.0)

    blob = {
        "dataset": pack.name,
        "model": model_name,
        "pollution": poll_meta,
        "sampled_periods": {
            "n_periods": n_periods,
            "period_len": period_len,
            "n_indices": int(len(idx)),
        },
        "rewards": {
            "clean": r_clean,
            "polluted": r_dirty,
            "replay_shuffled_path": r_replay_penalized,
        },
        "flywheel_clean": {k: clean[k] for k in clean if k != "note"},
        "flywheel_polluted": {k: dirty[k] for k in dirty if k != "note"},
        "dpo": dpo,
        "n_step_json": len(clean_steps),
        "step_json_head": clean_steps[:3],
        "replay_path_example": {
            "orig": (clean_steps[0].get("decision_path") if clean_steps else []),
            "shuffled": (replay_steps[0].get("decision_path") if replay_steps else []),
        },
        "reading": _reading(dpo, r_clean, r_dirty),
    }
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "steps_clean.jsonl").write_text(
            "\n".join(json.dumps(s) for s in clean_steps) + "\n"
        )
        (out_dir / "steps_polluted.jsonl").write_text(
            "\n".join(json.dumps(s) for s in dirty_steps) + "\n"
        )
        (out_dir / "steps_replay_shuffled.jsonl").write_text(
            "\n".join(json.dumps(s) for s in replay_steps) + "\n"
        )
        (out_dir / "dpo_pollute_replay.json").write_text(
            json.dumps(blob, indent=2) + "\n"
        )
    return blob


def _collect_steps(
    pack: Pack,
    model_name: str,
    warm: int,
    max_steps: int,
    seed: int,
    *,
    polluted: bool,
) -> List[Dict[str, Any]]:
    # Import private pieces by re-running flywheel with tick capture via monkey patch style:
    # duplicate minimal loop using public run_flywheel internals — call run_flywheel and
    # rebuild paths from returned aggregates is lossy. Instead instrument via copy of loop.
    from sandbox.forecast_flywheel import _ModelBox

    Z, y = make_supervised(pack.X, pack.y, n_lags=5)
    if len(y) < warm + 30:
        return []
    tr_end = warm
    scale = float(np.std(y[:warm]) + 1e-6)
    fallback = False
    steps_n = min(max_steps, len(y) - warm)
    box = _ModelBox(model_name, n_lags=5, seed=seed)
    box.fit(Z[:tr_end], y[:tr_end])
    naive_box = _ModelBox("naive_last", n_lags=5, seed=seed)
    naive_box.fit(Z[:tr_end], y[:tr_end])
    lift_window: List[float] = []
    out: List[Dict[str, Any]] = []
    n_retrain = 0
    for step in range(steps_n):
        i = warm + step
        use = "naive_last" if fallback else model_name
        predictor = naive_box if fallback else box
        y_hat = predictor.predict_one(Z[i])
        resid = float(y[i] - y_hat)
        surprise = abs(resid) / scale
        action = "idle"
        if surprise >= 2.5:
            action = "log_surprise"
        if (step + 1) % 50 == 0 and not fallback:
            action = "retrain"
            tr_end = i
            box.fit(Z[:tr_end], y[:tr_end])
            n_retrain += 1
            scale = float(np.std(y[:tr_end]) + 1e-6)
        naive_hat = naive_box.predict_one(Z[i])
        lift_window.append(abs(y[i] - naive_hat) - abs(y[i] - y_hat))
        if len(lift_window) >= 40:
            if float(np.mean(lift_window[-40:])) < 0 and not fallback and model_name == "hgb":
                action = "fallback_naive"
                fallback = True
        path = build_decision_path(action, use, surprise)
        out.append(
            {
                "step": step,
                "t": int(i),
                "y_true": float(y[i]),
                "y_hat": float(y_hat),
                "residual": resid,
                "surprise": float(surprise),
                "action": action,
                "model": use,
                "decision_path": path,
                "polluted": polluted,
            }
        )
    return out


def _path_shuffle_distance(
    a: Sequence[Dict[str, Any]], b: Sequence[Dict[str, Any]]
) -> float:
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    diffs = 0
    for i in range(n):
        if list(a[i].get("decision_path") or []) != list(b[i].get("decision_path") or []):
            diffs += 1
    return float(diffs / n)


def _reading(dpo: Dict[str, float], r_c: float, r_d: float) -> str:
    if not np.isfinite(dpo.get("mean_loss", np.nan)):
        return "no DPO pairs"
    if r_c > r_d and dpo.get("frac_chosen_better", 0) >= 0.5:
        return (
            "DPO prefers clean path over polluted-label path "
            f"(margin={dpo.get('mean_margin'):.4g}, loss={dpo.get('mean_loss'):.4g})"
        )
    return (
        "pollution did not clearly hurt path reward — check period size / model; "
        f"margin={dpo.get('mean_margin')}"
    )
