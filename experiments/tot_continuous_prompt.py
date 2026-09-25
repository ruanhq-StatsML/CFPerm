"""Continuous-time prompt on a ToT trajectory.

Clock T is one beat per committed search step, across episodes.
X is the step's scores. Y is the episode outcome, attached only after the
episode ends. The prompt for the current beat contains X and no Y.
The forest, the hold-out stack weight, and the feature split are fit on the
first 12 episodes and then frozen. Later beats only slide the window.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

import tot_agent_probe as probe
from online_llm_stack import (
    compute_vimp,
    holdout_stack_weight,
    hop_fires,
    make_rf_adapter,
    make_tabpfn_adapter,
    online_stack_weights,
    pval_vs_ref,
    render_sliding_prompt,
    split_by_vimp,
)

OUT = Path(__file__).resolve().parent / "results_prompt_prototype.json"
NAMES = ["judge", "stable", "spurious", "progress", "gap", "depth_frac"]
REF_EPISODES = 12
PROMPT_WINDOW = 8
SEED = 2026


def step_x(view, step: int, depth: int) -> np.ndarray:
    return np.array(
        [
            float(view.judge),
            float(view.stable),
            float(view.spurious),
            float(view.progress),
            abs(float(view.judge) - float(view.stable)),
            float(step) / float(max(depth, 1)),
        ],
        dtype=float,
    )


def collect_doorkey(n: int, drift_at: int):
    probe.CURRENT_MIX = 0.169
    probe.MIX_JITTER = 0.003
    dist = probe.dk_distances()
    start = (1, 1, 0, 0, 0)
    depth = 12

    def expand(state):
        return probe.dk_neighbors(state)

    def score_at(state, t, use_stable=False):
        d = dist.get(state, 30)
        progress = max(0.0, 1.0 - d / 20)
        stable = probe.dk_stable(state)
        spurious = probe.dk_spurious(state)
        calibrated = 0.85 * progress + 0.15 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, drift_at)
        steer = stable if use_stable else judge
        return probe.StepView(steer, judge, stable, spurious, progress)

    episodes = []
    for t in range(n):
        trace = []
        final, _, _ = probe.beam_search(
            start, expand, lambda state, t=t: score_at(state, t),
            beam=4, depth=depth, rng=random.Random(SEED + t), trace=trace,
        )
        xs = [step_x(view, step, depth) for step, view in trace]
        episodes.append((np.vstack(xs) if xs else np.zeros((0, len(NAMES))), 1.0 if final is not None else 0.0))
    return episodes


def collect_hotpot(n: int, drift_at: int):
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    bank = probe.hop_questions()
    depth = 2
    episodes = []
    for t in range(n):
        question, gold = bank[t % len(bank)]
        need = 2

        def done_fn(picked, gold=gold, need=need):
            return len(picked) >= 2 and len(picked & gold) >= need

        def expand(picked, gold=gold):
            if len(picked) >= 2:
                return []
            out = []
            for pid in range(len(probe.PASSAGES)):
                if pid in picked:
                    continue
                nxt = frozenset(set(picked) | {pid})
                out.append((nxt, done_fn(nxt)))
            return out

        def score_at(picked, t=t, question=question, gold=gold, need=need):
            if not picked:
                progress = stable = spurious = 0.0
            else:
                progress = len(set(picked) & gold) / need
                stable = float(np.mean([probe.hop_overlap(question, pid) for pid in picked]))
                spurious = float(np.mean([len(probe.PASSAGES[pid].split()) / 12 for pid in picked]))
            calibrated = 0.75 * progress + 0.25 * stable
            judge = probe.drifted_judge(calibrated, spurious, t, drift_at)
            return probe.StepView(judge, judge, stable, spurious, progress)

        trace = []
        final, _, _ = probe.beam_search(
            frozenset(), expand, score_at, beam=3, depth=depth, rng=random.Random(SEED + 1000 + t), trace=trace,
        )
        xs = [step_x(view, step, depth) for step, view in trace]
        episodes.append((np.vstack(xs) if xs else np.zeros((0, len(NAMES))), 1.0 if final is not None else 0.0))
    return episodes


def _rows(episodes):
    xs, ys, ep = [], [], []
    for i, (X, y) in enumerate(episodes):
        if len(X) == 0:
            continue
        xs.append(X)
        ys.append(np.full(len(X), y))
        ep.append(np.full(len(X), i))
    return np.vstack(xs), np.concatenate(ys), np.concatenate(ep)


def run_stream(name: str, episodes, drift_at: int):
    X, Y, ep = _rows(episodes)
    ref = ep < REF_EPISODES
    X_ref, Y_ref = X[ref], Y[ref]
    vimp = compute_vimp(X_ref, Y_ref, seed=SEED)
    if float(np.max(np.abs(vimp))) < 1e-8:
        high, low = list(range(len(NAMES))), []
    else:
        high, low = split_by_vimp(vimp, n_high=2)
    rf = make_rf_adapter()
    llm = make_tabpfn_adapter(n_context=8, window=64, seed=SEED, temperature=1.0)
    fit_rf = rf["fit"](X_ref[:, high], Y_ref, seed=SEED)
    y_rf_ref = np.asarray(rf["predict"](fit_rf, X_ref[:, high]), float)
    kind, obj = fit_rf
    y_null = y_rf_ref
    if kind == "rf" and hasattr(obj, "oob_prediction_"):
        oob = np.asarray(obj.oob_prediction_, float).ravel()
        if len(oob) == len(Y_ref) and np.isfinite(oob).mean() > 0.5:
            y_null = np.where(np.isfinite(oob), oob, y_rf_ref)
    e_ref = Y_ref - y_rf_ref
    n_hold = max(int(0.2 * len(Y_ref)), 8)
    n_fit = max(len(Y_ref) - n_hold, 1)
    if low:
        fit_llm = llm["fit"](X_ref[:n_fit][:, low], e_ref[:n_fit], seed=SEED)
        corr_hold, sig_hold = llm["predict_uncertainty"](fit_llm, X_ref[n_fit:][:, low])
        w_stack = holdout_stack_weight(e_ref[n_fit:], corr_hold)
        finite = sig_hold[np.isfinite(sig_hold)]
        sigma_med = float(np.median(finite)) if len(finite) else float(np.std(e_ref) + 1e-8)
    else:
        fit_llm = None
        w_stack = 0.0
        sigma_med = float(np.std(e_ref) + 1e-8)
    oob_mse = float(np.mean((Y_ref - y_null) ** 2))
    ref_batch = []
    for i in range(0, len(Y_ref), 8):
        sl = slice(i, i + 8)
        if len(Y_ref[sl]) == 0:
            break
        ref_batch.append(float(np.mean((Y_ref[sl] - y_null[sl]) ** 2)))

    window = []
    prompts = []
    pred_rows = []
    e_prev = None
    hops = []
    # Reference beats enter the window only after that episode's Y is known.
    for i, (Xe, ye) in enumerate(episodes[:REF_EPISODES]):
        if len(Xe) == 0:
            continue
        y_hat_rf = np.asarray(rf["predict"](fit_rf, Xe[:, high]), float)
        if fit_llm is None:
            w = np.zeros(len(Xe))
        else:
            _, sig = llm["predict_uncertainty"](fit_llm, Xe[:, low])
            w = online_stack_weights(sig, w_stack=w_stack, sigma_med=sigma_med)
        for s in range(len(Xe)):
            window.append({
                "t": int(i * 100 + s),
                "y": float(ye),
                "rf": float(y_hat_rf[s]),
                "e": float(ye - y_hat_rf[s]),
                "w": float(w[s]),
                "x": np.round(Xe[s], 4).tolist(),
            })
        window = window[-PROMPT_WINDOW:]

    beat = REF_EPISODES * 100
    for i, (Xe, ye) in enumerate(episodes[REF_EPISODES:], start=REF_EPISODES):
        if len(Xe) == 0:
            continue
        step_pred = []
        for s in range(len(Xe)):
            y_rf = float(np.asarray(rf["predict"](fit_rf, Xe[s:s + 1, high]), float)[0])
            if fit_llm is None:
                w, corr0 = 0.0, 0.0
            else:
                corr, sig = llm["predict_uncertainty"](fit_llm, Xe[s:s + 1, low])
                w = float(online_stack_weights(sig, w_stack=w_stack, sigma_med=sigma_med)[0])
                corr0 = float(corr[0])
            y_hat = y_rf + w * corr0
            prompt = render_sliding_prompt(window, Xe[s], w_stack=w_stack, sigma_med=sigma_med)
            prompt = "X columns: " + ",".join(NAMES) + "\n" + "high_vimp=" + ",".join(NAMES[j] for j in high) + " low_vimp=" + ",".join(NAMES[j] for j in low) + "\n" + prompt
            prompts.append({"episode": i, "step": s, "text": prompt, "y_hat": y_hat, "y_rf": y_rf})
            step_pred.append(y_hat)
            beat += 1
        mse = float(np.mean((ye - np.asarray(step_pred)) ** 2))
        hops.append(bool(hop_fires(mse, e_prev, e_floor=0.02)))
        e_prev = mse
        pred_rows.append({"episode": i, "y": ye, "mse": mse, "y_hat_last": step_pred[-1], "post": i >= drift_at})
        y_rf = np.asarray(rf["predict"](fit_rf, Xe[:, high]), float)
        if fit_llm is None:
            w = np.zeros(len(Xe))
        else:
            _, sig = llm["predict_uncertainty"](fit_llm, Xe[:, low])
            w = online_stack_weights(sig, w_stack=w_stack, sigma_med=sigma_med)
        for s in range(len(Xe)):
            window.append({
                "t": int(beat - len(Xe) + s),
                "y": float(ye),
                "rf": float(y_rf[s]),
                "e": float(ye - y_rf[s]),
                "w": float(w[s]),
                "x": np.round(Xe[s], 4).tolist(),
            })
        window = window[-PROMPT_WINDOW:]

    mses = np.array([r["mse"] for r in pred_rows], float)
    post = np.array([r["post"] for r in pred_rows])
    return {
        "task": name,
        "T": "one committed search step. The prompt index t is that beat. T is not a column of X.",
        "X": NAMES,
        "Y": "episode success, 1 if the search returns a solved state else 0. Written onto every step of that episode only after the episode ends. The query beat does not contain Y.",
        "feature_selection": {
            "vimp": {NAMES[j]: float(vimp[j]) for j in range(len(NAMES))},
            "rf_columns": [NAMES[j] for j in high],
            "residual_columns": [NAMES[j] for j in low],
            "note": (
                "reference Y does not vary, so every column has vimp 0 and no column is sent to the residual"
                if float(np.max(np.abs(vimp))) < 1e-8
                else "high-vimp columns stay in the frozen forest; low-vimp columns go to the residual"
            ),
        },
        "w_stack": float(w_stack),
        "sigma_med": float(sigma_med),
        "oob_mse": oob_mse,
        "trail_mse": float(np.mean(mses)) if len(mses) else None,
        "pval_vs_oob_batches": pval_vs_ref(mses[:: max(len(mses) // 4, 1)] if len(mses) else mses, np.asarray(ref_batch)).tolist(),
        "hop_rate_post": float(np.mean([h for h, r in zip(hops, pred_rows) if r["post"]])) if any(post) else None,
        "post_abs_err_last_step": float(np.mean([abs(r["y"] - r["y_hat_last"]) for r in pred_rows if r["post"]])) if any(post) else None,
        "prompt_example": prompts[0]["text"] if prompts else "",
        "n_prompts": len(prompts),
    }


def main():
    n, drift_at = 40, 16
    report = {
        "doorkey": run_stream("minigrid_doorkey", collect_doorkey(n, drift_at), drift_at),
        "hotpot": run_stream("hotpot_twohop", collect_hotpot(n, drift_at), drift_at),
    }
    OUT.write_text(json.dumps(report, indent=2))
    for key, block in report.items():
        print(key, "w_stack", round(block["w_stack"], 3), "oob", round(block["oob_mse"], 3), "trail", None if block["trail_mse"] is None else round(block["trail_mse"], 3))
        print(" rf", block["feature_selection"]["rf_columns"], "residual", block["feature_selection"]["residual_columns"])
        print("--- prompt ---")
        print(block["prompt_example"])
        print()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
