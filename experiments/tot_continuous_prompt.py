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
    clever_covariate,
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
        final, _, n_children = probe.beam_search(
            start, expand, lambda state, t=t: score_at(state, t),
            beam=4, depth=depth, rng=random.Random(SEED + t), trace=trace,
        )
        xs = [step_x(view, step, depth) for step, view in trace]
        episodes.append((np.vstack(xs) if xs else np.zeros((0, len(NAMES))), 1.0 if final is not None else 0.0, int(n_children)))
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
        final, _, n_children = probe.beam_search(
            frozenset(), expand, score_at, beam=3, depth=depth, rng=random.Random(SEED + 1000 + t), trace=trace,
        )
        xs = [step_x(view, step, depth) for step, view in trace]
        episodes.append((np.vstack(xs) if xs else np.zeros((0, len(NAMES))), 1.0 if final is not None else 0.0, int(n_children)))
    return episodes


def _rows(episodes):
    xs, ys, ep = [], [], []
    for i, item in enumerate(episodes):
        X, y = item[0], item[1]
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
    anomaly_rows = []
    e_prev = None
    hops = []
    # Reference beats enter the window only after that episode's Y is known.
    children = [int(item[2]) for item in episodes]
    for i, item in enumerate(episodes[:REF_EPISODES]):
        Xe, ye = item[0], item[1]
        if len(Xe) == 0:
            continue
        y_hat_rf = np.asarray(rf["predict"](fit_rf, Xe[:, high]), float)
        if fit_llm is None:
            w = np.zeros(len(Xe))
        else:
            _, sig = llm["predict_uncertainty"](fit_llm, Xe[:, low])
            w = online_stack_weights(sig, w_stack=w_stack, sigma_med=sigma_med)
        for s in range(len(Xe)):
            row = {
                "t": int(i * 100 + s),
                "y": float(ye),
                "rf": float(y_hat_rf[s]),
                "e": float(ye - y_hat_rf[s]),
                "w": float(w[s]),
                "x": np.round(Xe[s], 4).tolist(),
            }
            row.update(clever_covariate(ye, float(Xe[s, 1])))
            window.append(row)
            anomaly_rows.append({"episode": i, "post": False, "y": float(ye), "H": row["H"], "anomaly": row["anomaly"], "adjustment": row["adjustment"], "prior": row["prior"], "residual": row["e"]})
        window = window[-PROMPT_WINDOW:]

    beat = REF_EPISODES * 100
    for i, item in enumerate(episodes[REF_EPISODES:], start=REF_EPISODES):
        Xe, ye = item[0], item[1]
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
            prompt = render_sliding_prompt(
                window, Xe[s], w_stack=w_stack, sigma_med=sigma_med, query_prior=float(Xe[s, 1])
            )
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
            row = {
                "t": int(beat - len(Xe) + s),
                "y": float(ye),
                "rf": float(y_rf[s]),
                "e": float(ye - y_rf[s]),
                "w": float(w[s]),
                "x": np.round(Xe[s], 4).tolist(),
            }
            row.update(clever_covariate(ye, float(Xe[s, 1])))
            window.append(row)
            anomaly_rows.append({"episode": i, "post": i >= drift_at, "y": float(ye), "H": row["H"], "anomaly": row["anomaly"], "adjustment": row["adjustment"], "prior": row["prior"], "residual": row["e"]})
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
        "online_anomaly": anomaly_summary(anomaly_rows),
        "roi": yongzeng_ledger(anomaly_rows, pred_rows, children, hops, drift_at),
    }


def yongzeng_ledger(anomaly_rows, pred_rows, children, hops, drift_at: int) -> dict:
    """用增 stacks the prompt quantities on top of terminal Y.

    A channel is 用增 when its post-drift mean is above its pre-drift mean.
    Terminal Y uses that same rule. Child-node count is a cost, not an increment.
    The prompt width and the window length are fixed costs.
    """
    def channel(name, pre, post, role):
        delta = None if pre is None or post is None else float(post - pre)
        if role == "cost":
            kind = "成本"
        elif delta is None:
            kind = None
        elif delta > 1e-9:
            kind = "用增"
        elif delta < -1e-9:
            kind = "挽损" if role == "outcome" else "低于参考"
        else:
            kind = "持平"
        return {"name": name, "role": role, "pre": pre, "post": post, "delta": delta, "kind": kind}

    pre = [r for r in anomaly_rows if r["episode"] < drift_at]
    post = [r for r in anomaly_rows if r["episode"] >= drift_at]

    def avg(rows, key, absolute=False):
        if not rows:
            return None
        vals = np.array([r[key] for r in rows], float)
        if absolute:
            vals = np.abs(vals)
        return float(np.mean(vals))

    pre_hops = [h for h, r in zip(hops, pred_rows) if not r["post"]]
    post_hops = [h for h, r in zip(hops, pred_rows) if r["post"]]
    pre_err = [abs(r["y"] - r["y_hat_last"]) for r in pred_rows if not r["post"]]
    post_err = [abs(r["y"] - r["y_hat_last"]) for r in pred_rows if r["post"]]
    pre_kids = children[:drift_at]
    post_kids = children[drift_at:]

    channels = [
        channel("terminal_Y", avg(pre, "y"), avg(post, "y"), "outcome"),
        channel("abs_residual", avg(pre, "residual", True), avg(post, "residual", True), "monitor"),
        channel("online_anomaly_|H|", avg(pre, "anomaly"), avg(post, "anomaly"), "monitor"),
        channel("prior_e", avg(pre, "prior"), avg(post, "prior"), "monitor"),
        channel(
            "conservative_share",
            float(np.mean([r["adjustment"] == "conservative" for r in pre])) if pre else None,
            float(np.mean([r["adjustment"] == "conservative" for r in post])) if post else None,
            "monitor",
        ),
        channel(
            "hop_rate",
            float(np.mean(pre_hops)) if pre_hops else None,
            float(np.mean(post_hops)) if post_hops else None,
            "monitor",
        ),
        channel(
            "abs_error_last_step",
            float(np.mean(pre_err)) if pre_err else None,
            float(np.mean(post_err)) if post_err else None,
            "monitor",
        ),
        channel(
            "children_scored",
            float(np.mean(pre_kids)) if pre_kids else None,
            float(np.mean(post_kids)) if post_kids else None,
            "cost",
        ),
        channel("prompt_window_rows", float(PROMPT_WINDOW), float(PROMPT_WINDOW), "cost"),
        channel("prompt_columns", float(len(NAMES)), float(len(NAMES)), "cost"),
    ]
    raised = [c["name"] for c in channels if c["kind"] == "用增"]
    return {
        "rule": "post mean above pre mean is 用增. Terminal Y below its pre-drift mean is 挽损. Children, window length, and column count are costs.",
        "channels": channels,
        "yongzeng_channels": raised,
        "delta_children": channels[7]["delta"],
    }


def anomaly_summary(rows: list) -> dict:
    """|H| is the online anomaly score. The sign is the detector.

    Y in {0, 1} and prior e in (0, 1) give H = 1/e when Y=1 and H = -1/(1-e) when Y=0.
    If e < 1/2, a success has a larger |H| than a failure with the same prior.
    Magnitude alone is not a failure detector.
    """
    def pack(subset):
        if not subset:
            return None
        y = np.array([r["y"] for r in subset], float)
        h = np.array([r["H"] for r in subset], float)
        a = np.array([r["anomaly"] for r in subset], float)
        e = np.array([r["prior"] for r in subset], float)
        sign_ok = float(np.mean(np.sign(h) == np.where(y >= 0.5, 1.0, -1.0)))
        return {
            "n_steps": int(len(subset)),
            "mean_prior": float(np.mean(e)),
            "min_prior": float(np.min(e)),
            "mean_abs_H_when_Y1": float(np.mean(a[y >= 0.5])) if np.any(y >= 0.5) else None,
            "mean_abs_H_when_Y0": float(np.mean(a[y < 0.5])) if np.any(y < 0.5) else None,
            "sign_matches_outcome": sign_ok,
            "frac_aggressive_given_Y1": float(np.mean([r["adjustment"] == "aggressive" for r in subset if r["y"] >= 0.5])) if np.any(y >= 0.5) else None,
            "frac_conservative_given_Y0": float(np.mean([r["adjustment"] == "conservative" for r in subset if r["y"] < 0.5])) if np.any(y < 0.5) else None,
        }

    post = [r for r in rows if r["post"]]
    return {
        "score": "|H|",
        "formula": "H=(Y-e)/(e(1-e)), e=stable score of the committed step",
        "other_scores_not_this_one": [
            "residual |Y-RF(X)|",
            "successive MSE hop ratio",
            "mid-depth |judge-stable| versus frozen gaps",
            "path agreement with the first visit",
        ],
        "all_steps": pack(rows),
        "post_drift": pack(post),
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
        print(" anomaly", json.dumps(block["online_anomaly"]["post_drift"], ensure_ascii=False))
        print(" yongzeng", [(c["name"], None if c["delta"] is None else round(c["delta"], 4), c["kind"]) for c in block["roi"]["channels"]])
        print(" rf", block["feature_selection"]["rf_columns"], "residual", block["feature_selection"]["residual_columns"])
        print("--- prompt ---")
        print(block["prompt_example"])
        print()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
