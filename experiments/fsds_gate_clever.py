"""Gate the clever-covariate switch with FSDS.

H = 1/e still marks a step whose judge pick is unlikely under the stable softmax.
The switch fires only when that mark is set and the concept-drift share on
judge + spurious is at least half. Otherwise the step stays with the judge.

Selection efficiency is the AUC of the SSE drops against the label
"this column is spurious", the column the environment actually mixes in.
ROI uses solved-episode counts over the post-drift window. The denominator
is the number of steps that actually switched, not a token price.
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
import tot_agent_probe as probe
from clever_covariate_search import DRIFT_AT, N, REF, SEED, clever_h, softmax
from cot_fsds_next import concept_weights

OUT = Path(__file__).resolve().parent / "results_fsds_roi.json"
MIN_POST = 5


def _auc_spurious(concept: list[dict]) -> float | None:
    if not concept:
        return None
    scores = np.array([max(r["sse_drop"], 0.0) for r in concept], float)
    labels = np.array([1.0 if r["feature"] == "spurious" else 0.0 for r in concept])
    pos = scores[labels == 1.0]
    neg = scores[labels == 0.0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    # P(score_pos > score_neg) + 0.5 P(tie)
    wins = 0.0
    for p in pos:
        wins += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(wins / (len(pos) * len(neg)))


def _allow(feats: list[dict]) -> tuple[bool, float | None, str | None]:
    if len(feats) < REF + MIN_POST:
        return False, None, None
    weights, concept, _cov = concept_weights(feats)
    top = concept[0]["feature"] if concept else None
    return weights["stable"] >= 0.5, _auc_spurious(concept), top


def _choose(scored, policy: str, href: list[float], allow: bool):
    judges = [s[0] for s in scored]
    stables = [s[1] for s in scored]
    e = softmax(stables)
    j = int(np.argmax(judges))
    h = clever_h(float(e[j]))
    if policy == "judge" or len(href) < REF:
        return j, h, False
    thr = float(np.quantile(np.asarray(href[:REF], float), 0.95))
    hot = h > thr + 1e-9
    if policy == "clever" and hot:
        return int(np.argmax(stables)), h, True
    if policy == "gated" and hot and allow:
        return int(np.argmax(stables)), h, True
    return j, h, False


def run_hotpot(policy: str):
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    bank = probe.hop_questions()
    href, feats = [], []
    succ, switches, aucs, tops = [], [], [], []
    for t in range(N):
        allow, auc, top = _allow(feats)
        question, gold = bank[t % len(bank)]
        need = 2
        picked = frozenset()
        last = None
        switched = 0
        for _ in range(2):
            cands = []
            for pid in range(len(probe.PASSAGES)):
                if pid in picked:
                    continue
                nxt = frozenset(set(picked) | {pid})
                progress = len(set(nxt) & gold) / need
                stable = float(np.mean([probe.hop_overlap(question, i) for i in nxt]))
                spurious = float(np.mean([len(probe.PASSAGES[i].split()) / 12 for i in nxt]))
                calibrated = 0.75 * progress + 0.25 * stable
                judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
                cands.append((judge, stable, spurious, progress, nxt))
            idx, h, did = _choose(cands, policy, href, allow)
            if len(href) < REF:
                href.append(h)
            if did:
                switched += 1
            judge, stable, spurious, progress, nxt = cands[idx]
            picked = nxt
            last = (judge, stable, spurious, progress)
        y = 1.0 if len(picked & gold) >= need else 0.0
        succ.append(y)
        switches.append(switched)
        feats.append({
            "judge": last[0], "stable": last[1], "spurious": last[2], "progress": last[3], "y": y,
        })
        if auc is not None:
            aucs.append(auc)
            tops.append(top)
    return _pack(succ, switches, aucs, tops)


def run_game24(policy: str):
    probe.CURRENT_MIX = 0.28
    probe.MIX_JITTER = 0.0
    puzzles = probe.make_game24(random.Random(SEED), N)
    href, feats = [], []
    succ, switches, aucs, tops = [], [], [], []
    for t, puzzle in enumerate(puzzles):
        allow, auc, top = _allow(feats)
        state = puzzle
        last = None
        switched = 0
        y = 0.0
        for _ in range(3):
            kids = probe._g24_children(state)
            if not kids:
                break
            cands = []
            for nxt, done in kids:
                solv = 1.0 if probe._g24_solvable(nxt) else 0.0
                close = max(math.exp(-abs(x - 24) / 8) for x in nxt)
                stable = close if len(nxt) == 1 else 0.45 * close + 0.15
                spurious = sum(nxt) / (13 * len(nxt))
                calibrated = 0.85 * solv + 0.15 * stable
                judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
                cands.append((judge, stable, spurious, solv, nxt, done))
            idx, h, did = _choose(cands, policy, href, allow)
            if len(href) < REF:
                href.append(h)
            if did:
                switched += 1
            judge, stable, spurious, solv, nxt, done = cands[idx]
            state = nxt
            last = (judge, stable, spurious, solv)
            if done:
                y = 1.0
                break
        succ.append(y)
        switches.append(switched)
        if last is None:
            last = (0.0, 0.0, 0.0, 0.0)
        feats.append({
            "judge": last[0], "stable": last[1], "spurious": last[2], "progress": last[3], "y": y,
        })
        if auc is not None:
            aucs.append(auc)
            tops.append(top)
    return _pack(succ, switches, aucs, tops)


def _pack(succ, switches, aucs, tops):
    post = succ[DRIFT_AT:]
    post_sw = switches[DRIFT_AT:]
    counts = {}
    for name in tops:
        counts[name] = counts.get(name, 0) + 1
    return {
        "pre_success": float(np.mean(succ[:DRIFT_AT])),
        "post_success": float(np.mean(post)),
        "post_switches": int(np.sum(post_sw)),
        "post_episodes": int(len(post)),
        "selection_auc": None if not aucs else float(np.mean(aucs)),
        "concept_top_counts": counts,
    }


def roi(base, alt, v_task: float):
    extra = (alt["post_success"] - base["post_success"]) * alt["post_episodes"]
    switches = alt["post_switches"]
    revenue = extra * v_task
    return {
        "extra_solved": extra,
        "post_switches": switches,
        "revenue": revenue,
        "roi_per_switch": None if switches == 0 else extra / switches,
        "selection_auc": alt["selection_auc"],
    }


def main():
    report = {"tasks": {}}
    runners = (("hotpot", run_hotpot, 0.001), ("game24", run_game24, 0.10))
    for name, fn, price in runners:
        report["tasks"][name] = {}
        for policy in ("judge", "clever", "gated"):
            print(name, policy, flush=True)
            block = fn(policy)
            report["tasks"][name][policy] = block
            print(
                f"  pre={block['pre_success']:.2f} post={block['post_success']:.2f}"
                f" switches={block['post_switches']} auc={block['selection_auc']}"
            )
        report["tasks"][name]["roi_vs_judge"] = roi(report["tasks"][name]["judge"], report["tasks"][name]["gated"], price)
        report["tasks"][name]["roi_clever_vs_judge"] = roi(
            report["tasks"][name]["judge"], report["tasks"][name]["clever"], price
        )
        print(" ", name, "gated", report["tasks"][name]["roi_vs_judge"])
    OUT.write_text(json.dumps(report, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
