#!/usr/bin/env python3
"""Streaming (continuous) OnlineRFPerm tests on LLM infer streams.

Default batch clock uses n_per=20. This script makes detection *continuous*:

  1) ORF-freeze  — fit shallow RF on quiet reference; stream per-obs 0/1 error;
                   fire when rolling mean error / quiet mean ≥ γ
  2) ORF-slide   — slide a window of size ``win`` by 1; consecutive OOS
                   hop_fires (the OnlineRFPerm ratio gate, continuous)

Read (product):
  RF is the component model. If a quality-law hop is present, continuous ORF
  should fire. If it does not, the serving chain looks *smooth* under this
  probe (no usable P(Y|X) hop) — assuming labels/features are honest.

Usage::

  PYTHONPATH=. python3 scripts/agod/online_rfperm_streaming_test.py
  PYTHONPATH=. python3 scripts/agod/online_rfperm_streaming_test.py --win 20 --gate 1.25
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]
ROOT = Path(__file__).resolve().parents[2]
STREAM_ROOT = ROOT / "results" / "agod" / "online_rfperm_multi_datasets"
OUT = ROOT / "results" / "agod" / "online_rfperm_streaming_test"
DOCS = ROOT / "docs" / "biz"

sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402
from agod.online_rfperm import (  # noqa: E402
    error_floor,
    fit_online_probe,
    hop_fires,
    probe_err,
    shift_ratio,
)

DATASETS = ("halueval", "squad", "hotpotqa", "truthfulqa")


def first_k_consecutive(det, k: int) -> int | None:
    det = np.asarray(det, dtype=bool).ravel()
    if k <= 0 or len(det) < k:
        return None
    for i in range(len(det) - k + 1):
        if bool(det[i : i + k].all()):
            return int(i)
    return None


def load_records(name: str) -> list[dict]:
    path = STREAM_ROOT / name / "stream.jsonl"
    if not path.exists():
        raise SystemExit(f"missing {path}; run online_rfperm_multi_datasets.py first")
    rows = []
    with path.open() as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def featurize(records: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    X = live.featurize(
        [r["question"] for r in records],
        [r["answer"] for r in records],
        [r["rag_hit"] for r in records],
    )
    y = np.asarray([r["y_bad"] for r in records], dtype=int)
    return X, y


def rolling_mean(x: np.ndarray, win: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.full(len(x), np.nan)
    if win <= 0 or len(x) < win:
        return out
    c = np.cumsum(np.insert(x, 0, 0.0))
    out[win - 1 :] = (c[win:] - c[:-win]) / win
    return out


def orf_freeze_stream(
    X: np.ndarray,
    y: np.ndarray,
    *,
    ref_end: int,
    win: int,
    gate: float,
) -> dict:
    """Freeze RF on quiet ref; continuous rolling-error ratio on the trail."""
    probe = fit_online_probe(X[:ref_end], y[:ref_end], task="acc", seed=0)
    err = (probe.predict(X).astype(int) != y).astype(float)
    quiet_mean = float(np.mean(err[:ref_end]))
    fl = error_floor("acc", ref_end)
    denom = max(quiet_mean, fl)
    roll = rolling_mean(err, win)
    fires = np.zeros(len(err), dtype=bool)
    ratios = np.full(len(err), np.nan)
    for t in range(ref_end, len(err)):
        if not np.isfinite(roll[t]):
            continue
        ratio = float(roll[t]) / denom
        ratios[t] = ratio
        # continuous: fire when rolling hop-window error exceeds quiet by γ
        fires[t] = bool(np.isfinite(ratio) and ratio >= gate and roll[t] >= fl)
    trail = fires[ref_end:]
    return {
        "mode": "orf_freeze",
        "quiet_err_mean": quiet_mean,
        "hop_err_mean": float(np.mean(err[ref_end:])),
        "fires": fires,
        "trail_fires": trail,
        "ratios": ratios,
        "first1": first_k_consecutive(trail, 1),
        "first2": first_k_consecutive(trail, 2),
        "first3": first_k_consecutive(trail, 3),
        "n_trail_fires": int(trail.sum()),
        "detection_delay_obs": (
            None if first_k_consecutive(trail, 1) is None else int(first_k_consecutive(trail, 1))
        ),
    }


def orf_slide_stream(
    X: np.ndarray,
    y: np.ndarray,
    *,
    ref_end: int,
    win: int,
    gate: float,
) -> dict:
    """Slide window by 1; consecutive OOS OnlineRFPerm gate (continuous)."""
    n = len(y)
    fires = np.zeros(n, dtype=bool)
    ratios = np.full(n, np.nan)
    e_prev = None
    fl = error_floor("acc", win)
    # need two consecutive windows ending at t: [t-2w, t-w) and [t-w, t)
    t0 = max(ref_end, 2 * win)
    for t in range(t0, n + 1):
        prev = slice(t - 2 * win, t - win)
        now = slice(t - win, t)
        Xp, yp = X[prev], y[prev]
        Xn, yn = X[now], y[now]
        # skip degenerate constant-label windows that make RF vacuous
        if len(np.unique(yp)) < 2 and float(np.mean(yp)) in (0.0, 1.0):
            # still score e_now vs floor when entering hop
            probe = fit_online_probe(Xp, yp, task="acc", seed=t)
            e_now = probe_err(probe, Xn, yn, task="acc")
        else:
            probe = fit_online_probe(Xp, yp, task="acc", seed=t)
            e_now = probe_err(probe, Xn, yn, task="acc")
        hit = False
        if e_prev is not None:
            hit = bool(hop_fires(e_now, e_prev, gate=gate, e_floor=fl))
            if (
                not hit
                and float(e_prev) < fl
                and float(e_now) >= fl
                and float(e_now) / fl >= gate
            ):
                # quiet→hop escape when previous window was near-perfect
                hit = True
            ratios[t - 1] = float(shift_ratio(e_now, max(float(e_prev), fl), e_floor=fl))
        e_prev = float(e_now)
        # attribute fire to the last index of the current window
        fires[t - 1] = hit
    trail = fires[ref_end:]
    return {
        "mode": "orf_slide",
        "fires": fires,
        "trail_fires": trail,
        "ratios": ratios,
        "first1": first_k_consecutive(trail, 1),
        "first2": first_k_consecutive(trail, 2),
        "first3": first_k_consecutive(trail, 3),
        "n_trail_fires": int(trail.sum()),
        "detection_delay_obs": (
            None if first_k_consecutive(trail, 1) is None else int(first_k_consecutive(trail, 1))
        ),
    }


def smooth_control(
    X: np.ndarray,
    y: np.ndarray,
    *,
    ref_end: int,
    win: int,
    gate: float,
) -> dict:
    """All-quiet control: replace hop y with quiet resampling → expect no fire."""
    y_s = y.copy()
    quiet_y = y[:ref_end]
    rng = np.random.default_rng(0)
    y_s[ref_end:] = quiet_y[rng.integers(0, ref_end, size=len(y) - ref_end)]
    # also clone quiet-ish features on hop for a stronger smooth control
    X_s = X.copy()
    idx = rng.integers(0, ref_end, size=len(y) - ref_end)
    X_s[ref_end:] = X[idx]
    fr = orf_freeze_stream(X_s, y_s, ref_end=ref_end, win=win, gate=gate)
    sl = orf_slide_stream(X_s, y_s, ref_end=ref_end, win=win, gate=gate)
    return {
        "freeze_first1": fr["first1"],
        "slide_first1": sl["first1"],
        "freeze_n_fires": fr["n_trail_fires"],
        "slide_n_fires": sl["n_trail_fires"],
        "read": "smooth_control_should_stay_quiet",
    }


def run_one(name: str, *, win: int, gate: float) -> dict:
    records = load_records(name)
    X, y = featurize(records)
    cut_batch = next(r["batch"] for r in records if r.get("hopped"))
    n_per = sum(1 for r in records if r["batch"] == 0)
    ref_end = cut_batch * n_per
    hopped = np.asarray([bool(r["hopped"]) for r in records])

    freeze = orf_freeze_stream(X, y, ref_end=ref_end, win=win, gate=gate)
    slide = orf_slide_stream(X, y, ref_end=ref_end, win=win, gate=gate)
    control = smooth_control(X, y, ref_end=ref_end, win=win, gate=gate)

    # serialize without huge arrays
    def pack(d: dict) -> dict:
        return {
            "mode": d["mode"],
            "first1": d["first1"],
            "first2": d["first2"],
            "first3": d["first3"],
            "n_trail_fires": d["n_trail_fires"],
            "detection_delay_obs": d["detection_delay_obs"],
            **(
                {
                    "quiet_err_mean": d.get("quiet_err_mean"),
                    "hop_err_mean": d.get("hop_err_mean"),
                }
                if "quiet_err_mean" in d
                else {}
            ),
        }

    hop_y = float(np.mean(y[hopped])) if hopped.any() else float("nan")
    quiet_y = float(np.mean(y[~hopped])) if (~hopped).any() else float("nan")
    return {
        "dataset": name,
        "n": len(records),
        "n_per_source": n_per,
        "cut_batch": cut_batch,
        "ref_end": ref_end,
        "win": win,
        "gate": gate,
        "y_quiet": quiet_y,
        "y_hop": hop_y,
        "orf_freeze": pack(freeze),
        "orf_slide": pack(slide),
        "smooth_control": control,
        "read": _read(freeze, slide, control, quiet_y, hop_y),
    }


def _read(freeze, slide, control, y_q, y_h) -> str:
    caught = freeze["first1"] is not None or slide["first1"] is not None
    smooth_ok = control["freeze_first1"] is None and control["slide_first1"] is None
    if y_h - y_q < 0.2:
        return "weak_label_contrast_check_y"
    if caught and smooth_ok:
        return "hop_caught_smooth_control_quiet — RF component sees the law change"
    if caught and not smooth_ok:
        return "hop_caught_but_control_noisy — gate/window may be hot"
    if not caught and smooth_ok:
        return "no_fire_chain_looks_smooth_under_probe — or window/gate too strict"
    return "no_fire_and_control_noisy — revisit features/labels"


def write_report(summaries: list[dict], args: argparse.Namespace, out: Path) -> str:
    rows = []
    for s in summaries:
        fr, sl = s["orf_freeze"], s["orf_slide"]
        rows.append(
            f"| {s['dataset']} | {s['y_quiet']:.2f}→{s['y_hop']:.2f} | "
            f"{fr['detection_delay_obs']} | {sl['detection_delay_obs']} | "
            f"{s['smooth_control']['freeze_first1']}/{s['smooth_control']['slide_first1']} | "
            f"{s['read']} |"
        )
    body = "\n".join(rows)
    md = f"""# Streaming OnlineRFPerm test (continuous)

> RF shallow probe = component model on the tidy infer stream.  
> Continuous clock: freeze + rolling ratio, and sliding consecutive OOS.  
> **If a hop is present and RF does not fire → under this probe the chain looks smooth.**

## Setup

| item | value |
|------|-------|
| source streams | `results/agod/online_rfperm_multi_datasets/*/stream.jsonl` |
| faith / Y | answer-precision → binary $Y$ |
| win / gate | {args.win} / {args.gate} |
| ref | quiet window before cut (`n_ref=100` on default multi) |

## Modes

1. **orf_freeze** — fit RF on quiet; stream per-obs error; fire when rolling mean / quiet ≥ γ  
2. **orf_slide** — slide window by 1; OnlineRFPerm consecutive OOS gate  
3. **smooth_control** — rewrite hop rows as quiet resamples; expect **no** fire

Delay = first trail observation index with fire (`0` = first obs after cut).

## Results

| dataset | y quiet→hop | freeze delay | slide delay | smooth first1 (fr/sl) | read |
|---------|-------------|-------------:|------------:|-----------------------|------|
{body}

## How to read

- Hop injected + fire soon → RF component **caught** the quality-law change.  
- Hop injected + **no** fire → inference chain looks **丝滑 / smooth** under this probe (or labels/features hide the hop).  
- Smooth control should stay quiet; if it fires, the gate is too hot.

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_streaming_test.py --win {args.win} --gate {args.gate}
```
"""
    (out / "REPORT.md").write_text(md)
    (DOCS / "ONLINERFPERM_STREAMING_TEST.md").write_text(md)

    # latex one-pager
    tex_rows = []
    for s in summaries:
        fr, sl = s["orf_freeze"], s["orf_slide"]
        fd = "---" if fr["detection_delay_obs"] is None else str(fr["detection_delay_obs"])
        sd = "---" if sl["detection_delay_obs"] is None else str(sl["detection_delay_obs"])
        tex_rows.append(
            f"{s['dataset']} & {s['y_quiet']:.2f}$\\to${s['y_hop']:.2f} & {fd} & {sd} \\\\"
        )
    tex = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs,amsmath}",
            r"\title{Streaming OnlineRFPerm on LLM infer streams}",
            r"\author{CFPerm}",
            r"\date{\today}",
            r"\begin{document}",
            r"\maketitle",
            rf"Continuous clock: win$={args.win}$, $\gamma={args.gate}$. "
            r"Freeze = RF on quiet + rolling error ratio; "
            r"slide = consecutive OOS OnlineRFPerm with step 1.",
            r"\begin{table}[h]\centering",
            r"\caption{Streaming detection delay (trail observation index).}",
            r"\begin{tabular}{lrrr}",
            r"\toprule dataset & $y$ quiet$\to$hop & freeze delay & slide delay \\",
            r"\midrule",
            *tex_rows,
            r"\bottomrule\end{tabular}\end{table}",
            r"\paragraph{Read.}",
            r"RF component model should catch an injected quality hop. "
            r"No fire $\Rightarrow$ chain looks smooth under this probe.",
            r"\end{document}",
            "",
        ]
    )
    (out / "streaming_test.tex").write_text(tex)
    (DOCS / "ONLINERFPERM_STREAMING_TEST.tex").write_text(tex)
    return md


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--win", type=int, default=20, help="rolling / slide window size")
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--datasets", nargs="+", default=list(DATASETS))
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    summaries = [run_one(name, win=args.win, gate=args.gate) for name in args.datasets]
    (args.out / "summary.json").write_text(
        json.dumps(
            {
                "stance": "continuous streaming OnlineRFPerm; RF component; smooth=no-fire",
                "win": args.win,
                "gate": args.gate,
                "datasets": summaries,
            },
            indent=2,
        )
    )
    write_report(summaries, args, args.out)
    for s in summaries:
        print(
            f"[{s['dataset']}] freeze_delay={s['orf_freeze']['detection_delay_obs']} "
            f"slide_delay={s['orf_slide']['detection_delay_obs']} "
            f"smooth={s['smooth_control']['freeze_first1']}/{s['smooth_control']['slide_first1']} "
            f"| {s['read']}",
            flush=True,
        )
    print(f"[streaming] wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
