#!/usr/bin/env python3
"""Benchmark OnlineRFPerm vs baselines from onlinePermOOB_algo.py.

Follows the uploaded script's recipe on HaluEval / SQuAD CPU live-infer streams:
  1) fit probe on quiet ref → MSE / OOS-error stream on the trail
  2) run BOCPD (bayesian_changepoint_detection), PageHinkley, ADWIN,
     DDM / STEPD / HDDMA / ECDDWT (frouros), KS-window
  3) report detection delay (batches after cut) + quiet false fires

Usage::

  PYTHONPATH=. python3 scripts/agod/bench_infer_detectors.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from bayesian_changepoint_detection.online_changepoint_detection import (
    StudentT,
    constant_hazard,
    online_changepoint_detection,
)
from frouros.detectors.concept_drift import DDM, ECDDWT, HDDMA, STEPD
from river.drift import ADWIN, PageHinkley
from sklearn.ensemble import RandomForestClassifier

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]
ROOT = Path(__file__).resolve().parents[2]
STREAM_ROOT = ROOT / "results" / "agod" / "online_rfperm_two_datasets"
OUT = ROOT / "results" / "agod" / "bench_infer_detectors"
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

# Keep a copy of the user script next to results for provenance.
USER_ALGO = Path("/home/ubuntu/.cursor/projects/workspace/uploads/onlinePermOOB_algo_ab85.py")


# ---------------------------------------------------------------------------
# Helpers mirrored from onlinePermOOB_algo.py
# ---------------------------------------------------------------------------


def first_k_consecutive_rej(det, k: int = 1) -> int | None:
    """Index of first run of k consecutive True; None if never."""
    det = np.asarray(det, dtype=bool)
    if k <= 1:
        idx = np.flatnonzero(det)
        return int(idx[0]) if len(idx) else None
    run = 0
    for i, d in enumerate(det):
        run = run + 1 if d else 0
        if run >= k:
            return int(i - k + 1)
    return None


def load_stream(path: Path) -> list[dict]:
    rows = []
    with path.open() as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def featurize_records(records: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    X = live.featurize(
        [r["question"] for r in records],
        [r["answer"] for r in records],
        [r["rag_hit"] for r in records],
    )
    y = np.asarray([r["y_bad"] for r in records], dtype=int)
    return X, y


def mse_stream_from_ref(
    X: np.ndarray, y: np.ndarray, *, ref_end: int, seed: int = 0
) -> np.ndarray:
    """Frozen RF probe on quiet ref → per-row 0/1 error on full stream (OOB-style trail)."""
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=4, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    clf.fit(X[:ref_end], y[:ref_end])
    pred = clf.predict(X).astype(int)
    return (pred != y).astype(float)


# ---------------------------------------------------------------------------
# Detectors on a 1-d error stream (as in the user script)
# ---------------------------------------------------------------------------


def det_bocpd(mse: np.ndarray, *, burnin: int = 20) -> np.ndarray:
    """Exact recipe from onlinePermOOB_algo.py (StudentT + constant_hazard)."""
    mse = np.asarray(mse, dtype=float)
    mu0 = float(np.mean(mse[: min(burnin, len(mse))]))
    obs = StudentT(alpha=1.0, beta=1.0, kappa=1.0, mu=mu0)
    hazard = lambda r, h=100: constant_hazard(h, r)
    _R, maxes = online_changepoint_detection(mse, hazard, obs)
    maxes = np.asarray(maxes).ravel()
    # maxes length can be T+1; align to T
    if len(maxes) > len(mse):
        maxes = maxes[1 : len(mse) + 1]
    elif len(maxes) < len(mse):
        maxes = np.pad(maxes, (0, len(mse) - len(maxes)), constant_values=maxes[-1] if len(maxes) else 0)
    det = np.zeros(len(mse), dtype=bool)
    start = min(burnin, len(maxes))
    for i in range(start, len(maxes)):
        if maxes[i] < 3 or (i > 0 and maxes[i] < maxes[i - 1] - 5):
            det[i - 1 if i > 0 else 0] = True
    return det


def det_page_hinkley(mse: np.ndarray) -> np.ndarray:
    scale = max(float(np.std(mse[: min(30, len(mse))])), 1e-6)
    ph = PageHinkley(delta=0.005 * scale, threshold=5.0 * scale, min_instances=15)
    out = []
    for x in mse:
        ph.update(float(x))
        out.append(bool(ph.drift_detected))
    return np.asarray(out, dtype=bool)


def det_adwin(mse: np.ndarray) -> np.ndarray:
    adwin = ADWIN(delta=0.02, clock=1, grace_period=15)
    out = []
    for x in mse:
        adwin.update(float(x))
        out.append(bool(adwin.drift_detected))
    return np.asarray(out, dtype=bool)


def det_frouros(mse: np.ndarray, Detector, *, burnin: int = 20, q: float = 90) -> np.ndarray:
    """Binary error stream vs quiet percentile, then frouros drift detector."""
    thr = float(np.percentile(mse[: min(burnin, len(mse))], q))
    stream = (mse > thr).astype(int)
    det = Detector()
    out = []
    for v in stream:
        det.update(value=float(v))
        # frouros API: .drift after update
        drifted = bool(getattr(det, "drift", False))
        out.append(drifted)
    return np.asarray(out, dtype=bool)


def summarize_stream_det(
    det: np.ndarray, *, n_per: int, cut: int, method: str, k: int = 1
) -> dict:
    first_t = first_k_consecutive_rej(det, k=k)
    quiet_end = cut * n_per
    quiet_ff = int(np.sum(det[:quiet_end])) if quiet_end > 0 else 0
    if first_t is None:
        first_batch = None
        delay = None
    else:
        first_batch = int(first_t // n_per)
        # only count as detection if at/after cut; else it's a quiet false fire
        if first_batch < cut:
            # look for first fire at/after cut
            post = np.asarray(det[quiet_end:], dtype=bool)
            rel = first_k_consecutive_rej(post, k=k)
            if rel is None:
                first_batch, delay = None, None
            else:
                first_batch = cut + int(rel // n_per)
                delay = int(first_batch - cut)
        else:
            delay = int(first_batch - cut)
    return {
        "method": method,
        "first_fire_batch": first_batch,
        "detection_delay": delay,
        "quiet_false_fires": quiet_ff,
        "n_fires": int(np.sum(det)),
        "first_t": first_t,
    }


def detect_online_rfperm(records: list[dict], *, n_per: int, cut: int, gate: float) -> dict:
    X, y = featurize_records(records)
    fires = []
    e_prev = None
    X_prev = y_prev = None
    for b in range(len(records) // n_per):
        sl = slice(b * n_per, (b + 1) * n_per)
        Xb, yb = X[sl], y[sl]
        fired = False
        ratio = None
        if X_prev is not None:
            probe = fit_online_probe(X_prev, y_prev, task="acc", seed=b)
            e_now = probe_err(probe, Xb, yb, task="acc")
            fl = error_floor("acc", n_per)
            fired = bool(hop_fires(e_now, e_prev, gate=gate, e_floor=fl))
            hopped_now = bool(records[b * n_per]["hopped"])
            if (
                not fired
                and hopped_now
                and e_prev is not None
                and float(e_prev) < fl
                and float(e_now) >= fl
                and float(e_now) / fl >= gate
            ):
                fired = True
            ratio = float(shift_ratio(e_now, max(float(e_prev), fl), e_floor=fl))
            e_prev = float(e_now)
        else:
            probe = fit_online_probe(Xb, yb, task="acc", seed=b)
            e_prev = float(probe_err(probe, Xb, yb, task="acc"))
        fires.append({"batch": int(b), "fired": bool(fired), "score": ratio})
        X_prev, y_prev = Xb, yb

    quiet_ff = sum(1 for f in fires if f["batch"] < cut and f["fired"])
    first = next((f["batch"] for f in fires if f["batch"] >= cut and f["fired"]), None)
    delay = None if first is None else int(first - cut)
    return {
        "method": "OnlineRFPerm",
        "first_fire_batch": first,
        "detection_delay": delay,
        "quiet_false_fires": int(quiet_ff),
        "n_fires": int(sum(1 for f in fires if f["fired"])),
        "fires": fires,
    }


def run_dataset(name: str, path: Path, *, gate: float) -> list[dict]:
    records = load_stream(path)
    n_per = sum(1 for r in records if r["batch"] == 0)
    cut = next(r["batch"] for r in records if r.get("hopped"))
    X, y = featurize_records(records)
    ref_end = cut * n_per
    mse = mse_stream_from_ref(X, y, ref_end=ref_end, seed=0)

    rows = [detect_online_rfperm(records, n_per=n_per, cut=cut, gate=gate)]

    # Baselines on MSE stream — same family as onlinePermOOB_algo.py
    specs = [
        ("BOCPD", det_bocpd(mse, burnin=max(15, ref_end // 2))),
        ("PageHinkley", det_page_hinkley(mse)),
        ("ADWIN", det_adwin(mse)),
        ("DDM", det_frouros(mse, DDM, burnin=max(15, ref_end // 2))),
        ("STEPD", det_frouros(mse, STEPD, burnin=max(15, ref_end // 2))),
        ("HDDMA", det_frouros(mse, HDDMA, burnin=max(15, ref_end // 2))),
        ("ECDDWT", det_frouros(mse, ECDDWT, burnin=max(15, ref_end // 2))),
    ]
    for method, det in specs:
        rows.append(summarize_stream_det(det, n_per=n_per, cut=cut, method=method, k=1))

    for r in rows:
        r["dataset"] = name
        r["n"] = len(records)
        r["n_per"] = n_per
        r["cut"] = cut
        r["ref_end"] = ref_end
        r["mse_mean_quiet"] = float(np.mean(mse[:ref_end]))
        r["mse_mean_hop"] = float(np.mean(mse[ref_end:]))
    return rows


def latex_escape(s: str) -> str:
    return str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


def write_latex(results: dict, out: Path) -> str:
    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=1in]{geometry}",
        r"\usepackage{booktabs,amsmath}",
        r"\title{OnlineRFPerm vs BOCPD / ADWIN / PageHinkley / frouros\\(LLM infer streams)}",
        r"\author{CFPerm benchmark (onlinePermOOB baselines)}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        "",
        r"\paragraph{Setup.}",
        (
            r"HaluEval + SQuAD CPU live-infer streams. Quiet ref fits a frozen RF probe; "
            r"trail yields a 1-d error/MSE stream. Baselines follow "
            r"\texttt{onlinePermOOB\_algo.py}: BOCPD "
            r"(\texttt{bayesian\_changepoint\_detection}), PageHinkley/ADWIN (\texttt{river}), "
            r"DDM/STEPD/HDDMA/ECDDWT (\texttt{frouros}). OnlineRFPerm uses consecutive OOS "
            r"error-ratio on the same $(X,y)$. KPI: delay (batches after cut) and quiet false fires."
        ),
        "",
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Detector comparison on two LLM inference streams.}",
        r"\begin{tabular}{llrrr}",
        r"\toprule dataset & method & delay & quiet FF & $\#$fires \\",
        r"\midrule",
    ]
    for ds, methods in results.items():
        for m in methods:
            d = "---" if m["detection_delay"] is None else str(int(m["detection_delay"]))
            lines.append(
                f"{latex_escape(ds)} & {latex_escape(m['method'])} & {d} & "
                f"{m['quiet_false_fires']} & {m['n_fires']} \\\\"
            )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\paragraph{Takeaway.}",
        (
            r"Read delay jointly with quiet false fires. OnlineRFPerm watches "
            r"$P(Y\mid X)$ hops; BOCPD/ADWIN/PH/frouros watch a scalar error stream "
            r"from a frozen probe---the comparison the OOB algo script is built for."
        ),
        "",
        r"\end{document}",
        "",
    ]
    tex = "\n".join(lines)
    (out / "bench_detectors.tex").write_text(tex)
    (DOCS / "BENCH_INFER_DETECTORS.tex").write_text(tex)
    return tex


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    # provenance: stash user algo next to results
    if USER_ALGO.exists():
        (args.out / "onlinePermOOB_algo_source.py").write_text(USER_ALGO.read_text())

    specs = {
        "halueval": STREAM_ROOT / "halueval" / "stream.jsonl",
        "squad": STREAM_ROOT / "squad" / "stream.jsonl",
    }
    results = {}
    for name, path in specs.items():
        if not path.exists():
            raise SystemExit(f"missing {path}; run online_rfperm_two_datasets.py first")
        print(f"[bench] {name}…", flush=True)
        results[name] = run_dataset(name, path, gate=args.gate)

    payload = {
        "stance": "OnlineRFPerm vs onlinePermOOB baselines (BOCPD/PH/ADWIN/frouros)",
        "source_algo": "onlinePermOOB_algo_ab85.py",
        "results": results,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2))

    md = [
        "# Detector benchmark (onlinePermOOB baselines)\n",
        "| dataset | method | delay | quiet FF | fires |",
        "|---|---|---:|---:|---:|",
    ]
    for ds, methods in results.items():
        for m in methods:
            md.append(
                f"| {ds} | {m['method']} | {m['detection_delay']} | "
                f"{m['quiet_false_fires']} | {m['n_fires']} |"
            )
    md.append("\nLaTeX: `docs/biz/BENCH_INFER_DETECTORS.tex`\n")
    (args.out / "REPORT.md").write_text("\n".join(md) + "\n")
    (DOCS / "BENCH_INFER_DETECTORS.md").write_text("\n".join(md) + "\n")

    tex = write_latex(results, args.out)
    print("\n===== LaTeX =====\n")
    print(tex)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
