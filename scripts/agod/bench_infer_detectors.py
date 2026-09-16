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


def oob_stats(det: np.ndarray, *, n_per: int, cut: int, method: str) -> dict:
    """Stats matching onlinePermOOB_algo.py: SUM + first1/2/3 consecutive."""
    det = np.asarray(det, dtype=bool)
    quiet_end = cut * n_per
    first1 = first_k_consecutive_rej(det, 1)
    first2 = first_k_consecutive_rej(det, 2)
    first3 = first_k_consecutive_rej(det, 3)

    def _delay_batch(first_t: int | None) -> int | None:
        if first_t is None:
            return None
        # delay in batches relative to cut (negative ⇒ quiet false alarm)
        return int(first_t // n_per - cut)

    return {
        "method": method,
        "SUM": int(np.sum(det)),
        "first1": first1,
        "first2": first2,
        "first3": first3,
        "delay1": _delay_batch(first1),
        "delay2": _delay_batch(first2),
        "delay3": _delay_batch(first3),
        "quiet_false_fires": int(np.sum(det[:quiet_end])) if quiet_end > 0 else 0,
        # keep old keys for tests
        "first_fire_batch": None if first1 is None else int(first1 // n_per),
        "detection_delay": _delay_batch(first1),
        "n_fires": int(np.sum(det)),
        "first_t": first1,
    }


def detect_online_rfperm(records: list[dict], *, n_per: int, cut: int, gate: float) -> dict:
    """OnlineRFPerm batch fires → expand to per-t mask, then same first1/2/3 stats."""
    X, y = featurize_records(records)
    n_batches = len(records) // n_per
    batch_fired = np.zeros(n_batches, dtype=bool)
    e_prev = None
    X_prev = y_prev = None
    for b in range(n_batches):
        sl = slice(b * n_per, (b + 1) * n_per)
        Xb, yb = X[sl], y[sl]
        fired = False
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
            e_prev = float(e_now)
        else:
            probe = fit_online_probe(Xb, yb, task="acc", seed=b)
            e_prev = float(probe_err(probe, Xb, yb, task="acc"))
        batch_fired[b] = fired
        X_prev, y_prev = Xb, yb

    # expand: if batch fires, mark the whole batch True (for consecutive-k on t-grid)
    det = np.repeat(batch_fired, n_per)
    return oob_stats(det, n_per=n_per, cut=cut, method="OnlineRFPerm")


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
        rows.append(oob_stats(det, n_per=n_per, cut=cut, method=method))

    for r in rows:
        r["dataset"] = name
        r["n"] = len(records)
        r["n_per"] = n_per
        r["cut"] = cut
        r["cut_t"] = ref_end
        r["ref_end"] = ref_end
        r["mse_mean_quiet"] = float(np.mean(mse[:ref_end]))
        r["mse_mean_hop"] = float(np.mean(mse[ref_end:]))
    return rows


def _fmt(v) -> str:
    return "---" if v is None else str(int(v))


def write_latex(results: dict, out: Path) -> str:
    """Two tables: (i) first1/2/3 indices  (ii) delay1/2/3 in batches after cut."""
    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=1in]{geometry}",
        r"\usepackage{booktabs,amsmath}",
        r"\title{OnlineRFPerm vs OOB baselines\\first$_1$/first$_2$/first$_3$}",
        r"\author{CFPerm benchmark (onlinePermOOB\_algo stats)}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        "",
        r"\paragraph{Setup.}",
        (
            r"Same HaluEval/SQuAD CPU streams. Baselines from "
            r"\texttt{onlinePermOOB\_algo.py}. Stats: "
            r"\texttt{first}$_k=$ first index of $k$ consecutive rejections; "
            r"\texttt{delay}$_k=$ batch delay vs cut ($t/n_{\mathrm{per}}-\mathrm{cut}$; "
            r"negative $=$ quiet false alarm). Cut $t$ marked in caption."
        ),
        "",
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{first$_1$/first$_2$/first$_3$ (stream index) and SUM.}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule dataset & method & SUM & first$_1$ & first$_2$ & first$_3$ \\",
        r"\midrule",
    ]
    for ds, methods in results.items():
        cut_t = methods[0]["cut_t"]
        lines.append(rf"\multicolumn{{6}}{{l}}{{\emph{{{latex_escape(ds)}}}, cut $t={cut_t}$}} \\")
        for m in methods:
            lines.append(
                f"{latex_escape(ds)} & {latex_escape(m['method'])} & {m['SUM']} & "
                f"{_fmt(m['first1'])} & {_fmt(m['first2'])} & {_fmt(m['first3'])} \\\\"
            )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{delay$_1$/delay$_2$/delay$_3$ (batches after cut).}",
        r"\begin{tabular}{llrrr}",
        r"\toprule dataset & method & delay$_1$ & delay$_2$ & delay$_3$ \\",
        r"\midrule",
    ]
    for ds, methods in results.items():
        for m in methods:
            lines.append(
                f"{latex_escape(ds)} & {latex_escape(m['method'])} & "
                f"{_fmt(m['delay1'])} & {_fmt(m['delay2'])} & {_fmt(m['delay3'])} \\\\"
            )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\paragraph{Takeaway.}",
        (
            r"Same counting rule as the OOB script: first$_1$/$_2$/$_3$. "
            r"OnlineRFPerm expands a fired batch to $n_{\mathrm{per}}$ consecutive "
            r"trues so the consecutive-$k$ stats stay on the same $t$-grid."
        ),
        "",
        r"\end{document}",
        "",
    ]
    tex = "\n".join(lines)
    (out / "bench_detectors.tex").write_text(tex)
    (DOCS / "BENCH_INFER_DETECTORS.tex").write_text(tex)
    return tex


def latex_escape(s: str) -> str:
    return str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


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
        "# Detector benchmark (onlinePermOOB first1/2/3)\n",
        "| dataset | method | SUM | first1 | first2 | first3 | delay1 | delay2 | delay3 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for ds, methods in results.items():
        for m in methods:
            md.append(
                f"| {ds} | {m['method']} | {m['SUM']} | {m['first1']} | {m['first2']} | "
                f"{m['first3']} | {m['delay1']} | {m['delay2']} | {m['delay3']} |"
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
