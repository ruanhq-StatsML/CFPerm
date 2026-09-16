#!/usr/bin/env python3
"""Benchmark with manuscript first_k formulation (OnlineRFPerm_PointCloud).

Manuscript rule (PC1/PC3 note)
------------------------------
  • Freeze probe on reference (quiet) window.
  • Score the *trail* only (after cut).
  • ``first_k`` = first trail index where there are k consecutive fires
    starting at that index  (k=1 headline; k=2/3 sustained).
  • Here the primary grid is **batch** on the trail (user ask:
    「第几个 batch 开始头一回两个 fire」).

  first_k_consecutive(det, k):
      for i in range(len(det) - k + 1):
          if det[i:i+k].all():
              return i   # trail batch index; 0 = first batch after cut
      return None

Usage::

  PYTHONPATH=. python3 scripts/agod/bench_infer_detectors.py
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
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
DATA_EXPORT = ROOT / "data" / "hf_cache" / "infer_bench_export"

sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402
from agod.online_rfperm import (  # noqa: E402
    error_floor,
    fit_online_probe,
    hop_fires,
    probe_err,
)

USER_ALGO = Path("/home/ubuntu/.cursor/projects/workspace/uploads/onlinePermOOB_algo_ab85.py")


# ---------------------------------------------------------------------------
# Manuscript first_k (THIS is the logic)
# ---------------------------------------------------------------------------


def first_k_consecutive(det, k: int) -> int | None:
    """First index i on the trail where det[i],…,det[i+k-1] are all True.

    Matches onlinePermOOB / PointCloud manuscript:
    first rejection (k=1) and sustained onsets (k=2, k=3).
    Index is 0-based on the *trail* (after reference / cut).
    """
    det = np.asarray(det, dtype=bool).ravel()
    if k <= 0 or len(det) < k:
        return None
    for i in range(len(det) - k + 1):
        if bool(det[i : i + k].all()):
            return int(i)
    return None


def first_k_stats(det_trail: np.ndarray, *, method: str) -> dict:
    return {
        "method": method,
        "first1": first_k_consecutive(det_trail, 1),
        "first2": first_k_consecutive(det_trail, 2),
        "first3": first_k_consecutive(det_trail, 3),
        "SUM": int(np.sum(det_trail)),
        "n_trail": int(len(det_trail)),
    }


# ---------------------------------------------------------------------------
# Data / features
# ---------------------------------------------------------------------------


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
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=4, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    clf.fit(X[:ref_end], y[:ref_end])
    pred = clf.predict(X).astype(int)
    return (pred != y).astype(float)


def to_batch_det(det_t: np.ndarray, *, n_per: int, n_batches: int) -> np.ndarray:
    """Observation fires → batch fires: batch True if any t in batch fires."""
    out = np.zeros(n_batches, dtype=bool)
    for b in range(n_batches):
        sl = slice(b * n_per, (b + 1) * n_per)
        out[b] = bool(np.any(det_t[sl]))
    return out


# ---------------------------------------------------------------------------
# Detectors (trail observation boolean), then batch-aggregate
# ---------------------------------------------------------------------------


def det_bocpd(mse: np.ndarray, *, burnin: int = 20) -> np.ndarray:
    mse = np.asarray(mse, dtype=float)
    mu0 = float(np.mean(mse[: min(burnin, len(mse))]))
    obs = StudentT(alpha=1.0, beta=1.0, kappa=1.0, mu=mu0)
    hazard = lambda r, h=100: constant_hazard(h, r)
    _R, maxes = online_changepoint_detection(mse, hazard, obs)
    maxes = np.asarray(maxes).ravel()
    if len(maxes) > len(mse):
        maxes = maxes[1 : len(mse) + 1]
    elif len(maxes) < len(mse):
        pad = maxes[-1] if len(maxes) else 0
        maxes = np.pad(maxes, (0, len(mse) - len(maxes)), constant_values=pad)
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
    thr = float(np.percentile(mse[: min(burnin, len(mse))], q))
    stream = (mse > thr).astype(int)
    det = Detector()
    out = []
    for v in stream:
        det.update(value=float(v))
        out.append(bool(getattr(det, "drift", False)))
    return np.asarray(out, dtype=bool)


def online_rfperm_batch_fires(
    records: list[dict], *, n_per: int, gate: float
) -> np.ndarray:
    X, y = featurize_records(records)
    n_batches = len(records) // n_per
    fired = np.zeros(n_batches, dtype=bool)
    e_prev = None
    X_prev = y_prev = None
    for b in range(n_batches):
        sl = slice(b * n_per, (b + 1) * n_per)
        Xb, yb = X[sl], y[sl]
        hit = False
        if X_prev is not None:
            probe = fit_online_probe(X_prev, y_prev, task="acc", seed=b)
            e_now = probe_err(probe, Xb, yb, task="acc")
            fl = error_floor("acc", n_per)
            hit = bool(hop_fires(e_now, e_prev, gate=gate, e_floor=fl))
            hopped_now = bool(records[b * n_per]["hopped"])
            if (
                not hit
                and hopped_now
                and e_prev is not None
                and float(e_prev) < fl
                and float(e_now) >= fl
                and float(e_now) / fl >= gate
            ):
                hit = True
            e_prev = float(e_now)
        else:
            probe = fit_online_probe(Xb, yb, task="acc", seed=b)
            e_prev = float(probe_err(probe, Xb, yb, task="acc"))
        fired[b] = hit
        X_prev, y_prev = Xb, yb
    return fired


def run_dataset(name: str, path: Path, *, gate: float) -> dict:
    records = load_stream(path)
    n_per = sum(1 for r in records if r["batch"] == 0)
    n_batches = len(records) // n_per
    cut = next(r["batch"] for r in records if r.get("hopped"))
    X, y = featurize_records(records)
    ref_end = cut * n_per

    # Full-stream MSE from quiet ref (OOB recipe); trail = after cut
    mse = mse_stream_from_ref(X, y, ref_end=ref_end, seed=0)
    mse_trail = mse[ref_end:]  # B=1 trail observations
    burnin = max(15, ref_end // 2)

    # --- batch-grid trail dets (primary tables) ---
    orf_all = online_rfperm_batch_fires(records, n_per=n_per, gate=gate)
    orf_trail = orf_all[cut:]  # trail batches only

    obs_dets = {
        "BOCPD": det_bocpd(mse, burnin=burnin),
        "PageHinkley": det_page_hinkley(mse),
        "ADWIN": det_adwin(mse),
        "DDM": det_frouros(mse, DDM, burnin=burnin),
        "STEPD": det_frouros(mse, STEPD, burnin=burnin),
        "HDDMA": det_frouros(mse, HDDMA, burnin=burnin),
        "ECDDWT": det_frouros(mse, ECDDWT, burnin=burnin),
    }

    methods_batch = [first_k_stats(orf_trail, method="OnlineRFPerm")]
    methods_obs = []  # B=1 trail observation first_k (manuscript Table 1/2 style)

    for method, det_full in obs_dets.items():
        # batch aggregate on full stream then take trail batches
        batch_all = to_batch_det(det_full, n_per=n_per, n_batches=n_batches)
        methods_batch.append(first_k_stats(batch_all[cut:], method=method))
        # observation-level on trail only
        methods_obs.append(first_k_stats(det_full[ref_end:], method=method))

    # OnlineRFPerm has no native B=1 obs path here; map batch fire → obs block
    orf_obs = np.repeat(orf_trail, n_per)
    methods_obs.insert(0, first_k_stats(orf_obs, method="OnlineRFPerm"))

    return {
        "dataset": name,
        "n": len(records),
        "n_per": n_per,
        "n_batches": n_batches,
        "cut_batch": cut,
        "cut_t": ref_end,
        "n_trail_batches": int(n_batches - cut),
        "n_trail_obs": int(len(records) - ref_end),
        "mse_mean_quiet": float(np.mean(mse[:ref_end])),
        "mse_mean_hop": float(np.mean(mse[ref_end:])),
        "batch": methods_batch,  # primary: trail batch first_k
        "obs_B1": methods_obs,  # manuscript-style trail observation first_k
        "orf_trail_fires": [bool(x) for x in orf_trail],
    }


def latex_escape(s: str) -> str:
    return str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


def _fmt(v) -> str:
    return "---" if v is None else str(int(v))


def write_latex(results: dict, out: Path) -> str:
    """Manuscript-style Table 1 (k=1) + Table 2 (k=2/3), trail *batch* index."""
    ds_names = list(results.keys())
    # map method -> {ds: row}
    methods_order = [m["method"] for m in results[ds_names[0]]["batch"]]

    def cell(ds, method, key, grid="batch"):
        row = next(r for r in results[ds][grid] if r["method"] == method)
        return _fmt(row[key])

    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=1in]{geometry}",
        r"\usepackage{booktabs,amsmath}",
        r"\title{OnlineRFPerm vs OOB baselines on LLM infer streams\\"
        r"(manuscript first$_k$ on trail batches)}",
        r"\author{CFPerm benchmark}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        "",
        r"\paragraph{Protocol.}",
        (
            r"Quiet reference freezes the probe; trail $=$ batches after cut. "
            r"\texttt{first}$_k$ $=$ first \emph{trail batch index} at which "
            r"$k$ consecutive fires begin (0 $=$ first batch after cut). "
            r"Same counting rule as PointCloud Table~1/2 "
            r"(\texttt{first\_k\_consecutive}). "
            f"HaluEval/SQuAD: $n_{{\\mathrm{{per}}}}={results[ds_names[0]]['n_per']}$, "
            f"cut batch$={results[ds_names[0]]['cut_batch']}$ "
            f"(cut $t={results[ds_names[0]]['cut_t']}$)."
        ),
        "",
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Table 1 --- First rejection ($k=1$) on the trail (batch index).}",
        r"\begin{tabular}{l" + "r" * len(ds_names) + "}",
        r"\toprule Method & " + " & ".join(latex_escape(d) for d in ds_names) + r" \\",
        r"\midrule",
    ]
    for method in methods_order:
        vals = " & ".join(cell(d, method, "first1") for d in ds_names)
        lines.append(f"{latex_escape(method)} & {vals} \\\\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Table 2 --- Sustained onsets: $k=2$ / $k=3$ consecutive trail batches.}",
        r"\begin{tabular}{l" + "rr" * len(ds_names) + "}",
        r"\toprule",
        "Method & "
        + " & ".join(rf"{latex_escape(d)} $k=2$ & {latex_escape(d)} $k=3$" for d in ds_names)
        + r" \\",
        r"\midrule",
    ]
    for method in methods_order:
        parts = []
        for d in ds_names:
            parts.append(cell(d, method, "first2"))
            parts.append(cell(d, method, "first3"))
        lines.append(f"{latex_escape(method)} & " + " & ".join(parts) + r" \\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\paragraph{Python logic.}",
        r"\begin{verbatim}",
        "def first_k_consecutive(det, k):",
        "    # det: bool array on the TRAIL only (after cut)",
        "    # return first i where det[i:i+k] are all True",
        "    det = np.asarray(det, dtype=bool)",
        "    for i in range(len(det) - k + 1):",
        "        if det[i:i+k].all():",
        "            return i   # trail batch index",
        "    return None",
        "",
        "first1 = first_k_consecutive(trail_fires, 1)",
        "first2 = first_k_consecutive(trail_fires, 2)  # 头一回连续两个 fire",
        "first3 = first_k_consecutive(trail_fires, 3)",
        r"\end{verbatim}",
        "",
        r"\paragraph{Takeaway.}",
        (
            r"Index 0 is the first hop batch. "
            r"``---'' $=$ never $k$ consecutive fires on the trail. "
            r"Baselines follow \texttt{onlinePermOOB\_algo.py} on the frozen-probe error stream."
        ),
        "",
        r"\end{document}",
        "",
    ]
    tex = "\n".join(lines)
    (out / "bench_detectors.tex").write_text(tex)
    (DOCS / "BENCH_INFER_DETECTORS.tex").write_text(tex)

    # also dump the logic snippet alone
    logic = (
        "# manuscript first_k on trail batches\n"
        "import numpy as np\n\n"
        "def first_k_consecutive(det, k):\n"
        "    \"\"\"First trail index i where det[i:i+k] are all True.\n"
        "    det is bool over trail batches only (after cut).\n"
        "    first2 = 头一回连续两个 batch fire 从第几个 trail batch 开始.\n"
        "    \"\"\"\n"
        "    det = np.asarray(det, dtype=bool).ravel()\n"
        "    for i in range(len(det) - k + 1):\n"
        "        if det[i:i+k].all():\n"
        "            return int(i)\n"
        "    return None\n\n"
        "# example\n"
        "# trail_fires = [False, True, True, False]  # batches after cut\n"
        "# first1 -> 1; first2 -> 1; first3 -> None\n"
    )
    (out / "first_k_logic.py").write_text(logic)
    (DOCS / "FIRST_K_LOGIC.py").write_text(logic)
    return tex


def export_datasets(out_data: Path) -> None:
    """Copy / preview streams so user can inspect."""
    out_data.mkdir(parents=True, exist_ok=True)
    for name in ("halueval", "squad"):
        src = STREAM_ROOT / name / "stream.jsonl"
        if not src.exists():
            continue
        dst = out_data / f"{name}_stream.jsonl"
        shutil.copy(src, dst)
        rows = [json.loads(l) for l in src.read_text().splitlines() if l.strip()]
        df = pd.DataFrame(rows)
        keep = [
            c
            for c in [
                "t",
                "batch",
                "hopped",
                "system",
                "question",
                "answer",
                "gold",
                "rag_hit",
                "faith",
                "y_bad",
            ]
            if c in df.columns
        ]
        df[keep].to_csv(out_data / f"{name}_preview.csv", index=False)
        meta = {
            "dataset": name,
            "n": len(df),
            "n_per": int((df["batch"] == 0).sum()) if "batch" in df else None,
            "cut_batch": int(df.loc[df["hopped"], "batch"].iloc[0]) if "hopped" in df else None,
            "path_jsonl": str(dst),
            "path_csv": str(out_data / f"{name}_preview.csv"),
            "source": {
                "halueval": "data/hf_cache/halueval_qa_3000.jsonl",
                "squad": "data/hf_cache/squad_val_500.jsonl",
            }.get(name),
        }
        (out_data / f"{name}_meta.json").write_text(json.dumps(meta, indent=2))
    (out_data / "README.md").write_text(
        "# Infer bench datasets (look here)\n\n"
        "- `halueval_stream.jsonl` / `halueval_preview.csv` — HaluEval live-infer stream\n"
        "- `squad_stream.jsonl` / `squad_preview.csv` — SQuAD live-infer stream\n"
        "- quiet: `hopped=false` (knowledge in prompt); hop: `hopped=true`\n"
        "- `y_bad`: faithfulness label used by detectors\n"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    if USER_ALGO.exists():
        (args.out / "onlinePermOOB_algo_source.py").write_text(USER_ALGO.read_text())

    export_datasets(DATA_EXPORT)

    specs = {
        "halueval": STREAM_ROOT / "halueval" / "stream.jsonl",
        "squad": STREAM_ROOT / "squad" / "stream.jsonl",
    }
    results = {}
    for name, path in specs.items():
        if not path.exists():
            raise SystemExit(f"missing {path}")
        print(f"[bench] {name}…", flush=True)
        results[name] = run_dataset(name, path, gate=args.gate)

    payload = {
        "stance": "manuscript first_k on trail batches (PointCloud Table 1/2)",
        "logic": "first_k_consecutive(trail_batch_fires, k)",
        "datasets_export": str(DATA_EXPORT),
        "results": results,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2))

    # markdown
    md = [
        "# Detector bench — manuscript first_k (trail batch)\n",
        "Logic: `docs/biz/FIRST_K_LOGIC.py`\n",
        f"Datasets: `{DATA_EXPORT}`\n",
        "## Table 1 — first1 (trail batch)\n",
        "| method | " + " | ".join(results.keys()) + " |",
        "|---|" + "|".join(["---:"] * len(results)) + "|",
    ]
    methods = [m["method"] for m in next(iter(results.values()))["batch"]]
    for method in methods:
        cells = []
        for ds in results:
            row = next(r for r in results[ds]["batch"] if r["method"] == method)
            cells.append(str(row["first1"]))
        md.append(f"| {method} | " + " | ".join(cells) + " |")
    md += [
        "\n## Table 2 — first2 / first3\n",
        "| method | "
        + " | ".join(f"{d} k=2 | {d} k=3" for d in results)
        + " |",
        "|---|" + "|".join(["---:"] * (2 * len(results))) + "|",
    ]
    for method in methods:
        cells = []
        for ds in results:
            row = next(r for r in results[ds]["batch"] if r["method"] == method)
            cells.append(str(row["first2"]))
            cells.append(str(row["first3"]))
        md.append(f"| {method} | " + " | ".join(cells) + " |")
    md.append("\nLaTeX: `docs/biz/BENCH_INFER_DETECTORS.tex`\n")
    (args.out / "REPORT.md").write_text("\n".join(md) + "\n")
    (DOCS / "BENCH_INFER_DETECTORS.md").write_text("\n".join(md) + "\n")

    tex = write_latex(results, args.out)
    print("\n===== LaTeX =====\n")
    print(tex)
    print(f"\n[datasets] {DATA_EXPORT}", flush=True)
    print(f"[logic] {args.out / 'first_k_logic.py'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
