#!/usr/bin/env python3
"""Gradual concept in continuous time: hop vs trend vs area-latency."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gradual_latency import consecutive, summarize_detector  # noqa: E402
from logo_modality import (  # noqa: E402
    brier_or_mse,
    fit_serving,
    logo_batch,
    serving_mu,
)
from online_rfperm import FrozenRFPerm  # noqa: E402
from stream_dgps import TRIMODAL_GROUPS, make_trimodal_gradual_concept  # noqa: E402
from streaming_po_risk import (  # noqa: E402
    large_deviation,
    mmd_vs_reference,
    moving_average,
    rbf_bandwidth,
    ref_split_baseline,
    ref_split_mmd,
    streaming_po_and_mse,
)

OUT = ROOT / "results" / "gradual_concept_latency"
N_REF = 480
ONSET = 3
STREAM_OBS = 640


def _fmt(x, n=3):
    if x is None:
        return "∞"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v:
        return ""
    return f"{v:.{n}g}"


def _md_row(cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    return obj


def split_batches(X, Y, alpha, n_ref, n_new, n_batches):
    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    batches = []
    for t in range(n_batches):
        lo = n_ref + t * n_new
        hi = lo + n_new
        batches.append(
            {
                "t": t,
                "tau_end": hi - 1,
                "alpha": float(np.mean(alpha[lo:hi])),
                "X": X[lo:hi],
                "Y": Y[lo:hi],
            }
        )
    return X_ref, Y_ref, batches


def run_windows(n_new: int, *, with_logo: bool, seed: int = 2026) -> dict:
    n_batches = max(4, STREAM_OBS // int(n_new))
    X, Y, alpha, tau, meta = make_trimodal_gradual_concept(
        n_ref=N_REF,
        n_new=n_new,
        n_batches=n_batches,
        onset_batch=ONSET,
        seed=seed,
    )
    n_ref = meta["n_ref"]
    X_ref, Y_ref, batches = split_batches(X, Y, alpha, n_ref, n_new, n_batches)
    t0 = int(meta["onset_batch"])
    po_base = ref_split_baseline(X_ref, Y_ref, seed=seed)
    _, mse_base = streaming_po_and_mse(
        X_ref[: len(X_ref) // 2],
        Y_ref[: len(Y_ref) // 2],
        X_ref[len(X_ref) // 2 :],
        Y_ref[len(Y_ref) // 2 :],
        seed=seed,
    )
    sig = rbf_bandwidth(X_ref, seed=seed)
    mmd_base = ref_split_mmd(X_ref, sigma=sig, seed=seed)
    frozen, binary = fit_serving(X_ref, Y_ref, seed=seed)
    probe = FrozenRFPerm(X_ref, Y_ref, seed=seed)
    brier_ref = brier_or_mse(Y_ref, serving_mu(frozen, X_ref, binary), binary)

    serve_full = frozen
    rows = []
    loss_frozen, loss_full = [], []
    logo_ms = []
    X_seen, Y_seen = [X_ref], [Y_ref]
    for b in batches:
        t = b["t"]
        Xn, Yn = b["X"], b["Y"]
        rf = probe.step(Xn, Yn)
        po, mse = streaming_po_and_mse(X_ref, Y_ref, Xn, Yn, seed=seed)
        mmd = mmd_vs_reference(X_ref, Xn, sigma=sig, seed=seed)
        rec = {
            "t": t,
            "tau_end": b["tau_end"],
            "alpha": b["alpha"],
            "n_new": n_new,
            "po_stream": float(po),
            "mse_stream": float(mse),
            "mmd_vs_ref": mmd,
            "brier_frozen": brier_or_mse(Yn, serving_mu(frozen, Xn, binary), binary),
            "brier_full": brier_or_mse(Yn, serving_mu(serve_full, Xn, binary), binary),
            **rf,
            "pi_po_video": None,
            "mix_po_video": None,
            "global_action": None,
        }
        if with_logo:
            t0s = time.perf_counter()
            logo = logo_batch(
                X_ref,
                Y_ref,
                Xn,
                Yn,
                TRIMODAL_GROUPS,
                seed=seed,
                po_base=po_base,
                mse_base=mse_base,
                mmd_base=mmd_base,
            )
            logo_ms.append(1000.0 * (time.perf_counter() - t0s))
            rec["pi_po_video"] = float(logo["pi_po"]["video"])
            rec["pi_mmd_video"] = float(logo["pi_mmd"]["video"])
            rec["mix_po_video"] = float(logo["ratios"]["video"]["mix_po"])
            rec["global_action"] = logo["global_action"]
            rec["update"] = logo["plan"]["update"]
            rec["video_tower"] = logo["plan"]["towers"]["video"]["tower"]
        loss_frozen.append(rec["brier_frozen"])
        loss_full.append(rec["brier_full"])
        X_seen.append(Xn)
        Y_seen.append(Yn)
        serve_full, _ = fit_serving(np.vstack(X_seen), np.concatenate(Y_seen), seed=seed)
        rows.append(rec)

    po_ma = moving_average([r["po_stream"] for r in rows], window=max(3, 8 // max(n_new // 20, 1)))
    t_ma = moving_average([r["rfperm_T"] for r in rows], window=max(3, 8 // max(n_new // 20, 1)))
    brier_ma = moving_average(loss_frozen, window=max(3, 8 // max(n_new // 20, 1)))
    for r, p, tm, bm in zip(rows, po_ma, t_ma, brier_ma):
        r["po_ma"] = float(p)
        r["t_ma"] = float(tm)
        r["brier_ma"] = float(bm)

    pre_T = [r["rfperm_T"] for r in rows if r["t"] < t0]
    t_level = max(float(np.median(pre_T)) + 0.02, 0.02) if pre_T else 0.02
    brier_pre = float(np.mean([loss_frozen[i] for i in range(t0)])) if t0 else brier_ref

    hop = [bool(r["rfperm_hop"]) for r in rows]
    trend = [bool(r["t_ma"] >= t_level) for r in rows]
    po2 = [bool(large_deviation(r["po_ma"], po_base, ratio=2.0)) for r in rows]
    brier_walk = [bool(r["brier_ma"] >= 1.1 * brier_pre) for r in rows]
    share = [False] * len(rows)
    watch = [False] * len(rows)
    if with_logo:
        share = consecutive(
            [(r.get("pi_po_video") or 0.0) >= 0.5 for r in rows],
            k=2,
        )
        watch = [r.get("global_action") not in (None, "keep_training") for r in rows]

    detectors = [
        summarize_detector("last-two hop 1.5×", hop, t0=t0, n_new=n_new, alpha_batch=[r["alpha"] for r in rows], loss=loss_frozen, oracle=loss_full),
        summarize_detector("T MA vs pre-onset level", trend, t0=t0, n_new=n_new, alpha_batch=[r["alpha"] for r in rows], loss=loss_frozen, oracle=loss_full),
        summarize_detector("PO MA 2× baseline", po2, t0=t0, n_new=n_new, alpha_batch=[r["alpha"] for r in rows], loss=loss_frozen, oracle=loss_full),
        summarize_detector("Brier MA +10% vs pre", brier_walk, t0=t0, n_new=n_new, alpha_batch=[r["alpha"] for r in rows], loss=loss_frozen, oracle=loss_full),
    ]
    if with_logo:
        detectors += [
            summarize_detector("π_PO(video)≥0.5 twice", share, t0=t0, n_new=n_new, alpha_batch=[r["alpha"] for r in rows], loss=loss_frozen, oracle=loss_full),
            summarize_detector("board leaves keep", watch, t0=t0, n_new=n_new, alpha_batch=[r["alpha"] for r in rows], loss=loss_frozen, oracle=loss_full),
        ]

    compute = {
        "logo_ms_median": float(np.median(logo_ms)) if logo_ms else None,
        "logo_ms_mean": float(np.mean(logo_ms)) if logo_ms else None,
        "obs_per_call": int(n_new),
        "compute_lag_in_obs_if_1ms_per_obs": None,
    }
    if logo_ms:
        # If one observation arrives per 1 abstract time unit, how many obs arrive during one LOGO call?
        # Unknown wall arrival rate: report ms / n_new as cost per observation in the window.
        compute["ms_per_obs_in_window"] = float(np.median(logo_ms) / max(n_new, 1))

    return {
        "meta": meta,
        "t0": t0,
        "n_new": n_new,
        "n_batches": n_batches,
        "po_base": float(po_base),
        "mse_base": float(mse_base),
        "mmd_base": float(mmd_base),
        "brier_ref": float(brier_ref),
        "t_level": float(t_level),
        "rows": rows,
        "loss_frozen": loss_frozen,
        "loss_full": loss_full,
        "detectors": detectors,
        "compute": compute,
        "with_logo": with_logo,
    }


def plot_main(result: dict, dest: Path) -> None:
    import matplotlib.pyplot as plt

    rows = result["rows"]
    t = np.array([r["t"] for r in rows], dtype=float)
    t0 = result["t0"]
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 8.2), sharex=True)
    axes[0].plot(t, [r["alpha"] for r in rows], color="#555", label="α (true rotation)")
    axes[0].axvline(t0 - 0.5, color="#999", ls=":", lw=1)
    axes[0].set_ylabel("α")
    axes[0].legend(frameon=False, loc="upper left")
    axes[0].set_title("gradual concept · no jump")

    axes[1].plot(t, [r["rfperm_T"] for r in rows], color="#1f4e79", label="FrozenRF T")
    axes[1].plot(t, [r["t_ma"] for r in rows], color="#1f4e79", ls="--", label="T MA")
    hop_t = [r["t"] for r in rows if r.get("rfperm_hop")]
    if hop_t:
        axes[1].scatter(hop_t, [rows[i]["rfperm_T"] for i in hop_t], color="#b33", zorder=3, label="hop")
    axes[1].axvline(t0 - 0.5, color="#999", ls=":", lw=1)
    axes[1].set_ylabel("T")
    axes[1].legend(frameon=False, loc="upper left")

    axes[2].plot(t, result["loss_frozen"], color="#555", label="online Brier frozen")
    axes[2].plot(t, result["loss_full"], color="#2a7d4f", label="full update")
    if rows[0].get("pi_po_video") is not None:
        axes[2].plot(t, [r["pi_po_video"] or 0 for r in rows], color="#d48b16", label="π_PO video")
    axes[2].axvline(t0 - 0.5, color="#999", ls=":", lw=1)
    axes[2].set_xlabel("batch (window of n_new rows ≈ continuous time)")
    axes[2].set_ylabel("Brier / share")
    axes[2].legend(frameon=False, loc="upper left")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def write_markdown(main: dict, windows: list[dict], dest: Path) -> None:
    lines = [
        "# Gradual concept · continuous-time latency",
        "",
        "No jump. Last-two hop is the wrong WHEN. Evaluable latency is excess Brier until a detector fires.",
        "",
        f"Main window n_new={main['n_new']}, onset batch t0={main['t0']}, n_ref={N_REF}.",
        "",
        "## Detectors on the main gradual walk",
        "",
        _md_row(["detector", "hat", "delay batches", "delay obs", "α at hat", "FA pre", "area until hat", "never"]),
        _md_row(["---"] * 8),
    ]
    for d in main["detectors"]:
        lines.append(
            _md_row(
                [
                    d["detector"],
                    "∞" if d["never"] else d["hat_batch"],
                    _fmt(d["delay_batches"]),
                    _fmt(d["delay_obs"]),
                    _fmt(d["alpha_at_hat"]),
                    d["false_alarms_pre"],
                    _fmt(d["area_until_hat"]),
                    "yes" if d["never"] else "",
                ]
            )
        )
    lines += [
        "",
        "## Window size = statistical vs batching latency",
        "",
        _md_row(["n_new", "batches", "hop delay obs", "T-MA delay obs", "Brier+10% delay obs", "PO 2× delay obs", "hop FA"]),
        _md_row(["---"] * 7),
    ]
    for w in windows:
        by = {d["detector"]: d for d in w["detectors"]}
        lines.append(
            _md_row(
                [
                    w["n_new"],
                    w["n_batches"],
                    _fmt(by["last-two hop 1.5×"]["delay_obs"]),
                    _fmt(by["T MA vs pre-onset level"]["delay_obs"]),
                    _fmt(by["Brier MA +10% vs pre"]["delay_obs"]),
                    _fmt(by["PO MA 2× baseline"]["delay_obs"]),
                    by["last-two hop 1.5×"]["false_alarms_pre"],
                ]
            )
        )
    if main["compute"].get("logo_ms_median") is not None:
        lines += [
            "",
            f"LOGO wall time median {main['compute']['logo_ms_median']:.0f} ms / window of {main['n_new']} rows "
            f"({main['compute']['ms_per_obs_in_window']:.1f} ms per row in the window).",
        ]
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    main_run = run_windows(40, with_logo=True)
    windows = [run_windows(n, with_logo=False, seed=2026) for n in (20, 40, 80)]
    plot_main(main_run, OUT / "gradual_clocks.png")
    write_markdown(main_run, windows, OUT / "TABLES.md")
    slim = {
        "main": {
            "meta": main_run["meta"],
            "t0": main_run["t0"],
            "detectors": main_run["detectors"],
            "compute": main_run["compute"],
            "loss_frozen": main_run["loss_frozen"],
            "loss_full": main_run["loss_full"],
            "rows": [
                {
                    k: rec[k]
                    for k in (
                        "t",
                        "tau_end",
                        "alpha",
                        "rfperm_T",
                        "rfperm_hop",
                        "t_ma",
                        "po_stream",
                        "po_ma",
                        "mmd_vs_ref",
                        "brier_frozen",
                        "brier_full",
                        "pi_po_video",
                        "mix_po_video",
                        "global_action",
                        "update",
                        "video_tower",
                    )
                    if k in rec
                }
                for rec in main_run["rows"]
            ],
        },
        "windows": [
            {"n_new": w["n_new"], "n_batches": w["n_batches"], "detectors": w["detectors"]}
            for w in windows
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
