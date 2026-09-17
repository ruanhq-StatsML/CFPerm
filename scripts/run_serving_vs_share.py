#!/usr/bin/env python3
"""Serving error vs share error, with latency. Hop and continuous time share S(τ)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
SCRIPTS = ROOT / "scripts"
for p in (SRC, ROOT, SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from gradual_latency import (  # noqa: E402
    consecutive,
    disagreement_cell,
    first_index,
    lead_lag,
    share_error,
)
from logo_modality import brier_or_mse, fit_serving, serving_mu  # noqa: E402
from run_gradual_concept_latency import N_REF, ONSET, run_windows  # noqa: E402
from stream_dgps import make_trimodal_gradual_concept  # noqa: E402

OUT = ROOT / "results" / "serving_vs_share"


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


def replay_from(X_ref, Y_ref, batches, hat) -> list[float]:
    """Online Brier: frozen until hat (inclusive eval), then refit on all seen.

    hat=None never updates. hat=0 updates after the first scored batch.
    Eval is always before the refit of that batch.
    """
    model, binary = fit_serving(X_ref, Y_ref, seed=2026)
    Xs, Ys = [X_ref], [Y_ref]
    out = []
    for t, b in enumerate(batches):
        out.append(brier_or_mse(b["Y"], serving_mu(model, b["X"], binary), binary))
        Xs.append(b["X"])
        Ys.append(b["Y"])
        if hat is not None and t >= int(hat):
            model, binary = fit_serving(np.vstack(Xs), np.concatenate(Ys), seed=2026)
    return out


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


def plot_joint(rows, policies, dest: Path) -> None:
    import matplotlib.pyplot as plt

    t = [r["t"] for r in rows]
    t0 = ONSET
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.4), sharex=True)
    ax = axes[0]
    ax.plot(t, [r["serve_excess"] for r in rows], color="#1f4e79", label="serving excess")
    ax.plot(t, [r["share_err"] for r in rows], color="#d48b16", label="share error  (1−π_video after onset)")
    ax.axvline(t0 - 0.5, color="#999", ls=":", lw=1)
    ax.set_ylabel("error")
    ax.set_title("same clock · serving rent vs localization pointer")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.3)
    ax = axes[1]
    colors = {
        "never": "#555",
        "always": "#2a7d4f",
        "from_serve": "#1f4e79",
        "from_share": "#d48b16",
        "from_both": "#b33",
    }
    for name, loss in policies.items():
        ax.plot(t, loss, label=name, color=colors.get(name, "#333"))
    ax.axvline(t0 - 0.5, color="#999", ls=":", lw=1)
    ax.set_xlabel("batch")
    ax.set_ylabel("online Brier")
    ax.legend(frameon=False, loc="upper left", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    run = run_windows(40, with_logo=True)
    t0 = int(run["t0"])
    n_new = int(run["n_new"])
    rows_in = run["rows"]
    brier_pre = float(np.mean([r["brier_frozen"] for r in rows_in if r["t"] < t0]))
    serve_loud = [bool(r["brier_ma"] >= 1.1 * brier_pre) for r in rows_in]
    share_loud = consecutive([(r.get("pi_po_video") or 0.0) >= 0.5 for r in rows_in], k=2)
    hat_serve = first_index(serve_loud, start=t0)
    hat_share = first_index(share_loud, start=t0)
    hat_both = first_index(
        [a and b for a, b in zip(serve_loud, share_loud)],
        start=t0,
    )

    X, Y, alpha, _, meta = make_trimodal_gradual_concept(
        n_ref=N_REF, n_new=n_new, n_batches=run["n_batches"], onset_batch=ONSET, seed=2026
    )
    n_ref = meta["n_ref"]
    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    batches = []
    for rec in rows_in:
        lo = n_ref + rec["t"] * n_new
        hi = lo + n_new
        batches.append({"X": X[lo:hi], "Y": Y[lo:hi]})

    policies = {
        "never": replay_from(X_ref, Y_ref, batches, None),
        "always": replay_from(X_ref, Y_ref, batches, 0),
        "from_serve": replay_from(X_ref, Y_ref, batches, hat_serve),
        "from_share": replay_from(X_ref, Y_ref, batches, hat_share),
        "from_both": replay_from(X_ref, Y_ref, batches, hat_both),
    }

    joint_rows = []
    for rec, s_on, h_on, never_l, always_l in zip(
        rows_in, serve_loud, share_loud, policies["never"], policies["always"]
    ):
        after = rec["t"] >= t0
        pi = float(rec.get("pi_po_video") or 0.0)
        joint_rows.append(
            {
                "t": rec["t"],
                "alpha": rec["alpha"],
                "serve_excess": float(rec["brier_frozen"] - brier_pre),
                "share_err": share_error(pi, after_onset=after),
                "pi_po_video": pi,
                "serve_loud": bool(s_on),
                "share_loud": bool(h_on),
                "cell": disagreement_cell(s_on, h_on),
                "brier_never": never_l,
                "brier_always": always_l,
            }
        )

    lag = lead_lag(hat_share, hat_serve)
    area = {}
    for name, loss in policies.items():
        area[name] = float(np.sum(np.asarray(loss[t0:]) - np.asarray(policies["always"][t0:])))

    counts = {}
    post = [r for r in joint_rows if r["t"] >= t0]
    for r in post:
        counts[r["cell"]] = counts.get(r["cell"], 0) + 1
    mean_serve_by = {}
    mean_share_by = {}
    for cell in ("quiet", "share_only", "serve_only", "both"):
        sl = [r for r in post if r["cell"] == cell]
        mean_serve_by[cell] = float(np.mean([r["serve_excess"] for r in sl])) if sl else None
        mean_share_by[cell] = float(np.mean([r["share_err"] for r in sl])) if sl else None

    plot_joint(joint_rows, policies, OUT / "serving_vs_share.png")

    lines = [
        "# Serving error vs share error",
        "",
        "Hop and continuous time read the same S(τ). The coupling that needs latency is serving rent vs the localization pointer.",
        "",
        f"t0={t0}, n_new={n_new}. hat_serve={hat_serve}, hat_share={hat_share}, hat_both={hat_both}.",
        f"lead-lag (share − serve) = {_fmt(lag['lag_batches'])} batches.",
        "",
        "## Per-batch cells",
        "",
        _md_row(["t", "α", "serve excess", "share err", "π_PO video", "serve loud", "share loud", "cell"]),
        _md_row(["---"] * 8),
    ]
    for r in joint_rows:
        lines.append(
            _md_row(
                [
                    r["t"],
                    _fmt(r["alpha"]),
                    _fmt(r["serve_excess"]),
                    _fmt(r["share_err"]),
                    _fmt(r["pi_po_video"]),
                    "yes" if r["serve_loud"] else "",
                    "yes" if r["share_loud"] else "",
                    r["cell"],
                ]
            )
        )
    lines += [
        "",
        "## Policy online Brier after onset (vs always-update)",
        "",
        _md_row(["policy", "hat", "area vs always"]),
        _md_row(["---", "---", "---"]),
        _md_row(["never", "∞", _fmt(area["never"])]),
        _md_row(["from_serve", hat_serve if hat_serve is not None else "∞", _fmt(area["from_serve"])]),
        _md_row(["from_share", hat_share if hat_share is not None else "∞", _fmt(area["from_share"])]),
        _md_row(["from_both", hat_both if hat_both is not None else "∞", _fmt(area["from_both"])]),
        _md_row(["always", 0, _fmt(area["always"])]),
        "",
        "## Post-onset mean errors by cell",
        "",
        _md_row(["cell", "n", "mean serve excess", "mean share err"]),
        _md_row(["---"] * 4),
    ]
    for cell in ("quiet", "share_only", "serve_only", "both"):
        lines.append(
            _md_row(
                [
                    cell,
                    counts.get(cell, 0),
                    _fmt(mean_serve_by[cell]),
                    _fmt(mean_share_by[cell]),
                ]
            )
        )
    (OUT / "TABLES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "summary.json").write_text(
        json.dumps(
            jsonable(
                {
                    "hat_serve": hat_serve,
                    "hat_share": hat_share,
                    "hat_both": hat_both,
                    "lag": lag,
                    "area": area,
                    "counts": counts,
                    "mean_serve_by": mean_serve_by,
                    "mean_share_by": mean_share_by,
                    "rows": joint_rows,
                    "policies": policies,
                }
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
