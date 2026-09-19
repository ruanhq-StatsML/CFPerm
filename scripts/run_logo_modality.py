#!/usr/bin/env python3
"""Trimodal LOGO probe: two-layer shares, next-batch plan, regret / onlineMSE."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from logo_modality import (  # noqa: E402
    brier_or_mse,
    cumulative_regret,
    fit_serving,
    logo_batch,
    serving_mu,
    subset_excess,
)
from stream_dgps import TRIMODAL_GROUPS, make_trimodal_stream  # noqa: E402
from streaming_po_risk import (  # noqa: E402
    ref_split_baseline,
    ref_split_mmd,
    rbf_bandwidth,
    streaming_po_and_mse,
)

OUT = ROOT / "results" / "logo_modality"
KINDS = ("concept_video", "covariate_audio", "both")
N_REF = 640
N_NEW = 160
N_BATCHES = 5
ONSET = 2


def _md_row(cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def _fmt(x, n=3):
    if x is None:
        return ""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v:
        return ""
    return f"{v:.{n}g}"


def baselines(X_ref, Y_ref, seed=2026):
    po = ref_split_baseline(X_ref, Y_ref, seed=seed)
    _, mse = streaming_po_and_mse(
        X_ref[: len(X_ref) // 2],
        Y_ref[: len(Y_ref) // 2],
        X_ref[len(X_ref) // 2 :],
        Y_ref[len(Y_ref) // 2 :],
        seed=seed,
    )
    sig = rbf_bandwidth(X_ref, seed=seed)
    mmd = ref_split_mmd(X_ref, sigma=sig, seed=seed)
    return {"po_base": po, "mse_base": mse, "mmd_base": mmd, "sigma": sig}


def run_kind(kind: str, seed: int = 2026) -> dict:
    X, Y, sl, meta = make_trimodal_stream(
        n_ref=N_REF,
        n_new=N_NEW,
        n_batches=N_BATCHES,
        onset_batch=ONSET,
        kind=kind,
        seed=seed,
    )
    n_ref, n_new = meta["n_ref"], meta["n_new"]
    X_ref, Y_ref, sl_ref = X[:n_ref], Y[:n_ref], sl[:n_ref]
    base = baselines(X_ref, Y_ref, seed=seed)
    frozen, binary = fit_serving(X_ref, Y_ref, seed=seed)
    rows = []
    X_seen, Y_seen = [X_ref], [Y_ref]
    loss_frozen, loss_full, loss_logo = [], [], []
    serve_full = frozen
    serve_logo = frozen
    for t in range(N_BATCHES):
        lo = n_ref + t * n_new
        hi = lo + n_new
        Xn, Yn, sln = X[lo:hi], Y[lo:hi], sl[lo:hi]
        rec = logo_batch(
            X_ref,
            Y_ref,
            Xn,
            Yn,
            TRIMODAL_GROUPS,
            seed=seed,
            po_base=base["po_base"],
            mse_base=base["mse_base"],
            mmd_base=base["mmd_base"],
        )
        rec["t"] = t
        rec["onset"] = t >= ONSET
        rec["subset"] = subset_excess(
            X_ref, Y_ref, Xn, Yn, sln, labels_ref=sl_ref, seed=seed, min_n=20
        )
        rec["subset_logo"] = None
        if rec["subset"]:
            hot = rec["subset"][0]
            mask = sln == hot["subset"]
            if int(mask.sum()) >= 20:
                sub = logo_batch(
                    X_ref,
                    Y_ref,
                    Xn[mask],
                    Yn[mask],
                    TRIMODAL_GROUPS,
                    seed=seed,
                )
                rec["subset_logo"] = {
                    "subset": hot["subset"],
                    "n": int(mask.sum()),
                    "excess_loss": hot["excess_loss"],
                    "pi_loss": sub["pi_loss"],
                    "pi_po": sub["pi_po"],
                    "pi_mmd": sub["pi_mmd"],
                    "ratios": sub["ratios"],
                    "dominant": sub["plan"]["dominant"],
                }
        mu_f = serving_mu(frozen, Xn, binary)
        mu_u = serving_mu(serve_full, Xn, binary)
        mu_l = serving_mu(serve_logo, Xn, binary)
        lf = brier_or_mse(Yn, mu_f, binary)
        lu = brier_or_mse(Yn, mu_u, binary)
        ll = brier_or_mse(Yn, mu_l, binary)
        rec["online_brier_frozen"] = lf
        rec["online_brier_full"] = lu
        rec["online_brier_logo"] = ll
        loss_frozen.append(lf)
        loss_full.append(lu)
        loss_logo.append(ll)
        X_seen.append(Xn)
        Y_seen.append(Yn)
        Xc = np.vstack(X_seen)
        Yc = np.concatenate(Y_seen)
        serve_full, _ = fit_serving(Xc, Yc, seed=seed)
        upd = rec["plan"]["update"]
        if upd in ("full_train", "selective_freeze", "stem_adapt"):
            serve_logo, _ = fit_serving(Xc, Yc, seed=seed)
        rows.append(rec)

    oracle = loss_full if kind == "concept_video" else loss_frozen
    regret = cumulative_regret(loss_logo, oracle)
    return {
        "kind": kind,
        "meta": meta,
        "baselines": base,
        "rows": rows,
        "loss_frozen": loss_frozen,
        "loss_full": loss_full,
        "loss_logo": loss_logo,
        "regret": regret.tolist(),
        "oracle": "full_update" if kind == "concept_video" else "frozen",
        "shifted_slice": meta["shifted_slice"],
    }


def plot_curves(results: dict, dest: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.6), sharey=True)
    for ax, kind in zip(axes, KINDS):
        r = results[kind]
        t = np.arange(len(r["loss_frozen"]))
        ax.plot(t, r["loss_frozen"], label="frozen", color="#555")
        ax.plot(t, r["loss_full"], label="full update", color="#1f4e79")
        ax.plot(t, r["loss_logo"], label="LOGO policy", color="#d48b16", ls="--")
        ax.axvline(ONSET - 0.5, color="#999", ls=":", lw=1)
        ax.set_title(kind)
        ax.set_xlabel("batch")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("online Brier")
    axes[2].legend(frameon=False, fontsize=8)
    fig.suptitle("online Brier · gray dotted = labeled onset", fontsize=11)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def plot_shares(results: dict, dest: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(9.5, 8.2))
    mods = ("video", "audio", "text")
    colors = {"video": "#b33", "audio": "#6b4ea0", "text": "#2a7d4f"}
    for i, kind in enumerate(KINDS):
        rows = results[kind]["rows"]
        t = np.arange(len(rows))
        for j, key in enumerate(("pi_loss", "pi_po")):
            ax = axes[i, j]
            for g in mods:
                ax.plot(t, [r[key][g] for r in rows], label=g, color=colors[g])
            ax.axvline(ONSET - 0.5, color="#999", ls=":", lw=1)
            ax.set_ylim(-0.05, 1.05)
            ax.set_title(f"{kind} · {key}")
            ax.grid(alpha=0.3)
        axes[i, 0].set_ylabel("share")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[-1, 0].set_xlabel("batch")
    axes[-1, 1].set_xlabel("batch")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def write_markdown(results: dict, dest: Path) -> None:
    lines = [
        "# LOGO trimodal probe",
        "",
        "Shares are localization proxies, not a unique decomposition.",
        "Layer 1 = Brier LOGO. Layer 2 = PO-risk / MMD mix. T is the batch label.",
        "",
        f"n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}, onset={ONSET}.",
        "Shifted slice is always 3.",
        "",
        "## Per-batch two-layer shares and next-batch plan",
        "",
    ]
    header = [
        "kind",
        "t",
        "onset",
        "action",
        "update",
        "Brier",
        "π_loss v/a/t",
        "π_PO v/a/t",
        "π_MMD v/a/t",
        "mix_PO v/a/t",
        "video tower",
        "audio tower",
        "text tower",
        "fusion",
        "top subset",
        "subset π_MMD v/a/t",
        "subset dominant",
    ]
    lines.append(_md_row(header))
    lines.append(_md_row(["---"] * len(header)))

    def trip(src, key):
        return "/".join(_fmt(src[g][key] if isinstance(src.get(g), dict) else src.get(g, 0)) for g in ("video", "audio", "text"))

    def share_trip(src):
        return "/".join(_fmt(src.get(g, 0)) for g in ("video", "audio", "text"))

    for kind in KINDS:
        for rec in results[kind]["rows"]:
            r = rec["ratios"]
            plan = rec["plan"]
            sub = rec.get("subset_logo") or {}
            top = rec["subset"][0]["subset"] if rec["subset"] else ""
            lines.append(
                _md_row(
                    [
                        kind,
                        rec["t"],
                        "yes" if rec["onset"] else "",
                        rec["global_action"],
                        plan["update"],
                        _fmt(rec["loss_full"]),
                        trip(r, "pi_loss"),
                        trip(r, "pi_po"),
                        trip(r, "pi_mmd"),
                        trip(r, "mix_po"),
                        plan["towers"]["video"]["tower"],
                        plan["towers"]["audio"]["tower"],
                        plan["towers"]["text"]["tower"],
                        plan["fusion"],
                        top,
                        share_trip(sub.get("pi_mmd") or {}),
                        sub.get("dominant") or "",
                    ]
                )
            )
    lines += [
        "",
        "## online Brier and cumulative regret",
        "",
        _md_row(["kind", "oracle", "t", "frozen", "full update", "LOGO policy", "regret"]),
        _md_row(["---"] * 7),
    ]
    for kind in KINDS:
        r = results[kind]
        for t, (a, b, c, g) in enumerate(
            zip(r["loss_frozen"], r["loss_full"], r["loss_logo"], r["regret"])
        ):
            lines.append(
                _md_row([kind, r["oracle"], t, _fmt(a), _fmt(b), _fmt(c), _fmt(g)])
            )
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    return obj


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = {kind: run_kind(kind) for kind in KINDS}
    plot_curves(results, OUT / "online_brier.png")
    plot_shares(results, OUT / "logo_shares.png")
    write_markdown(results, OUT / "TABLES.md")
    slim = {}
    for kind, r in results.items():
        slim[kind] = {
            "meta": r["meta"],
            "oracle": r["oracle"],
            "loss_frozen": r["loss_frozen"],
            "loss_full": r["loss_full"],
            "loss_logo": r["loss_logo"],
            "regret": r["regret"],
            "rows": [
                {
                    "t": rec["t"],
                    "onset": rec["onset"],
                    "global_action": rec["global_action"],
                    "plan": rec["plan"],
                    "ratios": rec["ratios"],
                    "loss_full": rec["loss_full"],
                    "po_full": rec["po_full"],
                    "mmd_full": rec["mmd_full"],
                    "online_brier_frozen": rec["online_brier_frozen"],
                    "online_brier_full": rec["online_brier_full"],
                    "online_brier_logo": rec["online_brier_logo"],
                    "subset": rec["subset"][:4],
                    "subset_logo": rec.get("subset_logo"),
                }
                for rec in r["rows"]
            ],
        }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
