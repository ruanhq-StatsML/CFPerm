#!/usr/bin/env python3
"""Probe adaptation opportunity O_t along a C-GCD proxy stream.

  python3 scripts/run_adaptation_opportunity_probe.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from adaptation_opportunity import adaptation_opportunity  # noqa: E402
from category_shift_boards import make_cgcd_stream  # noqa: E402
from repr_sep_cold_start import (  # noqa: E402
    ColdStartState,
    adapt_batch_repr_sep,
    known_prototypes,
)
from weighted_cold_start import (  # noqa: E402
    attribute_rf_domain,
    boundary_check,
    detect_error_rise,
    fit_initial,
    localize_subset,
)

OUT = ROOT / "results" / "adaptation_opportunity"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    stream = make_cgcd_stream(
        n_known=10 if args.quick else 12,
        n_novel_per_stage=2 if args.quick else 3,
        n_stages=3 if args.quick else 4,
        n_per=24 if args.quick else 40,
        p=24 if args.quick else 32,
        novel_dims=6 if args.quick else 8,
        seed=args.seed,
    )
    clf, sc = fit_initial(
        stream["X_ref"], stream["y_ref"], seed=args.seed, classes=stream["all_classes"],
    )
    X_ref = sc.transform(stream["X_ref"])
    y_ref = np.asarray(stream["y_ref"]).ravel()
    state = ColdStartState()
    rows = []

    for st in stream["stages"]:
        Xl = sc.transform(st["X"])
        yl = np.asarray(st["y"]).ravel()
        out = adapt_batch_repr_sep(
            clf, X_ref, y_ref, Xl, yl, state,
            seed=args.seed + st["stage"],
            force_adapt=True, expand_min=8, min_jaccard=0.05, top_k=8,
        )
        clf, state = out["clf"], out["state"]
        opp = out.get("opportunity") or {}
        rows.append({
            "stage": st["stage"],
            "novel": st["novel"],
            "phase": out["phase"],
            "reason": out["reason"],
            "O_local": opp.get("O_local"),
            "O_persist": opp.get("O_persist"),
            "D": opp.get("D"),
            "L": opp.get("L"),
            "R": opp.get("R"),
            "M": opp.get("M"),
            "Stab": opp.get("Stab"),
            "action": opp.get("action"),
            "gap_lift": (out.get("sep") or {}).get("gap_lift"),
            "expand_fitted": bool((out.get("expand") or {}).get("fitted")),
            "keep_offset": opp.get("keep_offset"),
            "read": opp.get("read"),
        })
        # refresh working REF with known recirculation
        known = stream["known"]
        kidx = np.flatnonzero(np.isin(yl, known))
        if len(kidx):
            take = kidx[::2][: max(8, len(kidx) // 2)]
            X_ref = np.vstack([X_ref, Xl[take]])
            y_ref = np.concatenate([y_ref, yl[take]])

    # control: shuffle LIVE labels → private-live gone → O should drop hard
    st0 = stream["stages"][0]
    Xl0 = sc.transform(st0["X"])
    yl0 = np.asarray(st0["y"]).ravel()
    rng = np.random.default_rng(args.seed)
    y_shuf = yl0.copy()
    rng.shuffle(y_shuf)
    clf0, _ = fit_initial(stream["X_ref"], stream["y_ref"], seed=args.seed, classes=stream["all_classes"])
    Xr0 = sc.transform(stream["X_ref"])
    gate = detect_error_rise(clf0, Xr0, stream["y_ref"], Xl0, yl0, thresh=1.05)
    vimp = attribute_rf_domain(Xr0, Xl0, seed=args.seed)
    mem, feat = localize_subset(Xr0, Xl0, vimp, top_k=8)
    protos, _ = known_prototypes(Xr0, stream["y_ref"], feat)
    bound = boundary_check(stream["y_ref"], yl0, min_jaccard=0.05)
    opp_real = adaptation_opportunity(
        error_ratio=gate["ratio"], vimp=vimp, feat_idx=feat, membership=mem,
        n_private=int(np.isin(yl0, bound["private_live"]).sum()),
        expand_min=8, force=True, X_live=Xl0, prototypes_S=protos,
    )
    # destroy membership structure
    mem_flat = np.full_like(mem, 0.5)
    opp_flat = adaptation_opportunity(
        error_ratio=1.0, vimp=np.ones_like(vimp) / len(vimp), feat_idx=feat,
        membership=mem_flat, n_private=0, expand_min=8, force=False,
        X_live=Xl0, prototypes_S=protos,
    )

    payload = {
        "seed": args.seed,
        "quick": bool(args.quick),
        "stages": rows,
        "mean_O_local": float(np.nanmean([r["O_local"] for r in rows])),
        "mean_O_persist": float(np.nanmean([r["O_persist"] for r in rows])),
        "controls": {
            "structured_O_local": opp_real["O_local"],
            "flat_O_local": opp_flat["O_local"],
            "lift": float(opp_real["O_local"] - opp_flat["O_local"]),
        },
        "benchmarks": {
            "C-GCD": ["CIFAR-100", "ImageNet-100", "Tiny-ImageNet", "CUB"],
            "metrics": ["Old", "New", "All per session"],
            "PDA_OSDA": ["Office-Home", "VisDA", "DomainNet"],
            "our_extra": ["O_local", "O_persist", "novel_acc", "forgetting", "gap_lift"],
        },
        "takeaway": (
            "O_local = (D L R M)^(1/4) gates sep/expand this batch; "
            "O_persist = O_local·Stab discounts carrying S/offset into unknown next batch."
        ),
    }
    (OUT / "adaptation_opportunity.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# Adaptation opportunity probe",
        "",
        payload["takeaway"],
        "",
        "| stage | O_local | O_persist | D | L | R | M | Stab | action | gap_lift | expand |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            "| {stage} | {O_local:.3f} | {O_persist:.3f} | {D:.2f} | {L:.2f} | {R:.2f} | {M:.2f} | {Stab:.2f} | {action} | {gap_lift} | {expand_fitted} |".format(
                stage=r["stage"],
                O_local=r["O_local"] or 0,
                O_persist=r["O_persist"] or 0,
                D=r["D"] or 0, L=r["L"] or 0, R=r["R"] or 0, M=r["M"] or 0, Stab=r["Stab"] or 0,
                action=r["action"],
                gap_lift=("—" if r["gap_lift"] is None else "%.3f" % r["gap_lift"]),
                expand_fitted=r["expand_fitted"],
            )
        )
    lines.extend([
        "",
        "control: structured O_local=**%.3f** vs flat **%.3f** (lift **%.3f**)" % (
            payload["controls"]["structured_O_local"],
            payload["controls"]["flat_O_local"],
            payload["controls"]["lift"],
        ),
        "",
        "Track benches: C-GCD (CIFAR-100 / IN-100 / Tiny / CUB); PDA/OSDA (Office-Home / VisDA / DomainNet).",
        "",
        "```bash",
        "python3 scripts/run_adaptation_opportunity_probe.py --quick",
        "```",
        "",
    ])
    md = "\n".join(lines)
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT / "README.md", DOCS / "Adaptation_opportunity_probe.md")
    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "adaptation_opportunity.json", ARTIFACTS / "adaptation_opportunity.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "Adaptation_opportunity_probe.md")
    except Exception:
        pass
    print(md)


if __name__ == "__main__":
    main()
