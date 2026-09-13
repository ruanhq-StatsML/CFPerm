#!/usr/bin/env python3
"""PO-risk feature selection × heatmap hop-ridge (+ token/patch mask drill-down).

  python3 scripts/run_po_hop_attribution.py
  python3 scripts/run_po_hop_attribution.py --quick
  python3 scripts/run_po_hop_attribution.py --synthetic
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_multimodal_attribution import write_json  # noqa: E402
from po_hop_attribution import (  # noqa: E402
    METHODS,
    plot_po_suite,
    run_amazon_po_suite,
    run_msrvtt_po_suite,
    write_po_tex,
)

OUT = ROOT / "results" / "po_hop_attribution"
DOCS = ROOT / "docs" / "method"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=240)
    ap.add_argument("--skip-msrvtt", action="store_true")
    args = ap.parse_args()

    n_per = 80 if args.quick else int(args.n_per)
    n_seeds = 2 if args.quick else int(args.seeds)
    n_est = 20 if args.quick else 40
    source = "synthetic" if args.synthetic else "amazon"
    seeds = list(range(2026, 2026 + n_seeds))

    print("PO×hop Amazon", source, n_seeds, "seeds", n_per, "per", "n_est", n_est, flush=True)
    amazon = run_amazon_po_suite(
        seeds=seeds,
        n_per=n_per,
        n_batches=9,
        source=source,
        n_estimators=n_est,
    )

    print("PO×hop MSR-VTT", "synthetic" if args.synthetic else "live", flush=True)
    if args.skip_msrvtt:
        msrvtt = {
            "hop": {"online_mse": float("nan")},
            "hop_po_mod": {"online_mse": float("nan")},
            "mask": {"by_modality": {}},
            "lift_mod": float("nan"),
            "meta": {"source": "skipped"},
            "shares_last": {},
        }
    else:
        msrvtt = run_msrvtt_po_suite(
            n_batches=5,
            seed=2026,
            synthetic=bool(args.synthetic),
            n_estimators=n_est,
        )

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {
        "amazon": {
            "table": amazon["table"],
            "mask": amazon["mask"],
            "shapes": amazon["shapes"],
            "seeds": amazon["seeds"],
            "source": amazon["source"],
            "methods": amazon["methods"],
            "rows": amazon["rows"],
        },
        "msrvtt": {
            "hop_mse": msrvtt["hop"]["online_mse"],
            "hop_po_mod_mse": msrvtt["hop_po_mod"]["online_mse"],
            "lift_mod": msrvtt["lift_mod"],
            "mask": msrvtt["mask"],
            "meta": msrvtt["meta"],
            "shares_last": msrvtt.get("shares_last"),
            "hop_history": msrvtt["hop"].get("history"),
            "hop_po_mod_history": msrvtt["hop_po_mod"].get("history"),
        },
        "n_estimators": n_est,
    }
    write_json(OUT / "po_hop_attribution.json", payload)
    try:
        plot_po_suite(amazon, msrvtt, OUT / "po_hop_attribution.png")
    except Exception as exc:  # pragma: no cover
        print("plot failed:", exc, flush=True)
    write_po_tex(amazon, msrvtt, OUT / "PO_hop_attribution.tex")
    write_po_tex(amazon, msrvtt, DOCS / "PO_hop_attribution.tex")

    print("Amazon online MSE:", flush=True)
    for method in METHODS:
        cell = amazon["table"][method]["online_mse"]
        print("  %s  %.3f (%.3f)" % (method, cell["mean"], cell["sd"]), flush=True)
    print("Amazon token mask ΔMSE:", flush=True)
    for mode, cell in amazon["mask"].items():
        if mode == "none":
            continue
        d = cell["delta_mse"]
        print("  %s  %+.3f (%.3f)" % (mode, d["mean"], d["sd"]), flush=True)
    print(
        "MSR-VTT hop / hop+π_m  %.4f / %.4f  lift=%+.4f"
        % (msrvtt["hop"]["online_mse"], msrvtt["hop_po_mod"]["online_mse"], msrvtt["lift_mod"]),
        flush=True,
    )
    if msrvtt.get("shares_last"):
        print("  shares_last", msrvtt["shares_last"], flush=True)
    for g, cell in msrvtt["mask"].get("by_modality", {}).items():
        print(
            "  mask %s zero ΔMSE %+.4f  noise %+.4f"
            % (g, cell["zero"]["delta_mse"], cell["noise"]["delta_mse"]),
            flush=True,
        )


if __name__ == "__main__":
    main()
