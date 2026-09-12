#!/usr/bin/env python3
"""Conventional balance vs typed correlation vs effective-rank ranking.

  python3 scripts/run_typed_aux_losses.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_multimodal_attribution import write_json  # noqa: E402
from typed_aux_losses import (  # noqa: E402
    amazon_heatmap_ranks,
    plot_aux_comparison,
    plot_aux_risk,
    report_typed_streams,
    run_amazon_risk_suite,
    run_typed_risk_suite,
    write_risk_tex,
    write_tex,
)

OUT = ROOT / "results" / "typed_aux_losses"
DOCS = ROOT / "docs" / "method"
AMAZON_R = ROOT / "results" / "amazon_continuous_batches" / "amazon_batch_relationship.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per", type=int, default=48)
    ap.add_argument("--skip-risk", action="store_true")
    args = ap.parse_args()

    print("typed aux streams", flush=True)
    report = report_typed_streams(n_per=int(args.n_per), n_batches=8, seed=2026)
    amazon = None
    if AMAZON_R.exists():
        payload = json.loads(AMAZON_R.read_text())
        amazon = amazon_heatmap_ranks(
            payload["cosine"],
            labels=payload.get("labels") or payload.get("categories"),
            hops=payload.get("hops"),
        )
        print("amazon erank(R)", amazon["erank"], flush=True)
        print("amazon loo", [(r["label"], round(r["delta_erank"], 3)) for r in amazon["categories_by_loo_erank"]], flush=True)
        print("amazon hops", [(h["from"], h["to"], round(h["c_heat"], 3)) for h in amazon["hops_by_c"]], flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = dict(report)
    if amazon:
        payload = dict(report)
        payload["amazon"] = amazon
    write_json(OUT / "typed_aux_losses.json", payload)
    plot_aux_comparison(report, OUT / "typed_aux_losses.png", amazon=amazon)
    write_tex(report, amazon, OUT / "Typed_aux_losses.tex")
    write_tex(report, amazon, DOCS / "Typed_aux_losses.tex")

    for name, rec in report.items():
        print("==", name, "pair", rec["pair"], flush=True)
        print("  c", {g: round(rec["c"][g], 3) for g in rec["c"]}, flush=True)
        print("  uni CE", {g: round(rec["uni_ce"][g], 3) for g in rec["uni_ce"]}, flush=True)
        print("  inv-CE w", {g: round(rec["balance_inv"]["weights"][g], 3) for g in rec["balance_inv"]["weights"]}, flush=True)
        print("  typed L_corr", round(rec["typed_corr"]["value"], 3), flush=True)
        print("  Δerank(Cov)", {g: round(rec["erank_delta"][g], 3) for g in rec["erank_delta"]}, flush=True)
        print("  mean-Gram erank", {g: round(rec["erank_mean_gram"][g], 3) for g in rec["erank_mean_gram"]}, flush=True)
        print(
            "  rank c",
            rec["rank_c"],
            "rank inv",
            rec["rank_inv_weight"],
            "rank Δerank",
            rec["rank_erank_delta"],
            "rank mean-Gram",
            rec["rank_mean_gram"],
            flush=True,
        )
    print("wrote", OUT, flush=True)

    if not args.skip_risk:
        print("typed risk suite", flush=True)
        typed_risk = run_typed_risk_suite()
        print("amazon risk suite", flush=True)
        amazon_risk = run_amazon_risk_suite()
        payload["typed_risk"] = typed_risk
        payload["amazon_risk"] = amazon_risk
        write_json(OUT / "typed_aux_losses.json", payload)
        plot_aux_risk(typed_risk, amazon_risk, OUT / "typed_aux_risk.png")
        write_risk_tex(typed_risk, amazon_risk, OUT / "Typed_aux_risk.tex")
        write_risk_tex(typed_risk, amazon_risk, DOCS / "Typed_aux_risk.tex")
        for rname, block in typed_risk["table"].items():
            print("==", rname, flush=True)
            for aux, rec in block.items():
                print(
                    " ",
                    aux,
                    "BWT",
                    round(rec["bwt"]["mean"], 3),
                    "Brier",
                    round(rec["online_mse"]["mean"], 3),
                    "post-acc",
                    round(rec["post_acc"]["mean"], 3),
                    flush=True,
                )
        print("amazon", {k: {m: round(v["mean"], 3) for m, v in rec.items()} for k, rec in amazon_risk["table"].items()}, flush=True)


if __name__ == "__main__":
    main()
