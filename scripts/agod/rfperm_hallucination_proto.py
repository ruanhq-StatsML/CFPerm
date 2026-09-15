#!/usr/bin/env python3
"""RFPerm → LLM hallucination detection (prototype).

RFPerm does **not** prove a sentence is false. It does two things:

1. **Regime fire** — consecutive OOS error of a shallow RF probe jumps
   (e_now / e_prev ≥ γ). That means P(Y|X) for hallucination labels
   (or a proxy) shifted vs the previous batch.
2. **Instance po_risk0** — class: 1 − P_μ0(Y|X) on the new batch.
   High score = this answer breaks the *previous* probe map
   → hallucination *candidate under shift*, not a fact-check.

Mapping
-------
  batch     stream of LLM answers (time / domain / model version)
  X         answer+context features (overlap, length, self-agree, …)
  Y         hallucinated (1/0), or a cheap proxy (self-consistency fail)
  fire      hallucination *regime* alert
  po_risk0  per-answer alert rank

Not causal. Labels/proxies still required for Y.

  python3 scripts/agod/rfperm_hallucination_proto.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from agod.online_rfperm import (  # noqa: E402
    fit_online_rf,
    po_risk0_rows,
    probe_err,
    run_rfperm_stream,
)
from agod.po_refit import Stream  # noqa: E402

OUT = ROOT / "results/agod/rfperm_hallucination"
FEAT = [
    "ctx_overlap",
    "ans_len",
    "n_ents",
    "self_agree",
    "retr_hit",
    "hedge",
    "conf",
    "novel_ent",
]


def _sigmoid(z):
    z = np.clip(z, -20.0, 20.0)
    return 1.0 / (1.0 + np.exp(-z))


def make_hallucination_stream(
    *,
    n_batches: int = 8,
    n_per: int = 100,
    cut: int = 4,
    seed: int = 0,
) -> tuple[Stream, dict]:
    """Synthetic LLM stream with a sharp concept cut at ``cut``.

    Before: hallucination ≈ low ctx_overlap / low self_agree (easy probe).
    After:  high conf + novel_ent + long answers drive Y (old probe breaks).
    """
    rng = np.random.default_rng(seed)
    xs, ys, batches = [], [], []
    for b in range(n_batches):
        for _ in range(n_per):
            ctx = float(rng.uniform(0.05, 0.95))
            length = float(rng.normal(0.0, 1.0))
            n_ents = float(rng.uniform(0.0, 1.0))
            agree = float(rng.uniform(0.15, 1.0))
            retr = float(np.clip(ctx + rng.normal(0, 0.08), 0, 1))
            hedge = float(rng.uniform(0, 1))
            conf = float(rng.uniform(0.15, 1.0))
            novel = float(rng.uniform(0, 1))
            x = [ctx, length, n_ents, agree, retr, hedge, conf, novel]
            if b < cut:
                logit = (
                    5.0 * (0.5 - ctx)
                    + 3.5 * (0.5 - agree)
                    + 1.5 * (0.45 - retr)
                    + float(rng.normal(0, 0.2))
                )
            else:
                logit = (
                    5.0 * (conf - 0.5)
                    + 3.0 * (novel - 0.35)
                    + 2.5 * length
                    - 2.0 * ctx
                    + float(rng.normal(0, 0.2))
                )
            y = int(rng.random() < float(_sigmoid(logit)))
            xs.append(x)
            ys.append(y)
            batches.append(b)

    X = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=int)
    batch = np.asarray(batches, dtype=int)
    meta = {
        "feat": FEAT,
        "cut": int(cut),
        "n_batches": int(n_batches),
        "n_per": int(n_per),
        "rate_before": float(y[batch < cut].mean()),
        "rate_after": float(y[batch >= cut].mean()),
    }
    stream = Stream(
        X=X,
        y=y,
        batch=batch,
        name="llm_hallucination_synth",
        task="acc",
        meta=meta,
    )
    return stream, meta


def auroc(y_true: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y_true).astype(int)
    s = np.asarray(score, dtype=float)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = 0.0
    for p in pos:
        wins += float(np.mean(p > neg) + 0.5 * np.mean(p == neg))
    return float(wins / len(pos))


def eval_hop_ranking(stream: Stream, hop: dict) -> dict:
    if "po_risk0" not in hop:
        return {}
    idx = np.asarray(hop["index"], dtype=int)
    t1 = np.asarray(hop["treated"], dtype=int).astype(bool)
    po = np.asarray(hop["po_risk0"], dtype=float)
    rows = idx[t1]
    if rows.size == 0:
        return {}
    y = stream.y[rows].astype(int)
    scores = po[t1]
    top = y[np.argsort(-scores)][:10]
    return {
        "n_t1": int(rows.size),
        "halluc_rate": float(y.mean()),
        "auroc_po_risk0": auroc(y, scores),
        "precision_at_10": float(top.mean()) if rows.size >= 10 else float("nan"),
    }


def rank_alerts(stream: Stream, hop: dict, topk: int = 12) -> list[dict]:
    if "po_risk0" not in hop:
        return []
    idx = np.asarray(hop["index"], dtype=int)
    t1 = np.asarray(hop["treated"], dtype=int).astype(bool)
    po = np.asarray(hop["po_risk0"], dtype=float)
    local = np.flatnonzero(t1)
    if local.size == 0:
        return []
    order = local[np.argsort(-po[local])]
    out = []
    for j in order[:topk]:
        i = int(idx[j])
        out.append(
            {
                "row": i,
                "batch": int(stream.batch[i]),
                "po_risk0": float(po[j]),
                "y_true": int(stream.y[i]),
                "feat": {FEAT[k]: float(stream.X[i, k]) for k in range(len(FEAT))},
            }
        )
    return out


def run_proto(*, gate: float = 1.35, cut: int = 4, seed: int = 0) -> dict:
    """Default gate 1.35 for this synth cut; AGOD boards often use 1.5."""
    stream, meta = make_hallucination_stream(cut=cut, seed=seed)
    rec = run_rfperm_stream(stream, gate=gate, seed=seed, detail=True)
    fires = [h for h in rec["history"] if h.get("fired")]
    first_fire_t = int(fires[0]["t"]) if fires else None
    hop_at_cut = next((h for h in rec["history"] if h["t"] == cut), None)
    hop_quiet = next((h for h in rec["history"] if h["t"] == max(cut - 2, 1)), None)
    return {
        "meta": meta,
        "gate": float(gate),
        "fire_rate": float(rec["fire_rate"]),
        "online_mse": float(rec["online_mse"]),
        "first_fire_t": first_fire_t,
        "fires": [
            {
                "t": h["t"],
                "ratio": h["ratio"],
                "mean_r0": h["mean_r0"],
                "mean_r1": h["mean_r1"],
            }
            for h in fires
        ],
        "hop_at_cut": {
            "t": hop_at_cut["t"],
            "fired": hop_at_cut["fired"],
            "ratio": hop_at_cut["ratio"],
            "mean_r0": hop_at_cut["mean_r0"],
            "mean_r1": hop_at_cut["mean_r1"],
            "ranking": eval_hop_ranking(stream, hop_at_cut),
        }
        if hop_at_cut
        else None,
        "hop_quiet": {
            "t": hop_quiet["t"],
            "fired": hop_quiet["fired"],
            "ratio": hop_quiet["ratio"],
            "ranking": eval_hop_ranking(stream, hop_quiet),
        }
        if hop_quiet
        else None,
        "top_alerts_at_cut": rank_alerts(stream, hop_at_cut) if hop_at_cut else [],
        "note": (
            "fire = P(halluc|X) regime jump; po_risk0 = instance break of previous probe. "
            "Not a factuality oracle."
        ),
    }


def demo_instance_probe(*, cut: int = 4, seed: int = 1) -> dict:
    stream, meta = make_hallucination_stream(cut=cut, seed=seed)
    tr = stream.batch == (cut - 1)
    te = stream.batch == cut
    probe = fit_online_rf(stream.X[tr], stream.y[tr], seed=seed, task="acc")
    e = probe_err(probe, stream.X[te], stream.y[te], task="acc")
    po = po_risk0_rows(probe, stream.X[te], stream.y[te], task="acc")
    return {
        "trusted_batch": cut - 1,
        "new_batch": cut,
        "oos_err": float(e),
        "auroc": auroc(stream.y[te], po),
        "precision_at_10": float(stream.y[te][np.argsort(-po)][:10].mean()),
        "rate_new": float(stream.y[te].mean()),
        "meta_rates": {"before": meta["rate_before"], "after": meta["rate_after"]},
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report = run_proto(gate=1.35, cut=4, seed=0)
    inst = demo_instance_probe(cut=4, seed=1)

    print("RFPerm hallucination prototype")
    print("fire = regime jump in P(Y|X); po_risk0 = instance break of prior probe")
    print("not a fact-checker\n")
    print(
        f"cut={report['meta']['cut']}  gate={report['gate']}  "
        f"halluc_rate before/after="
        f"{report['meta']['rate_before']:.3f}/{report['meta']['rate_after']:.3f}"
    )
    print(f"fire_rate={report['fire_rate']:.3f}  first_fire_t={report['first_fire_t']}")
    print("fires", report["fires"])
    if report["hop_at_cut"]:
        h = report["hop_at_cut"]
        print(
            f"hop@cut t={h['t']} fired={h['fired']} ratio={h['ratio']:.3f} "
            f"e0={h['mean_r0']:.3f} e1={h['mean_r1']:.3f}"
        )
        print("  ranking", h["ranking"])
    if report["hop_quiet"]:
        q = report["hop_quiet"]
        print(f"hop@quiet t={q['t']} fired={q['fired']} ratio={q['ratio']:.3f}")
        print("  ranking", q["ranking"])

    print("\ntop po_risk0 alerts @ cut")
    for a in report["top_alerts_at_cut"][:8]:
        f = a["feat"]
        print(
            f"  row={a['row']} y={a['y_true']} po={a['po_risk0']:.3f} "
            f"ctx={f['ctx_overlap']:.2f} conf={f['conf']:.2f} novel={f['novel_ent']:.2f}"
        )
    print("\ninstance probe", inst)

    assert report["first_fire_t"] is not None, "expected regime fire near concept cut"
    assert report["first_fire_t"] >= report["meta"]["cut"] - 1
    auc = report["hop_at_cut"]["ranking"]["auroc_po_risk0"]
    assert auc == auc and auc > 0.55, auc
    assert inst["auroc"] > 0.55, inst["auroc"]

    (OUT / "report.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    (OUT / "instance_probe.json").write_text(json.dumps(inst, indent=2), encoding="utf-8")
    (OUT / "REPORT.md").write_text(
        "\n".join(
            [
                "# RFPerm hallucination prototype",
                "",
                "RFPerm **regime fire** = consecutive OOS jump in P(Y|X) for hallucination labels.",
                "Instance **po_risk0** = break of previous probe → candidate under shift.",
                "Not a fact-checker. Not causal.",
                "",
                f"- cut: {report['meta']['cut']}",
                f"- gate: {report['gate']}",
                f"- rate before/after: {report['meta']['rate_before']:.3f} / {report['meta']['rate_after']:.3f}",
                f"- fire_rate: {report['fire_rate']:.3f}",
                f"- first_fire_t: {report['first_fire_t']}",
                f"- hop@cut ranking: `{report['hop_at_cut']['ranking']}`",
                f"- instance probe AUROC: {inst['auroc']:.3f}",
                "",
                "Run: `python3 scripts/agod/rfperm_hallucination_proto.py`",
            ]
        ),
        encoding="utf-8",
    )
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
