#!/usr/bin/env python3
"""Multi-group attribution (CFPerm skeleton) on live two-stream tables.

This is not the serving hop. The hop gates ask whether P(Y|X) moved in time.
This object asks which X_j carry a group difference:

    T = group  (queue / channel / hop bucket)
    permute T, not X
    importance = how much group-to-group P(Y|X) hangs on X_j
    two-level threshold (feature-wise, then across-feature), reject iff ≥ top_k hits

Live tables already on disk:

    HH two-stream          T = helpful vs harmless queue
    real two-stream        T = BeaverTails vs ToxicChat judge
    Hotpot two-stream      T = sparse-usable vs dense-usable
    HH multi-step          T = hop 0 / hop 1 / hop 2+   (K=3)

Y is still the product label (pass/fail, channel usable). HH chosen is not Y.

The probe here is a linear group×X interaction (lstsq), so the permutation
loop is the CFPerm decision and you can see the ranking. Full permuCATE/LOCO/GRF
is the same loop with a heavier vimp.

Usage::

    PYTHONPATH=. python3 scripts/prototype_group_attribution.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.posthoc_localization import localize  # noqa: E402

AUDIT = ROOT / "results" / "manuscript" / "llm_audit"
HYBRID = ROOT / "results" / "manuscript" / "hybrid_retrieval"
OUT = ROOT / "results" / "manuscript" / "group_attribution"
LEAK = ("chosen", "rejected")
N_PERM = 40
LEVEL_FEATURE = 0.05
LEVEL_ACROSS = 0.05
TOP_K = 1
SEED = 2026


def x_columns(fieldnames) -> list[str]:
    return [c for c in fieldnames if c.startswith("x_")]


def zscore(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - mu) / sd


def interaction_vimp(X: np.ndarray, Y: np.ndarray, T: np.ndarray) -> np.ndarray:
    """|group × X| interaction energy. Reference group is the smallest label."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    T = np.asarray(T, dtype=int).ravel()
    n, p = X.shape
    groups = np.unique(T)
    if len(groups) < 2:
        return np.zeros(p)
    dummies = np.column_stack([(T == g).astype(float) for g in groups[1:]])
    k = dummies.shape[1]
    interact = (X[:, :, None] * dummies[:, None, :]).reshape(n, p * k)
    design = np.column_stack([np.ones(n), X, dummies, interact])
    beta, *_ = np.linalg.lstsq(design, Y, rcond=None)
    offset = 1 + p + k
    coef = np.asarray(beta[offset:], dtype=float).reshape(p, k)
    return np.sum(coef * coef, axis=1)


def cfperm_groups(
    X,
    Y,
    T,
    *,
    n_perm: int = N_PERM,
    level_feature: float = LEVEL_FEATURE,
    level_across_feature: float = LEVEL_ACROSS,
    top_k: int = TOP_K,
    seed: int = SEED,
    names: list[str] | None = None,
):
    """CFPerm decision: permute group labels T. Matches R/CFPerm_vimp.R."""
    X = zscore(np.asarray(X, dtype=float))
    Y = np.asarray(Y, dtype=float).ravel()
    T = np.asarray(T, dtype=int).ravel()
    p = X.shape[1]
    names = list(names) if names is not None else [f"x{j}" for j in range(p)]
    rng = np.random.default_rng(seed)
    imp = interaction_vimp(X, Y, T)
    perm = np.zeros((p, n_perm), dtype=float)
    for b in range(n_perm):
        perm[:, b] = interaction_vimp(X, Y, rng.permutation(T))
    pvals = (1.0 + np.sum(perm >= imp[:, None], axis=1)) / (1.0 + n_perm)
    q1_upper = np.quantile(perm, 1.0 - level_feature, axis=1)
    threshold = float(np.quantile(q1_upper, 1.0 - level_across_feature))
    hits = imp > threshold
    order = np.argsort(-imp)
    return {
        "names": names,
        "imp": imp,
        "pvals": pvals,
        "q1_upper": q1_upper,
        "threshold": threshold,
        "hits": hits.astype(int),
        "rejected": int(int(hits.sum()) >= int(top_k)),
        "n_hits": int(hits.sum()),
        "top": [names[int(j)] for j in order[: min(5, p)]],
        "n": int(len(Y)),
        "n_groups": int(len(np.unique(T))),
        "group_counts": {str(int(g)): int(np.sum(T == g)) for g in np.unique(T)},
        "y_by_group": {str(int(g)): float(np.mean(Y[T == g])) for g in np.unique(T)},
        "n_perm": int(n_perm),
        "level_feature": float(level_feature),
        "level_across_feature": float(level_across_feature),
        "top_k": int(top_k),
    }


def load_two_stream(path: Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty {path}")
    keys = {k.lower() for k in rows[0]}
    for bad in LEAK:
        if bad in keys:
            raise RuntimeError(f"{path} still has {bad} — that is not Y")
    cols = x_columns(rows[0].keys())
    T = np.asarray([int(r["T"]) for r in rows], dtype=int)
    Y = np.asarray([int(float(r["Y"])) for r in rows], dtype=int)
    X = np.asarray([[float(r[c]) for c in cols] for r in rows], dtype=float)
    return X, Y, T, cols


def load_multistep_groups(path: Path):
    """K=3: first hop / second hop / later hops. Same Y=pass/fail on this hop."""
    with path.open() as f:
        rows = list(csv.DictReader(f))
    keys = {k.lower() for k in rows[0]}
    for bad in LEAK:
        if bad in keys:
            raise RuntimeError(f"{path} still has {bad}")
    cols = x_columns(rows[0].keys())
    steps = np.asarray([int(r["step"]) for r in rows], dtype=int)
    T = np.clip(steps, 0, 2)
    Y = np.asarray([int(r["y"]) for r in rows], dtype=int)
    X = np.asarray([[float(r[c]) for c in cols] for r in rows], dtype=float)
    return X, Y, T, cols


def planted_dgp(n=400, p=6, seed=0):
    """Three groups. Only group 2 has a CATE on feature 1."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    T = rng.integers(0, 3, size=n)
    Y = 0.3 * X[:, 0] + 2.4 * (T == 2).astype(float) * X[:, 1] + 0.25 * rng.normal(size=n)
    names = [f"x{j}" for j in range(p)]
    return X, Y, T, names


def null_dgp(n=400, p=6, seed=1):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    T = rng.integers(0, 3, size=n)
    Y = 0.8 * X[:, 0] - 0.4 * X[:, 2] + 0.25 * rng.normal(size=n)
    names = [f"x{j}" for j in range(p)]
    return X, Y, T, names


STREAMS = [
    {
        "name": "hh_queues",
        "title": "HH helpful vs harmless",
        "t_meaning": "T=0 helpful queue, T=1 harmless queue",
        "y_meaning": "auditor pass/fail (not chosen)",
        "kind": "two_stream",
        "path": AUDIT / "xy_hh_two_stream.csv",
    },
    {
        "name": "real_judges",
        "title": "BeaverTails vs ToxicChat",
        "t_meaning": "T=0 BeaverTails is_safe, T=1 ToxicChat human",
        "y_meaning": "real auditor pass/fail",
        "kind": "two_stream",
        "path": AUDIT / "xy_real_two_stream.csv",
    },
    {
        "name": "hybrid_channels",
        "title": "sparse vs dense usable",
        "t_meaning": "T=0 sparse-only usable, T=1 dense-only usable",
        "y_meaning": "that channel covers gold titles",
        "kind": "two_stream",
        "path": HYBRID / "xy_hotpot_two_stream.csv",
    },
    {
        "name": "audit_hop_buckets",
        "title": "multi-step hop 0 / 1 / 2+",
        "t_meaning": "T=0 first hop, T=1 second, T=2 later",
        "y_meaning": "this hop pass/fail",
        "kind": "multistep",
        "path": AUDIT / "xy_hh_multistep_consistent.csv",
    },
]


def compact(rec: dict, spec: dict) -> dict:
    loc = rec.get("localization") or {}
    groups = loc.get("groups") or {}
    return {
        "name": spec["name"],
        "title": spec["title"],
        "t_meaning": spec["t_meaning"],
        "y_meaning": spec["y_meaning"],
        "n": rec["n"],
        "n_groups": rec["n_groups"],
        "group_counts": rec["group_counts"],
        "y_by_group": rec["y_by_group"],
        "rejected": rec["rejected"],
        "n_hits": rec["n_hits"],
        "threshold": rec["threshold"],
        "top": rec["top"],
        "hits": [n for n, h in zip(rec["names"], rec["hits"]) if h],
        "path": str(spec["path"].relative_to(ROOT)) if spec.get("path") else None,
        "loc_sig_pairs": groups.get("n_sig_pairs"),
        "loc_pairwise": groups.get("pairwise"),
    }


def plot_importance(runs: list[tuple[dict, dict]], path: Path) -> None:
    n = len(runs)
    fig, axes = plt.subplots(n, 1, figsize=(10.0, 2.4 * n), sharex=False)
    if n == 1:
        axes = [axes]
    for ax, (rec, spec) in zip(axes, runs):
        names = rec["names"]
        imp = rec["imp"]
        colors = ["#9b2c2c" if h else "#1f4e79" for h in rec["hits"]]
        ax.bar(range(len(names)), imp, color=colors)
        ax.axhline(rec["threshold"], color="#9b2c2c", ls="--", lw=1.0, label="across-feature threshold")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=50, ha="right", fontsize=8)
        flag = "REJECT" if rec["rejected"] else "quiet"
        ax.set_title(f"{spec['title']}  ·  {flag}  ·  top {', '.join(rec['top'][:3])}", fontsize=10)
        ax.set_ylabel("imp")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def render_report(rows: list[dict], extras: list[dict]) -> str:
    lines = [
        "# Group attribution — permute T, not X",
        "",
        "Serving gates ask whether P(Y|X) hopped in time. This is the other object:",
        "which X_j carry a **group** difference. T is the group. CFPerm permutes T.",
        "The 13 audit x_* are attribution coordinates for the group contrast, not for a time hop.",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/prototype_group_attribution.py",
        "```",
        "",
        "Decision (same as `R/CFPerm_vimp.R`): fit group x X importance, permute T, B times,",
        "feature p-value = fraction of nulls ≥ observed, across-feature threshold from the",
        f"{1-LEVEL_FEATURE:.2f} quantile of those upper tails, reject iff ≥ {TOP_K} features clear it.",
        "",
        "| Stream | T | groups | reject | hits | top |",
        "|---|---|---:|---|---|---|",
    ]
    for r in rows:
        hits = ", ".join(r["hits"]) if r["hits"] else "—"
        top = ", ".join(r["top"][:3])
        lines.append(
            f"| {r['title']} | {r['t_meaning']} | {r['n_groups']} | "
            f"{'yes' if r['rejected'] else 'no'} | {hits} | {top} |"
        )
    lines += [
        "",
        "Post-hoc localization: pull subset indices (T groups, quartiles of the top \(X_j\)),",
        "then pairwise **MMD** and **PO-risk**. Conditional means stay in the JSON; they are not the test.",
        "",
        "## Pairwise subset MMD / PO-risk",
        "",
        "| Stream | pair | n | MMD | MMD p | PO-risk | PO p | mean Y |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        for p in r.get("loc_pairwise") or []:
            lines.append(
                "| {title} | {a} vs {b} | {na}/{nb} | {mmd:.3g} | {mp:.3g} | {po:.3g} | {pp:.3g} | {ya:.3f}/{yb:.3f} |".format(
                    title=r["title"],
                    a=p["a"],
                    b=p["b"],
                    na=p["n_a"],
                    nb=p["n_b"],
                    mmd=p["mmd"],
                    mp=p["mmd_p"],
                    po=p.get("po_risk", float("nan")),
                    pp=p.get("po_p", float("nan")),
                    ya=p.get("mean_Y_a", float("nan")),
                    yb=p.get("mean_Y_b", float("nan")),
                )
            )
    lines += [
        "",
        "## Synthetic check",
        "",
        "| DGP | reject | top |",
        "|---|---|---|",
    ]
    for e in extras:
        lines.append(
            f"| {e['title']} | {'yes' if e['rejected'] else 'no'} | {', '.join(e['top'][:3])} |"
        )
    lines += [
        "",
        "Planted: three groups, only group 2 depends on `x1`. Null: Y depends on X, not on T.",
        "",
        "Not this object: refusal rate, HH chosen, episode success, single-channel Recall.",
        "",
    ]
    return "\n".join(lines) + "\n"


def run_stream(spec: dict) -> dict:
    if spec["kind"] == "two_stream":
        X, Y, T, cols = load_two_stream(spec["path"])
    elif spec["kind"] == "multistep":
        X, Y, T, cols = load_multistep_groups(spec["path"])
    else:
        raise ValueError(spec["kind"])
    rec = cfperm_groups(X, Y, T, names=cols, seed=SEED)
    rec["name"] = spec["name"]
    rec["title"] = spec["title"]
    top_name = rec["top"][0] if rec["top"] else None
    feat = X[:, cols.index(top_name)] if top_name in cols else None
    rec["localization"] = localize(X, Y, group_labels=T, feature=feat, n_perm=25, seed=SEED)
    return rec


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    runs = []
    skipped = []
    for spec in STREAMS:
        if not spec["path"].exists():
            skipped.append(str(spec["path"].relative_to(ROOT)))
            continue
        rec = run_stream(spec)
        runs.append((rec, spec))
        (OUT / f"{spec['name']}.json").write_text(
            json.dumps(jsonable(rec), indent=2) + "\n", encoding="utf-8"
        )
    extras = []
    for title, fn, seed in (
        ("planted CATE on x1 in group 2", planted_dgp, 0),
        ("null: Y depends on X, not T", null_dgp, 1),
    ):
        X, Y, T, names = fn(n=500, seed=seed)
        rec = cfperm_groups(X, Y, T, names=names, n_perm=25, seed=seed)
        rec["localization"] = localize(X, Y, group_labels=T, n_perm=20, seed=seed)
        extras.append({"title": title, **compact(rec, {"name": title, "title": title, "t_meaning": "", "y_meaning": "", "path": None})})
        extras[-1]["top"] = rec["top"]
        extras[-1]["rejected"] = rec["rejected"]
        extras[-1]["hits"] = [n for n, h in zip(rec["names"], rec["hits"]) if h]
        (OUT / ("planted.json" if "planted" in title else "null.json")).write_text(
            json.dumps(jsonable(rec), indent=2) + "\n", encoding="utf-8"
        )
    summary = [compact(rec, spec) for rec, spec in runs]
    if runs:
        plot_importance(runs, OUT / "importance_by_stream.png")
    report = render_report(summary, extras)
    if skipped:
        report += "Skipped missing files:\n\n" + "\n".join(f"- `{p}`" for p in skipped) + "\n"
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(report)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
