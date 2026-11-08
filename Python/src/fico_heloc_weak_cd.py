"""FICO HELOC: weak multi-feature + highly nonlinear concept drift.

Justification design vs median-split covariate shift:
  - HELOC risk is driven by *several weak* features jointly, not one strong axis.
  - Randomly sample K=5 features; keep P(X) identical across batches (random split).
  - On the new batch only, redraw Y from a highly nonlinear function of those 5
    features (interactions / trig / tanh) → pure concept drift with weak drivers.
  - CFPerm power + driver recovery (recall@5 among top VIMP) vs the earlier
    single-feature median-split experiment.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from fico_heloc_cfperm import (
    VIMP_NOTIONS,
    _encode_y,
    cfperm_test,
    clean_heloc,
    load_heloc,
    power_summary,
)
@dataclass
class WeakCDDraw:
    seed: int
    driver_features: List[str]
    driver_idx: List[int]
    intensity: float
    n0: int
    n1: int
    mean_y0: float
    mean_y1: float
    delta_y: float


def _standardize(X: np.ndarray) -> np.ndarray:
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - mu) / sd


def nonlinear_concept_score(Z: np.ndarray) -> np.ndarray:
    """Highly nonlinear score on K standardized weak drivers (K>=5).

    Interactions + trig + tanh — each coordinate alone is weak; the joint
    map drives P(Y|X) under the alternative.
    """
    Z = np.asarray(Z, dtype=float)
    if Z.shape[1] < 5:
        raise ValueError("need at least 5 driver columns")
    z0, z1, z2, z3, z4 = (Z[:, j] for j in range(5))
    return (
        np.sin(np.pi * z0 * z1)
        + np.tanh(z2**2 * z3)
        - np.cos(z4) * z0
        + 0.5 * z1 * z2 * z3
        + 0.35 * (z4**3)
        + 0.25 * np.sin(z0 + z2 - z4) * z1
    )


def sample_weak_drivers(
    feature_names: Sequence[str], k: int = 5, seed: int = 0
) -> List[str]:
    rng = np.random.default_rng(seed)
    names = list(feature_names)
    if k > len(names):
        raise ValueError(f"k={k} > p={len(names)}")
    idx = rng.choice(len(names), size=k, replace=False)
    return [names[int(i)] for i in idx]


def make_weak_cd_batches(
    df: pd.DataFrame,
    *,
    y_col: str = "RiskPerformance",
    k_drivers: int = 5,
    intensity: float = 2.5,
    max_n_per_arm: Optional[int] = 1500,
    seed: int = 2026,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], WeakCDDraw]:
    """Random split + nonlinear CD on new arm only; P(X) shared."""
    feature_cols = [c for c in df.columns if c != y_col]
    drivers = sample_weak_drivers(feature_cols, k=k_drivers, seed=seed)
    driver_idx = [feature_cols.index(d) for d in drivers]

    rng = np.random.default_rng(seed + 7)
    X_all = df[feature_cols].to_numpy(dtype=float)
    Y_all = _encode_y(df[y_col])
    n = len(Y_all)
    perm = rng.permutation(n)
    half = n // 2
    idx0, idx1 = perm[:half], perm[half:]
    if max_n_per_arm is not None:
        if len(idx0) > max_n_per_arm:
            idx0 = rng.choice(idx0, size=max_n_per_arm, replace=False)
        if len(idx1) > max_n_per_arm:
            idx1 = rng.choice(idx1, size=max_n_per_arm, replace=False)

    X0, Y0 = X_all[idx0], Y_all[idx0].astype(float)
    X1 = X_all[idx1]
    # exist keeps observational labels; new redraws Y from nonlinear drivers
    Z1 = _standardize(X1[:, driver_idx])
    g = nonlinear_concept_score(Z1)
    # scale g to unit-ish then apply intensity
    g = (g - g.mean()) / (g.std() + 1e-8)
    p = 1.0 / (1.0 + np.exp(-intensity * g))
    Y1 = rng.binomial(1, p).astype(float)

    X = np.vstack([X0, X1])
    Y = np.concatenate([Y0, Y1])
    W = np.concatenate([np.zeros(len(Y0), dtype=int), np.ones(len(Y1), dtype=int)])
    meta = WeakCDDraw(
        seed=seed,
        driver_features=drivers,
        driver_idx=driver_idx,
        intensity=float(intensity),
        n0=int(len(Y0)),
        n1=int(len(Y1)),
        mean_y0=float(np.mean(Y0)),
        mean_y1=float(np.mean(Y1)),
        delta_y=float(np.mean(Y1) - np.mean(Y0)),
    )
    return X, Y, W, feature_cols, meta


def _recovery_metrics(
    top_features: Sequence[str],
    rejected: Sequence[str],
    drivers: Sequence[str],
    k: int = 5,
) -> Dict[str, float]:
    drivers_set = set(drivers)
    topk = list(top_features)[:k]
    rej = list(rejected)
    hit_top = [f for f in topk if f in drivers_set]
    hit_rej = [f for f in rej if f in drivers_set]
    return {
        "recall_top5": float(len(hit_top) / max(len(drivers_set), 1)),
        "precision_top5": float(len(hit_top) / max(len(topk), 1)),
        "n_drivers_in_top5": float(len(hit_top)),
        "n_drivers_rejected": float(len(hit_rej)),
        "any_driver_rejected": float(1.0 if hit_rej else 0.0),
        "any_driver_top5": float(1.0 if hit_top else 0.0),
    }


def run_weak_cd_experiment(
    *,
    cache_path: Optional[Path] = None,
    n_draws: int = 5,
    k_drivers: int = 5,
    intensity: float = 2.5,
    vimps: Sequence[str] = VIMP_NOTIONS,
    n_perm: int = 200,
    n_estimators: int = 150,
    n_perm_cate: int = 50,
    max_n_per_arm: Optional[int] = 1500,
    seed: int = 2026,
    n_jobs: int = 2,
    level_feature: float = 0.99,
    level_across_feature: float = 0.95,
    run_posthoc: bool = True,
) -> Dict[str, object]:
    raw = load_heloc(cache_path=cache_path)
    clean = clean_heloc(raw)
    cells = []
    draws_meta = []
    posthoc_summaries = []
    posthoc_feat = []
    posthoc_sub = []

    for d in range(n_draws):
        draw_seed = seed + 101 * d
        X, Y, W, names, meta = make_weak_cd_batches(
            clean,
            k_drivers=k_drivers,
            intensity=intensity,
            max_n_per_arm=max_n_per_arm,
            seed=draw_seed,
        )
        draws_meta.append(asdict(meta))
        print(
            f"[weak-cd] draw={d} drivers={meta.driver_features} "
            f"ΔȲ={meta.delta_y:.3f} n0={meta.n0} n1={meta.n1}",
            flush=True,
        )
        for vimp in vimps:
            print(f"  vimp={vimp} B={n_estimators} n_perm={n_perm}", flush=True)
            # recovered_split flag unused for multi-driver; use "" split label
            res = cfperm_test(
                X,
                Y,
                W,
                names,
                vimp=vimp,
                split_feature="",  # multi-driver; recovery via custom metrics
                n_perm=n_perm,
                n_estimators=n_estimators,
                n_perm_cate=n_perm_cate,
                level_feature=level_feature,
                level_across_feature=level_across_feature,
                seed=draw_seed + 17,
                n_jobs=n_jobs,
            )
            rec = _recovery_metrics(
                res.top_features, res.rejected_features, meta.driver_features, k=5
            )
            cell = asdict(res)
            cell.update(
                {
                    "draw": d,
                    "draw_seed": draw_seed,
                    "drivers": meta.driver_features,
                    "intensity": meta.intensity,
                    "delta_y": meta.delta_y,
                    **rec,
                }
            )
            # multi-driver "recovered" = any driver rejected
            cell["recovered_split"] = int(rec["any_driver_rejected"])
            cells.append(cell)

        if run_posthoc:
            X0, Y0 = X[W == 0], Y[W == 0]
            X1, Y1 = X[W == 1], Y[W == 1]
            # true_idx = all 5 drivers for hit@k
            true_idx = meta.driver_idx
            from posthoc_subset_localization import feature_localize, subset_on_driver

            feat_df, fmeta, _ = feature_localize(
                X0,
                Y0,
                X1,
                Y1,
                seed=draw_seed,
                names=names,
                true_idx=true_idx,
                top_k=5,
                n_perm=25,
            )
            sub_df, contrast, flag = subset_on_driver(
                X0,
                Y0,
                X1,
                Y1,
                B=30,
                j_star=int(fmeta["top_j"]),
                driver_idx=list(fmeta["driver_idx"]),
                seed=draw_seed + 3,
                n_perm=25,
            )
            posthoc_summaries.append(
                {
                    "draw": d,
                    "drivers": meta.driver_features,
                    "pfi_top": fmeta["top_name"],
                    "hit_at_k": fmeta["hit_at_k"],
                    "auc_abs": fmeta["auc_abs"],
                    "contrast": contrast,
                    "concentrated": flag,
                    "delta_y": meta.delta_y,
                }
            )
            posthoc_feat.append(feat_df.assign(draw=d))
            posthoc_sub.append(
                sub_df.assign(
                    draw=d,
                    drivers=",".join(meta.driver_features),
                    pfi_top=fmeta["top_name"],
                )
            )

    # power by vimp across draws
    by_vimp: Dict[str, Dict[str, float]] = {}
    for v in vimps:
        sub = [c for c in cells if c["vimp"] == v]
        if not sub:
            continue
        by_vimp[v] = {
            "n_cells": float(len(sub)),
            "reject_rate": float(np.mean([c["rejected"] for c in sub])),
            "any_driver_rejected_rate": float(
                np.mean([c["any_driver_rejected"] for c in sub])
            ),
            "any_driver_top5_rate": float(np.mean([c["any_driver_top5"] for c in sub])),
            "mean_recall_top5": float(np.mean([c["recall_top5"] for c in sub])),
            "mean_n_drivers_in_top5": float(
                np.mean([c["n_drivers_in_top5"] for c in sub])
            ),
        }
    overall = {
        "n_cells": float(len(cells)),
        "reject_rate": float(np.mean([c["rejected"] for c in cells])) if cells else 0.0,
        "any_driver_top5_rate": float(np.mean([c["any_driver_top5"] for c in cells]))
        if cells
        else 0.0,
        "mean_recall_top5": float(np.mean([c["recall_top5"] for c in cells]))
        if cells
        else 0.0,
    }
    return {
        "meta": {
            "n_raw": int(len(raw)),
            "n_clean": int(len(clean)),
            "n_draws": int(n_draws),
            "k_drivers": int(k_drivers),
            "intensity": float(intensity),
            "n_perm": int(n_perm),
            "n_estimators": int(n_estimators),
            "n_perm_cate": int(n_perm_cate),
            "max_n_per_arm": max_n_per_arm,
            "vimps": list(vimps),
            "seed": seed,
            "design": "weak5_nonlinear_concept_drift",
        },
        "draws": draws_meta,
        "results": cells,
        "power": {"overall": overall, "by_vimp": by_vimp},
        "posthoc_summary": posthoc_summaries,
        "posthoc_feature_rank": pd.concat(posthoc_feat, ignore_index=True)
        if posthoc_feat
        else pd.DataFrame(),
        "posthoc_subset_driver": pd.concat(posthoc_sub, ignore_index=True)
        if posthoc_sub
        else pd.DataFrame(),
    }


def load_median_split_power(path: Path) -> Optional[pd.DataFrame]:
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_csv(p)


def write_weak_cd_outputs(
    payload: Dict[str, object],
    out_dir: Path,
    *,
    baseline_power_csv: Optional[Path] = None,
) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "weak_cd_results.json", "w", encoding="utf-8") as f:
        # drop big frames from json
        slim = {
            k: v
            for k, v in payload.items()
            if k not in ("posthoc_feature_rank", "posthoc_subset_driver")
        }
        json.dump(slim, f, indent=2, default=str)

    pd.DataFrame(payload["results"]).to_csv(out_dir / "weak_cd_cells.csv", index=False)
    pd.DataFrame(payload["draws"]).to_csv(out_dir / "weak_cd_draws.csv", index=False)
    power = payload["power"]
    rows = [{"scope": "overall", **power["overall"]}]
    for v, d in power["by_vimp"].items():
        rows.append({"scope": v, **d})
    power_df = pd.DataFrame(rows)
    power_df.to_csv(out_dir / "weak_cd_power.csv", index=False)

    if payload.get("posthoc_summary"):
        pd.DataFrame(payload["posthoc_summary"]).to_csv(
            out_dir / "weak_cd_posthoc_summary.csv", index=False
        )
    fr = payload.get("posthoc_feature_rank")
    if isinstance(fr, pd.DataFrame) and len(fr):
        fr.to_csv(out_dir / "weak_cd_posthoc_feature_rank.csv", index=False)
    sd = payload.get("posthoc_subset_driver")
    if isinstance(sd, pd.DataFrame) and len(sd):
        sd.to_csv(out_dir / "weak_cd_posthoc_subset_driver.csv", index=False)

    # comparison vs median-split baseline
    base = load_median_split_power(baseline_power_csv) if baseline_power_csv else None
    cmp_rows = []
    if base is not None:
        for _, r in base.iterrows():
            cmp_rows.append(
                {
                    "design": "median_split_covshift",
                    "scope": r["scope"],
                    "reject_rate": r.get("reject_rate", np.nan),
                    "recover_rate": r.get("recover_rate", np.nan),
                    "top1_hit_rate": r.get("top1_hit_rate", np.nan),
                    "top3_hit_rate": r.get("top3_hit_rate", np.nan),
                    "mean_recall_top5": np.nan,
                    "any_driver_top5_rate": np.nan,
                }
            )
    for r in rows:
        cmp_rows.append(
            {
                "design": "weak5_nonlinear_cd",
                "scope": r["scope"],
                "reject_rate": r.get("reject_rate", np.nan),
                "recover_rate": r.get("any_driver_rejected_rate", np.nan),
                "top1_hit_rate": np.nan,
                "top3_hit_rate": np.nan,
                "mean_recall_top5": r.get("mean_recall_top5", np.nan),
                "any_driver_top5_rate": r.get(
                    "any_driver_top5_rate", r.get("any_driver_top5_rate", np.nan)
                ),
            }
        )
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(out_dir / "power_comparison.csv", index=False)

    meta = payload["meta"]
    lines = [
        "# FICO HELOC — weak-5 + highly nonlinear concept drift",
        "",
        "Design justification: HELOC is driven by **several weak features**, not one",
        "strong median axis. We randomly sample `K=5` features, keep `P(X)` identical",
        "across batches, and redraw `Y` on the new batch from a highly nonlinear map",
        "of those five (sin/tanh/interactions/cubic) → pure concept drift.",
        "",
        f"- Clean n=`{meta['n_clean']}`, draws=`{meta['n_draws']}`, "
        f"intensity=`{meta['intensity']}`, `B={meta['n_estimators']}`, "
        f"`n_perm={meta['n_perm']}`, `max_n_per_arm={meta['max_n_per_arm']}`.",
        "",
        "## Power (this design)",
        "",
        "| scope | reject | any_driver@top5 | mean recall@5 | n |",
        "|---|---:|---:|---:|---:|",
    ]
    ov = power["overall"]
    lines.append(
        f"| overall | {ov['reject_rate']:.3f} | {ov['any_driver_top5_rate']:.3f} | "
        f"{ov['mean_recall_top5']:.3f} | {int(ov['n_cells'])} |"
    )
    for v, d in power["by_vimp"].items():
        lines.append(
            f"| {v} | {d['reject_rate']:.3f} | {d['any_driver_top5_rate']:.3f} | "
            f"{d['mean_recall_top5']:.3f} | {int(d['n_cells'])} |"
        )

    lines += ["", "## vs median-split covariate shift (supplement)", "",]
    if base is not None:
        lines += [
            "| design | scope | reject | recover / driver-hit |",
            "|---|---|---:|---:|",
        ]
        for _, r in cmp_df.iterrows():
            hit = r["recover_rate"]
            if r["design"] == "weak5_nonlinear_cd":
                hit = r["any_driver_top5_rate"]
            lines.append(
                f"| {r['design']} | {r['scope']} | {r['reject_rate']:.3f} | "
                f"{hit if hit == hit else float('nan'):.3f} |"
            )
    lines += [
        "",
        "### Takeaway",
        "",
        "- Median-split on one LR feature = **single strong covariate-shift** axis;",
        "  power concentrates on GRF recovery of that axis.",
        "- Weak-5 nonlinear CD = **several weak drivers** of `P(Y|X)`; rejection +",
        "  recall@5 justify CFPerm under the multi-weak-feature regime that matches HELOC.",
        "",
        "Files: `weak_cd_cells.csv`, `weak_cd_power.csv`, `power_comparison.csv`,",
        "`weak_cd_posthoc_*.csv`.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
