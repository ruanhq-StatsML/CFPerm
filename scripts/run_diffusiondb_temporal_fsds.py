#!/usr/bin/env python3
"""DiffusionDB temporal attribution (prompt tokens → image attribute) via FSDS.

Definition (lightweight prototype — feature selection, not graph tip):
  X = [prompt embed (TF-IDF/hash)] ⊕ optional [cfg, step, sampler one-hot]
  Y = image_nsfw (continuous float in metadata; binarized for FSDS)
  T = time window (equal-count | equal-time | fixed calendar width)

Pipeline:
  1. Sample subset from ``metadata.parquet``
  2. Build prompt token X; optionally concatenate generation hyperparams
  3. Assign T by chosen window scheme; report signed ΔȲ
  4. RF Domain VIMP (X → T)
  5. FSDS: Scaler→Var→SelectKBest→HGB/LR on high-Y (early→late)
  6. Write ranking + report

Window schemes matter: equal-count vs equal-time vs width=Nd change
which tokens look like drift vs high-Y predictors.

  PYTHONPATH=. python3 scripts/run_diffusiondb_temporal_fsds.py \\
    --n-sample 4000 --concat-hyperparams --window-scheme equal_count \\
    --out results/diffusiondb_temporal_fsds

  PYTHONPATH=. python3 scripts/run_diffusiondb_window_sweep.py \\
    --n-sample 3000 --out results/diffusiondb_window_sweep
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_extraction.text import HashingVectorizer, TfidfVectorizer
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_META = ROOT / "data" / "diffusiondb" / "metadata.parquet"


def _tok_clean(s: str) -> str:
    s = (s or "").lower()
    s = re.sub(r"[^\w\s\-]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def load_subset(
    meta_path: Path,
    *,
    n_sample: int,
    seed: int,
) -> pd.DataFrame:
    cols = ["prompt", "timestamp", "image_nsfw", "prompt_nsfw", "cfg", "step", "sampler"]
    df = pd.read_parquet(meta_path, columns=cols)
    df = df.dropna(subset=["prompt", "timestamp", "image_nsfw"]).copy()
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True)
    df = df.sort_values("timestamp").reset_index(drop=True)
    if n_sample < len(df):
        q = pd.qcut(np.arange(len(df)), q=4, labels=False)
        rng = np.random.default_rng(seed)
        idx = []
        per = max(1, n_sample // 4)
        for g in range(4):
            pool = np.where(q == g)[0]
            take = min(per, len(pool))
            idx.append(rng.choice(pool, size=take, replace=False))
        idx = np.unique(np.concatenate(idx))
        if len(idx) < n_sample:
            rest = np.setdiff1d(np.arange(len(df)), idx)
            need = n_sample - len(idx)
            idx = np.concatenate(
                [idx, rng.choice(rest, size=min(need, len(rest)), replace=False)]
            )
        df = df.iloc[np.sort(idx)].reset_index(drop=True)
    df["prompt_clean"] = df["prompt"].map(_tok_clean)
    return df


def assign_time_windows(
    df: pd.DataFrame,
    *,
    n_windows: int = 2,
    scheme: str = "equal_count",
    width_hours: Optional[float] = None,
) -> pd.DataFrame:
    """Assign integer window id ``T`` (0 = earliest).

    Schemes
    -------
    equal_count : qcut on row rank — each window ≈ same n (calendar width varies).
    equal_time  : equal calendar spans — n per window may imbalance.
    width       : fixed ``width_hours`` bins from ts_min; early=min T, late=max T.
    """
    out = df.sort_values("timestamp").reset_index(drop=True).copy()
    ts = pd.to_datetime(out["timestamp"], utc=True)
    scheme = (scheme or "equal_count").lower()
    if scheme == "equal_count":
        ranks = np.arange(len(out))
        out["T"] = pd.qcut(ranks, q=n_windows, labels=False).astype(int)
        out.attrs["window_meta"] = {
            "scheme": "equal_count",
            "n_windows": int(n_windows),
            "note": "balanced n; calendar width varies by window",
        }
    elif scheme == "equal_time":
        t0, t1 = ts.min(), ts.max()
        edges = pd.date_range(t0, t1, periods=n_windows + 1)
        cats = pd.cut(ts, bins=edges, labels=False, include_lowest=True)
        out["T"] = cats.astype(int)
        out.attrs["window_meta"] = {
            "scheme": "equal_time",
            "n_windows": int(n_windows),
            "edges": [str(e) for e in edges],
            "note": "equal calendar span; n per window may imbalance",
        }
    elif scheme == "width":
        if width_hours is None or width_hours <= 0:
            raise ValueError("scheme=width requires width_hours > 0")
        t0 = ts.min()
        hours = (ts - t0).dt.total_seconds() / 3600.0
        out["T"] = np.floor(hours / float(width_hours)).astype(int)
        out.attrs["window_meta"] = {
            "scheme": "width",
            "width_hours": float(width_hours),
            "n_bins": int(out["T"].nunique()),
            "note": "fixed calendar width; early=min T, late=max T",
        }
    else:
        raise ValueError(f"unknown window scheme: {scheme}")
    return out


def build_light_token_X(
    prompts: List[str],
    *,
    method: str = "tfidf",
    max_features: int = 512,
    seed: int = 0,
) -> Tuple[np.ndarray, List[str], object]:
    """Lightweight prompt token embedding → dense X and feature names."""
    if method == "hash":
        vec = HashingVectorizer(
            n_features=max_features,
            alternate_sign=False,
            norm="l2",
            token_pattern=r"(?u)\b\w[\w\-]+\b",
            ngram_range=(1, 2),
        )
        X = vec.transform(prompts).astype(np.float64)
        names = [f"hash_{i}" for i in range(max_features)]
        return X.toarray(), names, vec
    vec = TfidfVectorizer(
        max_features=max_features,
        min_df=3,
        ngram_range=(1, 2),
        token_pattern=r"(?u)\b\w[\w\-]+\b",
        sublinear_tf=True,
    )
    X = vec.fit_transform(prompts).astype(np.float64).toarray()
    names = list(vec.get_feature_names_out())
    return X, names, vec


def concat_hyperparams(
    X_tok: np.ndarray,
    tok_names: List[str],
    df: pd.DataFrame,
) -> Tuple[np.ndarray, List[str]]:
    """Concatenate generation hyperparams onto prompt token X."""
    cfg = (
        pd.to_numeric(df.get("cfg"), errors="coerce")
        .fillna(0.0)
        .to_numpy(float)
        .reshape(-1, 1)
    )
    step = (
        pd.to_numeric(df.get("step"), errors="coerce")
        .fillna(0.0)
        .to_numpy(float)
        .reshape(-1, 1)
    )
    samp = df.get("sampler")
    if samp is None:
        samp = pd.Series(["unk"] * len(df))
    samp = samp.fillna("unk").astype(str)
    enc = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    samp_oh = enc.fit_transform(samp.to_numpy().reshape(-1, 1))
    samp_names = [f"sampler={c}" for c in enc.categories_[0]]

    def _z(a: np.ndarray) -> np.ndarray:
        mu, sd = float(a.mean()), float(a.std())
        return (a - mu) / (sd + 1e-8)

    Xp = np.hstack([X_tok, _z(cfg), _z(step), samp_oh])
    names = list(tok_names) + ["hp_cfg", "hp_step"] + samp_names
    return Xp, names


def rf_domain_vimp(
    X: np.ndarray, T: np.ndarray, *, seed: int, n_trees: int = 80
) -> np.ndarray:
    rf = RandomForestClassifier(
        n_estimators=n_trees,
        max_depth=10,
        min_samples_leaf=5,
        random_state=seed,
        n_jobs=1,
        class_weight="balanced_subsample",
    )
    rf.fit(X, T)
    return rf.feature_importances_.astype(float)


def run_fsds_y(
    X_tr: np.ndarray,
    y_tr: np.ndarray,
    X_te: np.ndarray,
    y_te: np.ndarray,
    feat_names: List[str],
    *,
    select_k: int,
    seed: int,
) -> Dict:
    """FSDS on binary high-Y: Scaler → Var → SelectKBest → HGB/LR."""
    if len(np.unique(y_tr)) < 2:
        return {"ok": False, "reason": "train needs both classes"}
    k = min(select_k, X_tr.shape[1], max(1, X_tr.shape[0] - 1))
    pipe = Pipeline(
        [
            ("sc", StandardScaler()),
            ("var", VarianceThreshold(1e-10)),
            ("sel", SelectKBest(f_classif, k=k)),
        ]
    )
    t0 = time.time()
    Xt = pipe.fit_transform(X_tr, y_tr)
    var_mask = pipe.named_steps["var"].get_support()
    cols_var = [c for c, m in zip(feat_names, var_mask) if m]
    sel = pipe.named_steps["sel"]
    selected = [c for c, m in zip(cols_var, sel.get_support()) if m]
    ranking = (
        pd.DataFrame({"feature": cols_var, "f_score": sel.scores_})
        .sort_values("f_score", ascending=False)
        .reset_index(drop=True)
    )
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["selected"] = ranking["feature"].isin(selected).astype(int)
    Xte = pipe.transform(X_te)
    out: Dict = {
        "ok": True,
        "n_selected": len(selected),
        "selected": selected,
        "ranking": ranking,
        "n_train": int(len(y_tr)),
        "n_test": int(len(y_te)),
        "pos_train": float(y_tr.mean()),
        "pos_test": float(y_te.mean()) if len(y_te) else float("nan"),
        "sec_fit": float(time.time() - t0),
        "models": {},
    }
    if len(y_te) == 0 or len(np.unique(y_te)) < 2:
        return out
    hgb = HistGradientBoostingClassifier(
        max_depth=4, learning_rate=0.08, max_iter=80, random_state=seed
    )
    hgb.fit(Xt, y_tr)
    ph = hgb.predict_proba(Xte)[:, 1]
    out["models"]["hgb"] = {
        "auc": float(roc_auc_score(y_te, ph)),
        "ap": float(average_precision_score(y_te, ph)),
    }
    lr = LogisticRegression(max_iter=400, random_state=seed)
    lr.fit(Xt, y_tr)
    pl = lr.predict_proba(Xte)[:, 1]
    out["models"]["logreg"] = {
        "auc": float(roc_auc_score(y_te, pl)),
        "ap": float(average_precision_score(y_te, pl)),
    }
    return out


def mean_shift_by_feature(
    X0: np.ndarray, X1: np.ndarray, names: List[str]
) -> pd.DataFrame:
    m0 = X0.mean(axis=0)
    m1 = X1.mean(axis=0)
    d = m1 - m0
    return (
        pd.DataFrame(
            {
                "feature": names,
                "mean_early": m0,
                "mean_late": m1,
                "delta": d,
                "abs_delta": np.abs(d),
                "sign": np.sign(d).astype(int),
            }
        )
        .sort_values("abs_delta", ascending=False)
        .reset_index(drop=True)
    )


def run_once(
    df: pd.DataFrame,
    *,
    embed: str,
    max_features: int,
    select_k: int,
    y_quantile: float,
    seed: int,
    concat_hyperparams_flag: bool,
    n_windows: int,
    window_scheme: str,
    width_hours: Optional[float],
) -> Dict:
    df_w = assign_time_windows(
        df, n_windows=n_windows, scheme=window_scheme, width_hours=width_hours
    )
    X, feat_names, _ = build_light_token_X(
        df_w["prompt_clean"].tolist(),
        method=embed,
        max_features=max_features,
        seed=seed,
    )
    d_tok = X.shape[1]
    if concat_hyperparams_flag:
        X, feat_names = concat_hyperparams(X, feat_names, df_w)
    T = df_w["T"].to_numpy(int)
    Y_cont = df_w["image_nsfw"].to_numpy(float)
    occupied = sorted(np.unique(T))
    if len(occupied) < 2:
        return {
            "ok": False,
            "reason": "need ≥2 occupied windows",
            "window_meta": df_w.attrs.get("window_meta"),
        }
    t_early, t_late = int(occupied[0]), int(occupied[-1])
    y_by_t = {int(t): float(Y_cont[T == t].mean()) for t in occupied}
    n_by_t = {int(t): int((T == t).sum()) for t in occupied}
    # calendar span per occupied window (shows scheme difference)
    span_h = {}
    for t in occupied:
        sub = df_w.loc[df_w["T"] == t, "timestamp"]
        span_h[int(t)] = float((sub.max() - sub.min()).total_seconds() / 3600.0)
    delta_y = y_by_t[t_late] - y_by_t[t_early]
    direction = (
        "positive/正向"
        if delta_y > 0
        else ("negative/负向" if delta_y < 0 else "flat")
    )
    thr_y = float(np.quantile(Y_cont[T == t_early], y_quantile))
    y_bin = (Y_cont >= thr_y).astype(int)
    vimp = rf_domain_vimp(X, T, seed=seed)
    cov = (
        pd.DataFrame({"feature": feat_names, "vimp_cov": vimp})
        .sort_values("vimp_cov", ascending=False)
        .reset_index(drop=True)
    )
    cov["rank"] = np.arange(1, len(cov) + 1)
    cmean = mean_shift_by_feature(X[T == t_early], X[T == t_late], feat_names)
    res = run_fsds_y(
        X[T == t_early],
        y_bin[T == t_early],
        X[T == t_late],
        y_bin[T == t_late],
        feat_names,
        select_k=select_k,
        seed=seed,
    )
    ranking = res.get("ranking")
    if ranking is None:
        return {"ok": False, "reason": res.get("reason", "fsds failed")}
    merged = ranking.merge(cov[["feature", "vimp_cov"]], on="feature", how="left")
    merged = merged.merge(
        cmean[["feature", "delta", "abs_delta", "sign", "mean_early", "mean_late"]],
        on="feature",
        how="left",
    )
    merged["score_blend"] = (
        merged["f_score"].rank(pct=True).fillna(0)
        + merged["vimp_cov"].rank(pct=True).fillna(0)
        + merged["abs_delta"].rank(pct=True).fillna(0)
    ) / 3.0
    merged = merged.sort_values("score_blend", ascending=False).reset_index(drop=True)
    merged["blend_rank"] = np.arange(1, len(merged) + 1)
    hp_mask = cov["feature"].astype(str).str.startswith(("hp_", "sampler="))
    hp_vimp_share = float(cov.loc[hp_mask, "vimp_cov"].sum()) if hp_mask.any() else 0.0
    return {
        "ok": True,
        "window_meta": dict(df_w.attrs.get("window_meta") or {}),
        "n": int(len(df_w)),
        "d_tokens": int(d_tok),
        "d_total": int(X.shape[1]),
        "concat_hyperparams": bool(concat_hyperparams_flag),
        "embed": embed,
        "t_early": t_early,
        "t_late": t_late,
        "n_by_T": n_by_t,
        "span_hours_by_T": span_h,
        "y_by_T": y_by_t,
        "delta_Y_image_nsfw": float(delta_y),
        "direction": direction,
        "y_threshold": thr_y,
        "hp_vimp_share": hp_vimp_share,
        "fsds": {
            "ok": res["ok"],
            "n_selected": res.get("n_selected"),
            "models": res.get("models"),
            "top_fsds": ranking.head(15)[["feature", "f_score"]].to_dict(
                orient="records"
            ),
        },
        "top_covariate": cov.head(15)[["feature", "vimp_cov"]].to_dict(
            orient="records"
        ),
        "top_cmean": cmean.head(15)[["feature", "delta", "sign"]].to_dict(
            orient="records"
        ),
        "top_blend": merged.head(20)[
            ["feature", "f_score", "vimp_cov", "delta", "sign", "score_blend"]
        ].to_dict(orient="records"),
        "_tables": {
            "cov": cov,
            "ranking": ranking,
            "cmean": cmean,
            "merged": merged,
            "df": df_w,
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", type=Path, default=DEFAULT_META)
    ap.add_argument("--n-sample", type=int, default=8000)
    ap.add_argument("--n-windows", type=int, default=2)
    ap.add_argument(
        "--window-scheme",
        choices=["equal_count", "equal_time", "width"],
        default="equal_count",
    )
    ap.add_argument("--width-hours", type=float, default=None)
    ap.add_argument("--max-features", type=int, default=512)
    ap.add_argument("--embed", choices=["tfidf", "hash"], default="tfidf")
    ap.add_argument(
        "--concat-hyperparams",
        action="store_true",
        help="concat cfg/step/sampler onto prompt TF-IDF/hash X",
    )
    ap.add_argument("--select-k", type=int, default=40)
    ap.add_argument("--y-quantile", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out", type=Path, default=ROOT / "results" / "diffusiondb_temporal_fsds"
    )
    args = ap.parse_args()

    if not args.meta.is_file():
        raise SystemExit(
            f"missing {args.meta}; download with:\n"
            f"  hf download poloclub/diffusiondb --type dataset "
            f"--include metadata.parquet --local-dir data/diffusiondb"
        )

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"load subset n={args.n_sample} from {args.meta}", flush=True)
    df = load_subset(args.meta, n_sample=args.n_sample, seed=args.seed)
    print(
        f"build X ({args.embed}"
        + (", ⊕ hyperparams" if args.concat_hyperparams else "")
        + f")  window={args.window_scheme}",
        flush=True,
    )
    out = run_once(
        df,
        embed=args.embed,
        max_features=args.max_features,
        select_k=args.select_k,
        y_quantile=args.y_quantile,
        seed=args.seed,
        concat_hyperparams_flag=args.concat_hyperparams,
        n_windows=args.n_windows,
        window_scheme=args.window_scheme,
        width_hours=args.width_hours,
    )
    if not out.get("ok"):
        raise SystemExit(out)
    tabs = out.pop("_tables")
    tabs["cov"].to_csv(args.out / "covariate_vimp_tokens.csv", index=False)
    tabs["ranking"].to_csv(args.out / "fsds_fscore_tokens.csv", index=False)
    tabs["cmean"].to_csv(args.out / "cmean_token_delta.csv", index=False)
    tabs["merged"].to_csv(args.out / "temporal_token_ranking.csv", index=False)
    tabs["df"][["prompt", "timestamp", "T", "image_nsfw", "prompt_nsfw"]].to_parquet(
        args.out / "subset_meta.parquet", index=False
    )
    summary = {
        **{k: v for k, v in out.items()},
        "ts_min": str(tabs["df"]["timestamp"].min()),
        "ts_max": str(tabs["df"]["timestamp"].max()),
        "sec": float(time.time() - t0),
    }
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n"
    )

    lines = [
        "# DiffusionDB temporal FS (prompt → image_nsfw)",
        "",
        "## Definition",
        f"- **X** = `{args.embed}` prompt tokens"
        + (" **⊕** cfg/step/sampler" if args.concat_hyperparams else "")
        + f" (d_tok={summary['d_tokens']}, d_total={summary['d_total']})",
        "- **Y** = continuous `image_nsfw` → early-quantile binary for FSDS",
        f"- **T** = `{args.window_scheme}` "
        f"(early T={summary['t_early']} → late T={summary['t_late']})",
        f"- window_meta: `{summary.get('window_meta')}`",
        "",
        f"Span: `{summary['ts_min']}` → `{summary['ts_max']}`.",
        "",
        "## Signed direction",
        f"- Ȳ_early={summary['y_by_T'][summary['t_early']]:.4f}, "
        f"Ȳ_late={summary['y_by_T'][summary['t_late']]:.4f}, "
        f"ΔȲ={summary['delta_Y_image_nsfw']:.4f} → **{summary['direction']}**",
        f"- n_by_T={summary['n_by_T']}",
        f"- span_hours_by_T={summary['span_hours_by_T']}",
        f"- hp_vimp_share={summary['hp_vimp_share']:.3f}",
        "",
        "## Methods (feature selection)",
        "1. RF Domain VIMP (X→T)",
        "2. cmean signed μ_late−μ_early",
        "3. FSDS SelectKBest→HGB/LR on high-Y; train early / test late",
        "",
        "Qwen embedding: optional later concat into X (deps not required here).",
        "",
    ]
    if summary["fsds"].get("models"):
        m = summary["fsds"]["models"]
        if "hgb" in m:
            lines.append(f"- HGB AUC(late)={m['hgb']['auc']:.3f}")
        if "logreg" in m:
            lines.append(f"- LogReg AUC(late)={m['logreg']['auc']:.3f}")
    lines += [
        "",
        "## Top blend",
        "",
        "| rank | feature | f | vimp | Δ |",
        "|---:|---|---:|---:|---:|",
    ]
    for i, r in enumerate(summary["top_blend"][:15], 1):
        lines.append(
            f"| {i} | `{r['feature']}` | {r['f_score']:.2f} | "
            f"{r['vimp_cov']:.4f} | {r['delta']:.4f} |"
        )
    (args.out / "DIFFUSIONDB_TEMPORAL_FSDS_REPORT.md").write_text(
        "\n".join(lines) + "\n"
    )
    print(
        json.dumps(
            {k: summary[k] for k in summary if not str(k).startswith("top_")},
            indent=2,
            default=str,
        )
    )
    print("top blend:", [r["feature"] for r in summary["top_blend"][:10]])
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
