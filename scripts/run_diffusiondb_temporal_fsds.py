#!/usr/bin/env python3
"""DiffusionDB temporal attribution (prompt tokens → image attribute) via FSDS.

Definition (lightweight prototype):
  X = prompt token features (Tfidf / hashed bag-of-tokens — light embedding)
  Y = image_nsfw (image-side attribute from metadata; no PNG download)
  T = time window from timestamp (early vs late within the 2M span)

Pipeline:
  1. Sample a small subset from ``metadata.parquet``
  2. Build lightweight token embedding X
  3. Assign T windows; report signed ΔȲ
  4. Covariate drift: RF Domain VIMP (X → T)
  5. Concept / outcome: FSDS (Scaler→Var→SelectKBest→HGB/LR) predicting
     high-Y from X on early window; score tokens on late holdout
  6. Write ranking + report

  PYTHONPATH=. python3 scripts/run_diffusiondb_temporal_fsds.py \\
    --n-sample 8000 --out results/diffusiondb_temporal_fsds
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
from sklearn.preprocessing import StandardScaler


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
        # stratified by time quartile indices for coverage
        q = pd.qcut(np.arange(len(df)), q=4, labels=False)
        rng = np.random.default_rng(seed)
        idx = []
        per = max(1, n_sample // 4)
        for g in range(4):
            pool = np.where(q == g)[0]
            take = min(per, len(pool))
            idx.append(rng.choice(pool, size=take, replace=False))
        idx = np.unique(np.concatenate(idx))
        # top-up if short
        if len(idx) < n_sample:
            rest = np.setdiff1d(np.arange(len(df)), idx)
            need = n_sample - len(idx)
            idx = np.concatenate([idx, rng.choice(rest, size=min(need, len(rest)), replace=False)])
        df = df.iloc[np.sort(idx)].reset_index(drop=True)
    df["prompt_clean"] = df["prompt"].map(_tok_clean)
    return df


def assign_time_windows(df: pd.DataFrame, *, n_windows: int = 2) -> pd.DataFrame:
    """Equal-count time windows along sorted timestamps (T=0 early …)."""
    out = df.sort_values("timestamp").reset_index(drop=True).copy()
    ranks = np.arange(len(out))
    # qcut on ranks for balanced windows
    out["T"] = pd.qcut(ranks, q=n_windows, labels=False).astype(int)
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
        # hashing has no vocabulary — synthetic names
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


def rf_domain_vimp(X: np.ndarray, T: np.ndarray, *, seed: int, n_trees: int = 80) -> np.ndarray:
    """Covariate-drift importance: RF predicting time window from X."""
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
    out: Dict = {
        "ok": True,
        "n_train": int(len(y_tr)),
        "n_test": int(len(y_te)),
        "pos_train": float(y_tr.mean()),
        "pos_test": float(y_te.mean()) if len(y_te) else float("nan"),
        "n_selected": len(selected),
        "selected": selected,
        "ranking": ranking,
        "sec_select": float(time.time() - t0),
        "models": {},
    }
    if len(y_te) == 0 or len(np.unique(y_te)) < 2:
        out["models"]["note"] = "test missing both classes"
        return out
    Xv = pipe.transform(X_te)
    hgb = HistGradientBoostingClassifier(
        max_depth=6, learning_rate=0.08, max_iter=100, random_state=seed
    )
    hgb.fit(Xt, y_tr)
    ph = hgb.predict_proba(Xv)[:, 1]
    out["models"]["hgb"] = {
        "auc": float(roc_auc_score(y_te, ph)),
        "ap": float(average_precision_score(y_te, ph)),
    }
    lr = LogisticRegression(max_iter=400, C=0.5, class_weight="balanced", random_state=seed)
    lr.fit(Xt, y_tr)
    pl = lr.predict_proba(Xv)[:, 1]
    out["models"]["logreg"] = {
        "auc": float(roc_auc_score(y_te, pl)),
        "ap": float(average_precision_score(y_te, pl)),
    }
    return out


def mean_shift_by_feature(X0: np.ndarray, X1: np.ndarray, names: List[str]) -> pd.DataFrame:
    """Signed tip-cmean style: μ_late − μ_early per token feature."""
    m0, m1 = X0.mean(axis=0), X1.mean(axis=0)
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", type=Path, default=DEFAULT_META)
    ap.add_argument("--n-sample", type=int, default=8000)
    ap.add_argument("--n-windows", type=int, default=2, help="2 = early vs late")
    ap.add_argument("--max-features", type=int, default=512)
    ap.add_argument("--embed", choices=["tfidf", "hash"], default="tfidf")
    ap.add_argument("--select-k", type=int, default=40)
    ap.add_argument("--y-quantile", type=float, default=0.7, help="high-Y threshold on early window")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "diffusiondb_temporal_fsds")
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
    df = assign_time_windows(df, n_windows=args.n_windows)
    print(
        f"  rows={len(df)}  T spans: "
        + ", ".join(
            f"T={t} [{df.loc[df['T']==t,'timestamp'].min()} .. {df.loc[df['T']==t,'timestamp'].max()}] n={int((df['T']==t).sum())}"
            for t in sorted(df["T"].unique())
        ),
        flush=True,
    )

    print(f"build light token X ({args.embed}, d≤{args.max_features})", flush=True)
    X, feat_names, _vec = build_light_token_X(
        df["prompt_clean"].tolist(),
        method=args.embed,
        max_features=args.max_features,
        seed=args.seed,
    )
    T = df["T"].to_numpy(int)
    Y_cont = df["image_nsfw"].to_numpy(float)

    # signed outcome direction early → late
    y_by_t = {int(t): float(Y_cont[T == t].mean()) for t in np.unique(T)}
    t_early, t_late = int(T.min()), int(T.max())
    delta_y = y_by_t[t_late] - y_by_t[t_early]
    direction = "positive/正向" if delta_y > 0 else ("negative/负向" if delta_y < 0 else "flat")

    # binary Y from early-window quantile (avoid late leakage into threshold)
    thr_y = float(np.quantile(Y_cont[T == t_early], args.y_quantile))
    y_bin = (Y_cont >= thr_y).astype(int)

    # --- covariate drift: X → T ---
    print("RF domain VIMP (covariate drift X→T)", flush=True)
    vimp = rf_domain_vimp(X, T, seed=args.seed)
    cov = (
        pd.DataFrame({"feature": feat_names, "vimp_cov": vimp})
        .sort_values("vimp_cov", ascending=False)
        .reset_index(drop=True)
    )
    cov["rank"] = np.arange(1, len(cov) + 1)

    # --- tip-cmean signed token shift ---
    cmean = mean_shift_by_feature(X[T == t_early], X[T == t_late], feat_names)

    # --- FSDS: train on early, holdout late ---
    print("FSDS train early → test late (predict high image_nsfw)", flush=True)
    res = run_fsds_y(
        X[T == t_early],
        y_bin[T == t_early],
        X[T == t_late],
        y_bin[T == t_late],
        feat_names,
        select_k=args.select_k,
        seed=args.seed,
    )
    ranking = res.get("ranking")
    if ranking is None:
        raise SystemExit(f"FSDS failed: {res}")

    # merge views
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

    cov.to_csv(args.out / "covariate_vimp_tokens.csv", index=False)
    ranking.to_csv(args.out / "fsds_fscore_tokens.csv", index=False)
    cmean.to_csv(args.out / "cmean_token_delta.csv", index=False)
    merged.to_csv(args.out / "temporal_token_ranking.csv", index=False)
    df[["prompt", "timestamp", "T", "image_nsfw", "prompt_nsfw"]].to_parquet(
        args.out / "subset_meta.parquet", index=False
    )

    summary = {
        "n": int(len(df)),
        "d_tokens": len(feat_names),
        "embed": args.embed,
        "n_windows": args.n_windows,
        "t_early": t_early,
        "t_late": t_late,
        "ts_min": str(df["timestamp"].min()),
        "ts_max": str(df["timestamp"].max()),
        "y_by_T": y_by_t,
        "delta_Y_image_nsfw": float(delta_y),
        "direction": direction,
        "y_threshold": thr_y,
        "fsds": {
            "ok": res["ok"],
            "n_selected": res.get("n_selected"),
            "models": res.get("models"),
            "top_fsds": ranking.head(15)[["feature", "f_score"]].to_dict(orient="records"),
        },
        "top_covariate": cov.head(15)[["feature", "vimp_cov"]].to_dict(orient="records"),
        "top_cmean": cmean.head(15)[["feature", "delta", "sign"]].to_dict(orient="records"),
        "top_blend": merged.head(20)[
            ["feature", "f_score", "vimp_cov", "delta", "sign", "score_blend"]
        ].to_dict(orient="records"),
        "sec": float(time.time() - t0),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))

    lines = [
        "# DiffusionDB temporal attribution (prompt tokens → image_nsfw)",
        "",
        "## Definition",
        "- **X** = lightweight prompt token embedding "
        f"(`{args.embed}`, d={len(feat_names)})",
        "- **Y** = `image_nsfw` (image attribute from metadata; no PNG download)",
        f"- **T** = {args.n_windows} equal-count time windows "
        f"(early T={t_early} → late T={t_late})",
        "",
        f"Span: `{summary['ts_min']}` → `{summary['ts_max']}` (2M gallery; short calendar span).",
        "",
        "## Signed direction (outcome)",
        f"- Ȳ_early={y_by_t[t_early]:.4f}, Ȳ_late={y_by_t[t_late]:.4f}, "
        f"ΔȲ={delta_y:.4f} → **{direction}**",
        "",
        "## Methods",
        "1. **Covariate drift**: RF Domain VIMP (X→T)",
        "2. **Tip-cmean**: signed μ_late−μ_early per token feature",
        "3. **FSDS**: Scaler→Var→SelectKBest→HGB/LR predicting high-Y "
        "(threshold from early quantile); train early / test late",
        "4. **Blend rank**: mean of percentile ranks (f_score, vimp_cov, |Δ|)",
        "",
    ]
    if res.get("models"):
        m = res["models"]
        if "hgb" in m:
            lines.append(f"- FSDS HGB AUC(late)={m['hgb']['auc']:.3f} AP={m['hgb']['ap']:.3f}")
        if "logreg" in m:
            lines.append(
                f"- FSDS LogReg AUC(late)={m['logreg']['auc']:.3f} AP={m['logreg']['ap']:.3f}"
            )
    lines += ["", "## Top blend tokens", "", "| rank | token | f_score | vimp_cov | Δ | sign |", "|---:|---|---:|---:|---:|---:|"]
    for i, r in enumerate(summary["top_blend"][:15], 1):
        lines.append(
            f"| {i} | `{r['feature']}` | {r['f_score']:.2f} | {r['vimp_cov']:.4f} | "
            f"{r['delta']:.4f} | {int(r['sign'])} |"
        )
    lines += [
        "",
        "## Note",
        "Full CLIP token grid (77×768) can replace the light embedding later; "
        "protocol (T windows, FSDS, cmean sign, domain VIMP) stays the same.",
        "",
        f"Artifacts under `{args.out}/`.",
    ]
    (args.out / "DIFFUSIONDB_TEMPORAL_FSDS_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k not in {"top_blend", "top_fsds", "top_covariate", "top_cmean"}}, indent=2))
    print("top blend:", [r["feature"] for r in summary["top_blend"][:10]])
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
