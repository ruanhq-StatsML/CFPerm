"""Post-hoc subset localization after CFPerm / OOD-VIMP feature selection.

Matches the notebook / GitHub localization script:
  1. feature_localize — PFI(new) − PFI(exist) ranks shift drivers
  2. subset_on_driver — median-split on j★; bootstrap ΔVIMP + ΔMSE on low/high
  3. subset_localize  — cut on z-cols; PSI / TVD / Hellinger + SMD on selected x
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score

N_PERM = 25
TOP_K = 5
N_MIN = 20
N_CAP = 1500
Z95 = 1.96
B_BOOT = 30


def _dist(y0, y1, nbins: int = 20) -> Dict[str, float]:
    y0, y1 = np.asarray(y0, float), np.asarray(y1, float)
    y0, y1 = y0[~np.isnan(y0)], y1[~np.isnan(y1)]
    if len(y0) < 2 or len(y1) < 2:
        return dict(PSI=np.nan, TVD=np.nan, Hellinger=np.nan)
    edges = np.unique(np.quantile(np.r_[y0, y1], np.linspace(0, 1, nbins + 1)))
    if len(edges) < 2:
        return dict(PSI=np.nan, TVD=np.nan, Hellinger=np.nan)
    p, _ = np.histogram(y0, bins=edges)
    q, _ = np.histogram(y1, bins=edges)
    p = np.clip(p / max(p.sum(), 1), 1e-12, None)
    q = np.clip(q / max(q.sum(), 1), 1e-12, None)
    return dict(
        PSI=float(((p - q) * np.log(p / q)).sum()),
        TVD=float(0.5 * np.abs(p - q).sum()),
        Hellinger=float(np.sqrt(0.5 * ((np.sqrt(p) - np.sqrt(q)) ** 2).sum())),
    )


def _cut_ids(z0, z1, nbins: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """Fit cuts on z0 only; return group ids for both batches."""
    s = pd.Series(z0)
    if s.nunique(dropna=True) <= 6:
        levels = {v: i for i, v in enumerate(pd.unique(s.dropna()))}
        g0 = np.array([levels.get(v, -1) if pd.notna(v) else -1 for v in z0], int)
        g1 = np.array([levels.get(v, -1) if pd.notna(v) else -1 for v in z1], int)
        return g0, g1
    edges = np.unique(
        np.quantile(
            pd.to_numeric(s, errors="coerce").dropna(),
            np.linspace(0, 1, nbins + 1),
        )
    )
    if len(edges) < 2:
        return np.full(len(z0), -1, int), np.full(len(z1), -1, int)

    def dig(z):
        v = pd.to_numeric(pd.Series(z), errors="coerce").to_numpy(float)
        return np.where(np.isnan(v), -1, np.digitize(v, edges[1:-1])).astype(int)

    return dig(z0), dig(z1)


def pfi(X, Y, seed: int, n_perm: int = N_PERM) -> np.ndarray:
    Xs = np.asarray(X, float)
    Ys = np.asarray(Y, float).ravel()
    if len(Ys) < 16:
        return np.zeros(Xs.shape[1])
    m = Ridge(alpha=1.0, random_state=int(seed) % (2**31 - 1)).fit(Xs, Ys)
    base = float(np.mean((Ys - m.predict(Xs)) ** 2))
    out = np.zeros(Xs.shape[1])
    for j in range(Xs.shape[1]):
        rng = np.random.default_rng(seed + 1000 + j)
        stack = np.zeros(n_perm)
        for i in range(n_perm):
            Xp = Xs.copy()
            Xp[:, j] = Xp[rng.permutation(len(Xp)), j]
            stack[i] = float(np.mean((Ys - m.predict(Xp)) ** 2)) - base
        out[j] = float(stack.mean())
    return out


def mse_pair(X0, Y0, X1, Y1, seed: int) -> Tuple[float, float, float]:
    m = Ridge(alpha=1.0, random_state=seed).fit(X0, Y0)
    mse0 = float(np.mean((Y0 - m.predict(X0)) ** 2))
    mse1 = float(np.mean((Y1 - m.predict(X1)) ** 2))
    return mse1 - mse0, mse0, mse1


def feature_localize(
    X0,
    Y0,
    X1,
    Y1,
    seed: int,
    names: Sequence[str],
    true_idx: Optional[Sequence[int]] = None,
    top_k: int = TOP_K,
    n_perm: int = N_PERM,
) -> Tuple[pd.DataFrame, Dict, np.ndarray]:
    p0 = pfi(X0, Y0, seed + 11, n_perm=n_perm)
    p1 = pfi(X1, Y1, seed + 13, n_perm=n_perm)
    degree = p1 - p0
    order = np.argsort(-np.abs(degree))
    top = order[: min(top_k, len(order))]
    hit = np.nan
    auc = np.nan
    if true_idx is not None:
        true = set(int(i) for i in true_idx)
        hit = float(sum(int(j in true) for j in top)) / max(len(top), 1)
        y = np.array([1 if j in true else 0 for j in range(len(degree))])
        if y.min() != y.max():
            auc = float(roc_auc_score(y, np.abs(degree)))
    rows = []
    for r, j in enumerate(order):
        rows.append(
            dict(
                rank=r + 1,
                feature=str(names[j]) if j < len(names) else f"X{j}",
                j=int(j),
                degree=float(degree[j]),
                pfi_exist=float(p0[j]),
                pfi_new=float(p1[j]),
                is_true=int(
                    true_idx is not None and int(j) in set(int(i) for i in true_idx)
                ),
            )
        )
    meta = dict(
        top_j=int(top[0]),
        top_name=str(names[top[0]]) if top[0] < len(names) else f"X{top[0]}",
        top_degree=float(degree[top[0]]),
        energy=float(np.max(np.abs(degree))),
        hit_at_k=hit,
        auc_abs=auc,
        driver_idx=top.tolist(),
    )
    return pd.DataFrame(rows), meta, degree


def subset_on_driver(
    X0,
    Y0,
    X1,
    Y1,
    B: int,
    j_star: int,
    driver_idx: Sequence[int],
    seed: int,
    n_perm: int = N_PERM,
    n_min: int = N_MIN,
) -> Tuple[pd.DataFrame, float, str]:
    """Median split on driver feature → low/high; ΔVIMP and ΔMSE on each."""
    col = np.concatenate([X0[:, j_star], X1[:, j_star]])
    thr = float(np.median(col))
    rows = []
    driver_idx = list(driver_idx)
    for label, op in (("low", np.less), ("high", np.greater_equal)):
        m0 = op(X0[:, j_star], thr)
        m1 = op(X1[:, j_star], thr)
        n0, n1 = int(m0.sum()), int(m1.sum())
        if n0 < n_min or n1 < n_min:
            rows.append(
                dict(
                    subset=label,
                    n0=n0,
                    n1=n1,
                    thr=thr,
                    dvimp=np.nan,
                    dvimp_lo=np.nan,
                    dvimp_hi=np.nan,
                    dmse=np.nan,
                    skipped="too_small",
                )
            )
            continue
        boots = []
        for b in range(B):
            p0 = pfi(X0[m0], Y0[m0], seed + 200 + b, n_perm=n_perm)
            p1 = pfi(X1[m1], Y1[m1], seed + 400 + b, n_perm=n_perm)
            boots.append(float(np.mean((p1 - p0)[driver_idx])))
        boots = np.asarray(boots, float)
        dv = float(boots.mean())
        se = float(boots.std(ddof=1) / np.sqrt(max(len(boots), 1)))
        dmse, _, _ = mse_pair(X0[m0], Y0[m0], X1[m1], Y1[m1], seed + 7)
        rows.append(
            dict(
                subset=label,
                n0=n0,
                n1=n1,
                thr=thr,
                dvimp=dv,
                dvimp_lo=dv - Z95 * se,
                dvimp_hi=dv + Z95 * se,
                dmse=float(dmse),
                skipped="",
            )
        )
    sdf = pd.DataFrame(rows)
    if sdf.dvimp.notna().sum() == 2:
        hi = float(sdf.loc[sdf.subset == "high", "dvimp"].iloc[0])
        lo = float(sdf.loc[sdf.subset == "low", "dvimp"].iloc[0])
        contrast = abs(hi) - abs(lo)
        flag = "high" if abs(hi) >= abs(lo) else "low"
    else:
        contrast, flag = np.nan, ""
    return sdf, float(contrast) if contrast == contrast else np.nan, flag


def subset_localize(
    df0: pd.DataFrame,
    df1: pd.DataFrame,
    *,
    x_cols: Sequence[str],
    y_col: str,
    z_cols: Sequence[str],
    selected: Optional[Sequence[str]] = None,
    top_k: int = TOP_K,
    nbins: int = 5,
    min_n: int = N_MIN,
) -> Dict:
    """Cut on z_cols; report SMD + PSI/TVD/Hellinger for selected x (or top |Δmean|+|Δstd|)."""
    X0 = df0[list(x_cols)].to_numpy(float)
    X1 = df1[list(x_cols)].to_numpy(float)
    if selected is None:
        scores = np.abs(X0.mean(0) - X1.mean(0)) + np.abs(X0.std(0) - X1.std(0))
        selected = [x_cols[j] for j in np.argsort(-scores)[: min(top_k, len(x_cols))]]
    selected = list(selected)

    rows = []
    for zc in z_cols:
        g0, g1 = _cut_ids(df0[zc].to_numpy(), df1[zc].to_numpy(), nbins=nbins)
        for xc in selected:
            y0 = df0[xc].to_numpy(float)
            y1 = df1[xc].to_numpy(float)
            for grp in set(g0) | set(g1):
                if grp < 0:
                    continue
                a, b = y0[g0 == grp], y1[g1 == grp]
                a, b = a[~np.isnan(a)], b[~np.isnan(b)]
                if len(a) < min_n or len(b) < min_n:
                    continue
                rows.append(
                    dict(
                        z_cut=zc,
                        y_from_x=xc,
                        grp=int(grp),
                        mean0=float(a.mean()),
                        mean1=float(b.mean()),
                        SMD=float(a.mean() - b.mean()),
                        n0=len(a),
                        n1=len(b),
                        **_dist(a, b),
                    )
                )
    return {"selected": selected, "table": pd.DataFrame(rows)}


def run_posthoc_on_batches(
    X0: np.ndarray,
    Y0: np.ndarray,
    X1: np.ndarray,
    Y1: np.ndarray,
    names: Sequence[str],
    *,
    split_feature: str,
    seed: int = 2026,
    B: int = B_BOOT,
    top_k: int = TOP_K,
    n_perm: int = N_PERM,
    z_cols: Optional[Sequence[str]] = None,
    true_feature: Optional[str] = None,
) -> Dict:
    """Full post-hoc pipeline for one exist/new pair."""
    names = list(names)
    true_idx = None
    if true_feature is not None and true_feature in names:
        true_idx = [names.index(true_feature)]

    feat_df, meta, _ = feature_localize(
        X0, Y0, X1, Y1, seed, names, true_idx=true_idx, top_k=top_k, n_perm=n_perm
    )
    j_star = int(meta["top_j"])
    driver_idx = list(meta["driver_idx"])
    sub_df, contrast, flag = subset_on_driver(
        X0, Y0, X1, Y1, B, j_star, driver_idx, seed, n_perm=n_perm
    )
    sub_df = sub_df.assign(
        split_feature=split_feature,
        driver=meta["top_name"],
        contrast=contrast,
        concentrated=flag,
    )

    # z-cut table on top drivers (+ split feature)
    df0 = pd.DataFrame(X0, columns=names)
    df1 = pd.DataFrame(X1, columns=names)
    df0["_Y"] = Y0
    df1["_Y"] = Y1
    selected = [meta["top_name"]] + [
        r.feature for _, r in feat_df.head(top_k).iterrows() if r.feature != meta["top_name"]
    ]
    selected = selected[:top_k]
    if split_feature in names and split_feature not in selected:
        selected = [split_feature] + selected
        selected = selected[: top_k + 1]
    if z_cols is None:
        z_cols = [split_feature, meta["top_name"]]
        z_cols = list(dict.fromkeys([z for z in z_cols if z in names]))
    loc = subset_localize(
        df0,
        df1,
        x_cols=names,
        y_col="_Y",
        z_cols=z_cols,
        selected=selected,
    )
    if not loc["table"].empty:
        loc["table"] = loc["table"].assign(split_feature=split_feature)

    summary = dict(
        split_feature=split_feature,
        driver=meta["top_name"],
        top_degree=meta["top_degree"],
        energy=meta["energy"],
        hit_at_k=meta["hit_at_k"],
        auc_abs=meta["auc_abs"],
        contrast=contrast,
        concentrated=flag,
        dmse_low=float(sub_df.loc[sub_df.subset == "low", "dmse"].iloc[0])
        if (sub_df.subset == "low").any()
        else np.nan,
        dmse_high=float(sub_df.loc[sub_df.subset == "high", "dmse"].iloc[0])
        if (sub_df.subset == "high").any()
        else np.nan,
    )
    return {
        "feature_rank": feat_df.assign(split_feature=split_feature),
        "subset_driver": sub_df,
        "zcut": loc["table"],
        "summary": summary,
        "meta": meta,
    }
