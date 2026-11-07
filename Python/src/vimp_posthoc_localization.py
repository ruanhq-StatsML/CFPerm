#!/usr/bin/env python3
"""Post-hoc VIMP localization after the impossibility issue.

Impossibility: ΔVIMP cannot uniquely decompose CD vs CS vs ε.
Localization (this script): do not decompose type — find *which features*
moved, then *which observation subset* has the large two-batch VIMP gap.

Procedure (PDF: localize drivers, then compare cohorts)
-------------------------------------------------------
1. Feature rank: degree_j = PFI_new(j) − PFI_exist(j). Top-k = drivers.
2. Split observations on the top driver (median → two subsets: low / high).
3. On each subset, recompute two-batch ΔVIMP (mean degree on the driver block)
   and ΔMSE. The subset with the larger |ΔVIMP| is the localized cohort.

Repeats: injected CD/CS DGPs + VIMP preprocess datasets
(air_quality, bank, diabetes, daily_slice, m5).

  PYTHONPATH=src:scripts python scripts/vimp_posthoc_localization.py
  → results/method_vimp_posthoc_localization/
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path
from warnings import filterwarnings

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "scripts")]

filterwarnings("ignore")

OUT = ROOT / "results" / "method_vimp_posthoc_localization"
SEED = 1250
N_PERM = 16
N_CAP = 1200
TOP_K = 4
N_MIN = 80
Z95 = 1.96
B_BOOT = 12


def _load_imp():
    path = ROOT / "scripts" / "vimp_shift_impossibility.py"
    spec = importlib.util.spec_from_file_location("vimp_imp", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_real():
    packs = []
    ds = ROOT / "datasets"

    # air quality: time order
    try:
        aq = ds / "Air_Quality_Data"
        frames = [pd.read_csv(aq / f) for f in sorted(os.listdir(aq)) if f.endswith(".csv")]
        df = pd.concat(frames, ignore_index=True)
        for col in df.columns:
            if col not in ("Date", "Time"):
                df[col] = pd.to_numeric(
                    df[col].astype(str).str.replace(";", "").str.replace(",", "."),
                    errors="coerce",
                )
        df = df.rename(columns={"RH": "Y"})
        x_cols = [
            "CO(GT)", "PT08.S1(CO)", "NMHC(GT)", "C6H6(GT)", "PT08.S2(NMHC)",
            "NOx(GT)", "PT08.S3(NOx)", "NO2(GT)", "PT08.S4(NO2)", "PT08.S5(O3)",
        ]
        med = df.loc[df["Y"] != -200, "Y"].median()
        df["Y"] = np.where(df["Y"] == -200, med, df["Y"])
        X, Y, idx = cap_xy(df[x_cols].to_numpy(float), df["Y"].to_numpy(float))
        names = [x_cols[i] for i in idx]
        X0, Y0, X1, Y1 = split_mid(X, Y)
        packs.append(("air_quality", X0, Y0, X1, Y1, names, None))
        print("[load] air_quality", len(Y0), len(Y1), flush=True)
    except Exception as e:
        print("[skip] air_quality", e, flush=True)

    try:
        df = pd.read_csv(ds / "bank" / "bank-additional-full.csv", sep=";")
        y = df["duration"].to_numpy(float)
        feat = pd.get_dummies(df.drop(columns=["duration", "y"]), drop_first=True)
        names0 = feat.columns.tolist()
        X, Y, idx = cap_xy(feat.to_numpy(float), y)
        names = [names0[i] for i in idx]
        X0, Y0, X1, Y1 = split_mid(X, Y)
        packs.append(("bank", X0, Y0, X1, Y1, names, None))
        print("[load] bank", len(Y0), len(Y1), flush=True)
    except Exception as e:
        print("[skip] bank", e, flush=True)

    try:
        src = pd.read_csv(ds / "diabetes" / "source_DiabetesReadmission.csv")
        tgt = pd.read_csv(ds / "diabetes" / "target_DiabetesReadmission.csv")
        ycol = "readmitted"
        xcols = [c for c in src.columns if c != ycol]
        X0r = src[xcols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        Y0r = pd.to_numeric(src[ycol], errors="coerce").to_numpy(float)
        X1r = tgt[xcols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        Y1r = pd.to_numeric(tgt[ycol], errors="coerce").to_numpy(float)
        X0, Y0, idx = cap_xy(X0r, Y0r)
        X1 = np.asarray(X1r, float)
        X1 = X1[:, idx] if X1.shape[1] > max(idx) else X1
        X1, Y1, _ = cap_xy(X1, Y1r, p_max=X0.shape[1])
        p = min(X0.shape[1], X1.shape[1])
        X0, X1 = X0[:, :p], X1[:, :p]
        names = [xcols[int(i)] for i in idx[:p]]
        packs.append(("diabetes", X0, Y0, X1, Y1, names, None))
        print("[load] diabetes", len(Y0), len(Y1), flush=True)
    except Exception as e:
        print("[skip] diabetes", e, flush=True)

    try:
        slice_dir = ds / "Daily_Slice"
        files = sorted(f for f in os.listdir(slice_dir) if f.endswith(".csv"))
        pick = np.linspace(0, len(files) - 1, min(6, len(files)), dtype=int)
        frames = [pd.read_csv(slice_dir / files[i]) for i in pick]
        df = pd.concat(frames, ignore_index=True)
        df["tpep_pickup_datetime"] = pd.to_datetime(df["tpep_pickup_datetime"])
        df = df.sort_values("tpep_pickup_datetime", kind="mergesort")
        x_cols = ["passenger_count", "trip_distance", "PULocationID", "DOLocationID"]
        for col in x_cols + ["total_amount"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df[col] = df[col].fillna(df[col].median())
        X, Y, idx = cap_xy(df[x_cols].to_numpy(float), df["total_amount"].to_numpy(float))
        names = [x_cols[i] for i in idx]
        X0, Y0, X1, Y1 = split_mid(X, Y)
        packs.append(("daily_slice", X0, Y0, X1, Y1, names, None))
        print("[load] daily_slice", len(Y0), len(Y1), flush=True)
    except Exception as e:
        print("[skip] daily_slice", e, flush=True)

    try:
        m5 = ds / "reformulated_M5"
        files = sorted(f for f in os.listdir(m5) if f.endswith(".csv"))[:12]
        df = pd.concat([pd.read_csv(m5 / f) for f in files], ignore_index=True)
        df = df.dropna(axis=1, how="any")
        ycol = "demand" if "demand" in df.columns else df.columns[-1]
        x_cols = [c for c in df.columns if c not in ("date", ycol) and pd.api.types.is_numeric_dtype(df[c])]
        X, Y, idx = cap_xy(df[x_cols].to_numpy(float), df[ycol].to_numpy(float))
        names = [x_cols[i] for i in idx]
        X0, Y0, X1, Y1 = split_mid(X, Y)
        packs.append(("m5", X0, Y0, X1, Y1, names, None))
        print("[load] m5", len(Y0), len(Y1), flush=True)
    except Exception as e:
        print("[skip] m5", e, flush=True)

    return packs


def cap_xy(X, Y, n=N_CAP, p_max=16):
    X = np.asarray(X, float)
    Y = np.asarray(Y, float).ravel()
    ok = np.isfinite(Y) & np.all(np.isfinite(X), axis=1)
    X, Y = X[ok], Y[ok]
    if X.shape[1] > p_max:
        var = np.nanvar(X, axis=0)
        keep = np.argsort(var)[::-1][:p_max]
        X = X[:, keep]
        names_idx = keep
    else:
        names_idx = np.arange(X.shape[1])
    if len(Y) > 2 * n:
        X, Y = X[: 2 * n], Y[: 2 * n]
    return X, Y, names_idx


def split_mid(X, Y):
    m = len(Y) // 2
    return X[:m], Y[:m], X[m:], Y[m:]


def pfi(X, Y, seed, n_perm=N_PERM):
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


def mse_pair(X0, Y0, X1, Y1, seed):
    m = Ridge(alpha=1.0, random_state=seed).fit(X0, Y0)
    mse0 = float(np.mean((Y0 - m.predict(X0)) ** 2))
    mse1 = float(np.mean((Y1 - m.predict(X1)) ** 2))
    return mse1 - mse0, mse0, mse1


def feature_localize(X0, Y0, X1, Y1, seed, names, true_idx=None):
    p0 = pfi(X0, Y0, seed + 11)
    p1 = pfi(X1, Y1, seed + 13)
    degree = p1 - p0
    order = np.argsort(-np.abs(degree))
    top = order[:TOP_K]
    hit = np.nan
    auc = np.nan
    if true_idx is not None:
        true = set(int(i) for i in true_idx)
        hit = float(sum(int(j in true) for j in top)) / TOP_K
        y = np.array([1 if j in true else 0 for j in range(len(degree))])
        if y.min() != y.max():
            auc = float(roc_auc_score(y, np.abs(degree)))
    rows = []
    for r, j in enumerate(order):
        rows.append(dict(
            rank=r + 1,
            feature=str(names[j]) if j < len(names) else f"X{j}",
            j=int(j),
            degree=float(degree[j]),
            pfi_exist=float(p0[j]),
            pfi_new=float(p1[j]),
            is_true=int(true_idx is not None and int(j) in set(int(i) for i in true_idx)),
        ))
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


def subset_on_driver(X0, Y0, X1, Y1, B, j_star, driver_idx, seed):
    """Median split on driver feature → two subsets; ΔVIMP and ΔMSE on each."""
    #With the bootstrap procedure for the subset level mean difference
    col = np.concatenate([X0[:, j_star], X1[:, j_star]])
    thr = float(np.median(col))
    rows = []
    for label, op in (("low", np.less), ("high", np.greater_equal)):
        m0 = op(X0[:, j_star], thr)
        m1 = op(X1[:, j_star], thr)
        n0, n1 = int(m0.sum()), int(m1.sum())
        if n0 < N_MIN or n1 < N_MIN:
            rows.append(dict(
                subset=label, n0=n0, n1=n1, thr=thr,
                dvimp=np.nan, dvimp_lo=np.nan, dvimp_hi=np.nan,
                dmse=np.nan, skipped="too_small",
            ))
            continue
        boots = []
        for b in range(B):
            p0 = pfi(X0[m0], Y0[m0], seed + 200 + b)
            p1 = pfi(X1[m1], Y1[m1], seed + 400 + b)
            boots.append(float(np.mean((p1 - p0)[driver_idx])))
        boots = np.asarray(boots, float)
        dv = float(boots.mean())
        se = float(boots.std(ddof=1) / np.sqrt(max(len(boots), 1)))
        dmse, _, _ = mse_pair(X0[m0], Y0[m0], X1[m1], Y1[m1], seed + 7)
        rows.append(dict(
            subset=label, n0=n0, n1=n1, thr=thr,
            dvimp=dv, dvimp_lo=dv - Z95 * se, dvimp_hi=dv + Z95 * se,
            dmse=float(dmse), skipped="",
        ))
    sdf = pd.DataFrame(rows)
    if sdf.dvimp.notna().sum() == 2:
        hi = float(sdf.loc[sdf.subset == "high", "dvimp"].iloc[0])
        lo = float(sdf.loc[sdf.subset == "low", "dvimp"].iloc[0])
        contrast = abs(hi) - abs(lo)
        flag = "high" if abs(hi) >= abs(lo) else "low"
    else:
        contrast, flag = np.nan, ""
    return sdf, float(contrast) if contrast == contrast else np.nan, flag


def run_one(name, X0, Y0, X1, Y1, names, seed, true_idx=None):
    feat, meta, _ = feature_localize(X0, Y0, X1, Y1, seed, names, true_idx=true_idx)
    feat.insert(0, "dataset", name)
    sub, contrast, flag = subset_on_driver(
        X0, Y0, X1, Y1, meta["top_j"], meta["driver_idx"], seed
    )
    sub.insert(0, "dataset", name)
    summary = dict(
        dataset=name, n0=int(len(Y0)), n1=int(len(Y1)), p=int(X0.shape[1]),
        top_feature=meta["top_name"], top_degree=meta["top_degree"],
        energy=meta["energy"], hit_at_k=meta["hit_at_k"], auc_abs=meta["auc_abs"],
        subset_contrast=contrast, flagged_subset=flag,
    )
    if len(sub) == 2 and sub.dvimp.notna().all():
        for _, r in sub.iterrows():
            summary[f"dvimp_{r.subset}"] = float(r.dvimp)
            summary[f"dvimp_{r.subset}_lo"] = float(r.dvimp_lo)
            summary[f"dvimp_{r.subset}_hi"] = float(r.dvimp_hi)
            summary[f"dmse_{r.subset}"] = float(r.dmse)
            summary[f"n_{r.subset}_0"] = int(r.n0)
            summary[f"n_{r.subset}_1"] = int(r.n1)
    return feat, sub, summary


def plot_subset_bars(summ: pd.DataFrame, path: Path):
    sub = summ.dropna(subset=["dvimp_high", "dvimp_low"])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    x = np.arange(len(sub))
    ax.bar(x - 0.18, sub.dvimp_low, 0.36, label="subset low", color="0.65")
    ax.bar(x + 0.18, sub.dvimp_high, 0.36, label="subset high", color="C0")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(sub.dataset, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel(r"two-batch $\Delta$VIMP on driver block")
    ax.set_title("Post-hoc localization: VIMP gap on two subsets of the top driver")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_tex(summ: pd.DataFrame, feat: pd.DataFrame, path: Path):
    lines = [
        r"% auto: scripts/vimp_posthoc_localization.py",
        r"\begin{table}[H]",
        r"\centering\small",
        r"\caption{Feature localization: top driver by $|\Delta$VIMP$|$ (two batches).}",
        r"\begin{tabular}{@{}llrrrr@{}}",
        r"\toprule",
        r"Dataset & top feature & $\Delta$VIMP & energy & hit@4 & AUC$_\Delta$ \\",
        r"\midrule",
    ]
    for _, r in summ.iterrows():
        hit = "---" if pd.isna(r.hit_at_k) else f"{r.hit_at_k:.2f}"
        auc = "---" if pd.isna(r.auc_abs) else f"{r.auc_abs:.2f}"
        feat_name = str(r.top_feature).replace("_", r"\_")[:22]
        lines.append(
            f"{str(r.dataset).replace('_', r'\_')} & {feat_name} & "
            f"{r.top_degree:.3f} & {r.energy:.3f} & {hit} & {auc} \\\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\begin{table}[H]",
        r"\centering\small",
        r"\caption{Subset localization: two-batch $\Delta$VIMP on low vs high of the top driver.",
        r"Flagged subset = larger $|\Delta$VIMP$|$.}",
        r"\begin{tabular}{@{}llrrrr@{}}",
        r"\toprule",
        r"Dataset & flag & $\Delta$VIMP$_{\mathrm{low}}$ & $\Delta$VIMP$_{\mathrm{high}}$ "
        r"& contrast & $\Delta$MSE$_{\mathrm{flag}}$ \\",
        r"\midrule",
    ]
    for _, r in summ.iterrows():
        if pd.isna(r.get("dvimp_high", np.nan)):
            continue
        flag = str(r.flagged_subset)
        dmse_f = r.dmse_high if flag == "high" else r.dmse_low
        lines.append(
            f"{str(r.dataset).replace('_', r'\_')} & {flag} & "
            f"{r.dvimp_low:.3f} & {r.dvimp_high:.3f} & "
            f"{r.subset_contrast:.3f} & {dmse_f:.3f} \\\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ]
    # synth top-4 feature table
    synth = feat[feat.dataset.str.startswith("synth_") & (feat["rank"] <= 4)]
    if not synth.empty:
        lines += [
            r"\begin{table}[H]",
            r"\centering\small",
            r"\caption{Injected DGPs: top-4 features by $|\Delta$VIMP$|$ (true alphas $=0$--$9$).}",
            r"\begin{tabular}{@{}llrrr@{}}",
            r"\toprule",
            r"DGP & feature & rank & $\Delta$VIMP & true \\",
            r"\midrule",
        ]
        for _, r in synth.iterrows():
            lines.append(
                f"{str(r.dataset).replace('_', r'\_')} & {r.feature} & "
                f"{int(r['rank'])} & {r.degree:.3f} & {int(r.is_true)} \\\\"
            )
        lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    path.write_text("\n".join(lines) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    imp = _load_imp()
    packs = []
    true = list(imp.ACTIVE)
    names = [f"X{j}" for j in range(imp.P)]
    X0, Y0, X1, Y1 = imp.generate_shift_pair(delta_beta=0.30, delta_gamma=0.0, seed=SEED)
    packs.append(("synth_cd", X0, Y0, X1, Y1, names, true))
    print("[load] synth_cd", flush=True)
    X0, Y0, X1, Y1 = imp.generate_shift_pair(delta_beta=0.0, delta_gamma=0.30, seed=SEED + 1)
    packs.append(("synth_cs", X0, Y0, X1, Y1, names, true))
    print("[load] synth_cs", flush=True)
    packs.extend(load_real())

    feat_all, sub_all, summaries = [], [], []
    for i, (name, X0, Y0, X1, Y1, names, true_idx) in enumerate(packs):
        print(f"[{name}] n0={len(Y0)} n1={len(Y1)} p={X0.shape[1]}", flush=True)
        feat, sub, summ = run_one(
            name, X0, Y0, X1, Y1, names, SEED + 17 * i, true_idx=true_idx
        )
        feat_all.append(feat)
        sub_all.append(sub)
        summaries.append(summ)
        print(
            f"  top={summ['top_feature']} ΔVIMP={summ['top_degree']:.3f} "
            f"flag={summ['flagged_subset']} contrast={summ['subset_contrast']}",
            flush=True,
        )

    feat_df = pd.concat(feat_all, ignore_index=True)
    sub_df = pd.concat(sub_all, ignore_index=True)
    summ_df = pd.DataFrame(summaries)
    feat_df.to_csv(OUT / "feature_rank.csv", index=False)
    sub_df.to_csv(OUT / "subset_vimp.csv", index=False)
    summ_df.to_csv(OUT / "localization_summary.csv", index=False)
    plot_subset_bars(summ_df, OUT / "subset_vimp_bars.png")
    write_tex(summ_df, feat_df, OUT / "localization_tables.tex")

    lines = [
        "# Post-hoc VIMP localization",
        "",
        "After impossibility (cannot decompose CD vs CS), localize:",
        "1. features with large two-batch ΔVIMP",
        "2. the observation subset (low vs high of the top driver) where |ΔVIMP| is larger",
        "",
        summ_df.to_string(index=False),
        "",
        "One-click:",
        "```bash",
        "PYTHONPATH=src:scripts python scripts/vimp_posthoc_localization.py",
        "```",
    ]
    (OUT / "LOCALIZATION_NOTE.md").write_text("\n".join(lines) + "\n")
    payload = summ_df.to_dict(orient="records")
    (OUT / "payload.json").write_text(json.dumps(payload, indent=2, default=str))
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
