#!/usr/bin/env python3
"""AFFEC multimodal FSDS · feature selection for distribution shift.

Unzipped streams under data/affec/ → align modality channel indices →
subsample ~10k windows → CFPerm-style batch VIMP (permute-W null) →
emit selected feature-indices per modality.

  python3 scripts/run_affec_fsds.py
"""
from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "affec"
OUT = ROOT / "results" / "affec_fsds"
DOCS = ROOT / "docs" / "method"

SEED = 2026
N_SUBSAMPLE = 10_000
N_PERM = 20  # permute-W null for MMD-LOCO (board; not a method paper)
ALPHA = 0.05
WIN_SEC = 1.0  # window length for EEG / physio aggregation
MAX_WIN_PHYSIO = 24
MAX_WIN_EEG = 16
MMD_MAX_N = 400
# RF Domain Classifier board: ranked indices per modality
TOP_K = 8
SKIP_MMD = os.environ.get("AFFEC_FSDS_SKIP_MMD", "0") == "1"


# ---------------------------------------------------------------------------
# Channel Name Attribution (user-provided index map)
# ---------------------------------------------------------------------------
CHANNEL_NAMES = {
    "eye_tracking": {
        0: ("fixation_x", "fixation point X coordinate"),
        1: ("fixation_y", "fixation point Y coordinate"),
        2: ("gaze_x_left", "left eye gaze direction X"),
        3: ("gaze_y_left", "left eye gaze direction Y"),
        4: ("gaze_x_right", "right eye gaze direction X"),
        5: ("gaze_y_right", "right eye gaze direction Y"),
        6: ("fixation_duration", "fixation duration"),
        7: ("blink_left", "left eye blink flag"),
        8: ("blink_right", "right eye blink flag"),
        9: ("eye_open_left", "left eye openness"),
        10: ("eye_open_right", "right eye openness"),
        11: ("pupil_x_left", "left pupil X position"),
        12: ("pupil_y_left", "left pupil Y position"),
        13: ("pupil_x_right", "right pupil X position"),
        14: ("pupil_y_right", "right pupil Y position"),
        15: ("validity", "eye-tracking validity flag"),
    },
    "pupil": {
        i: (name, desc)
        for i, (name, desc) in {
            0: ("pupil_diameter", "pupil diameter"),
            1: ("pupil_diameter_raw", "raw pupil diameter"),
            2: ("pupil_diameter_filt", "filtered pupil diameter"),
            3: ("pupil_baseline", "pupil baseline diameter"),
            4: ("pupil_peak", "pupil peak diameter"),
            5: ("pupil_latency", "time to pupil peak"),
            6: ("pupil_velocity", "pupil change velocity"),
            7: ("pupil_acceleration", "pupil change acceleration"),
            8: ("pupil_variability", "pupil variability"),
            9: ("pupil_slope", "pupil change slope"),
            10: ("eye_pos_x", "eye position X"),
            11: ("eye_pos_y", "eye position Y"),
            12: ("eye_pos_x_raw", "raw eye position X"),
            13: ("eye_pos_y_raw", "raw eye position Y"),
            14: ("eye_pos_velocity_x", "eye position velocity X"),
            15: ("eye_pos_velocity_y", "eye position velocity Y"),
            16: ("eye_pos_dispersion", "eye position dispersion"),
            17: ("eye_pos_range_x", "eye position range X"),
            18: ("eye_pos_range_y", "eye position range Y"),
            19: ("eye_pos_slope", "eye position slope"),
            20: ("pupil_validity", "pupil data validity flag"),
        }.items()
    },
    "cursor": {
        0: ("cursor_x", "cursor X coordinate"),
        1: ("cursor_y", "cursor Y coordinate"),
        2: ("cursor_velocity", "cursor movement velocity"),
        3: ("cursor_state", "cursor state (click/hover)"),
    },
    "gsr_eda": {
        **{
            0: ("gsr_raw", "raw skin conductance"),
            1: ("gsr_filtered", "filtered skin conductance"),
            2: ("gsr_phasic", "phasic component (instant response)"),
            3: ("gsr_tonic", "tonic component (background stress)"),
            4: ("gsr_scr_count", "skin conductance response count"),
            5: ("gsr_peaks", "skin conductance peak count"),
            6: ("gsr_slope", "skin conductance slope"),
            7: ("gsr_mean", "skin conductance mean"),
            8: ("gsr_std", "skin conductance std"),
            9: ("gsr_range", "skin conductance range"),
            10: ("body_temp", "body temperature"),
            11: ("temp_mean", "body temperature mean"),
            12: ("temp_std", "body temperature std"),
            13: ("temp_slope", "body temperature slope"),
            14: ("acc_x", "accelerometer X"),
            15: ("acc_y", "accelerometer Y"),
            16: ("acc_z", "accelerometer Z"),
            17: ("acc_magnitude", "accelerometer magnitude"),
            18: ("acc_std", "accelerometer std"),
            19: ("acc_slope", "accelerometer slope"),
        },
        **{i: (f"gsr_feat_{i}", f"GSR feature {i}") for i in range(20, 40)},
    },
    "eeg": {
        **{
            0: ("Fp1", "prefrontal left (executive function / emotion regulation)"),
            1: ("Fp2", "prefrontal right (executive function / emotion regulation)"),
            2: ("F3", "left frontal (decision / planning)"),
            3: ("F4", "right frontal (decision / planning)"),
            4: ("C3", "left central (motor cortex)"),
            5: ("C4", "right central (motor cortex)"),
            6: ("P3", "left parietal (spatial attention)"),
            7: ("P4", "right parietal (spatial attention)"),
            8: ("O1", "left occipital (visual processing)"),
            9: ("O2", "right occipital (visual processing)"),
        },
        **{i: (f"EEG_{i}", f"EEG channel {i}") for i in range(10, 63)},
    },
}

# AFFEC raw column → user feature index (within modality)
GAZE_COL_TO_IDX = {
    "FPOGX": 0,
    "FPOGY": 1,
    "LPOGX": 2,
    "LPOGY": 3,
    "RPOGX": 4,
    "RPOGY": 5,
    "FPOGD": 6,
    "LPOGV": 9,  # openness/validity proxy
    "RPOGV": 10,
    "BPOGX": 11,  # best-eye as pupil-x proxy slot when pupil stream separate
    "BPOGY": 12,
    "FPOGV": 15,
}
PUPIL_COL_TO_IDX = {
    "LPD": 0,
    "LPUPILD": 1,
    "RPD": 2,
    "RPUPILD": 3,
    "LPS": 4,
    "RPS": 5,
    "LPCX": 10,
    "LPCY": 11,
    "RPCX": 12,
    "RPCY": 13,
    "LEYEX": 14,
    "LEYEY": 15,
    "REYEX": 16,
    "REYEY": 17,
    "LPV": 19,
    "RPV": 20,
}
# GSR: pack calibrated / primary signals into first slots, rest by order
GSR_PRIORITY = [
    "GSR_raw",
    "GSR_cal",
    "GSR_Conductance_cal",
    "Temperature_cal",
    "Temperature_raw",
    "Low_Noise_Accelerometer_X_cal",
    "Low_Noise_Accelerometer_Y_cal",
    "Low_Noise_Accelerometer_Z_cal",
    "Wide_Range_Accelerometer_X_cal",
    "Wide_Range_Accelerometer_Y_cal",
    "Wide_Range_Accelerometer_Z_cal",
    "Gyroscope_X_cal",
    "Gyroscope_Y_cal",
    "Gyroscope_Z_cal",
    "Pressure_cal",
    "VSenseBatt_cal",
]


def _read_physio(tsv_gz: Path, json_path: Path) -> tuple[list[str], np.ndarray]:
    meta = json.loads(json_path.read_text())
    cols = list(meta["Columns"])
    rows = []
    with gzip.open(tsv_gz, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != len(cols):
                continue
            try:
                rows.append([float(x) if x not in ("", "n/a", "NA") else np.nan for x in parts])
            except ValueError:
                continue
    arr = np.asarray(rows, dtype=float) if rows else np.zeros((0, len(cols)))
    return cols, arr


def _window_means(arr: np.ndarray, onset: np.ndarray, win: float, rng: np.random.Generator, max_win: int = 8):
    """Sample up to max_win non-overlapping windows; return list of mean vectors over feature cols."""
    if arr.size == 0 or onset.size == 0:
        return []
    t0, t1 = float(np.nanmin(onset)), float(np.nanmax(onset))
    if not np.isfinite(t0) or t1 - t0 < win:
        mu = np.nanmean(arr, axis=0)
        return [mu] if np.isfinite(mu).any() else []
    starts = np.arange(t0, t1 - win, win)
    if len(starts) == 0:
        return []
    if len(starts) > max_win:
        starts = rng.choice(starts, size=max_win, replace=False)
    out = []
    for s in starts:
        m = (onset >= s) & (onset < s + win)
        if m.sum() < 3:
            continue
        out.append(np.nanmean(arr[m], axis=0))
    return out


def _pack_by_map(cols: list[str], arr: np.ndarray, col_to_idx: dict[str, int], n_out: int) -> np.ndarray:
    """Map named columns into fixed-length feature vector (NaN fill)."""
    out = np.full((arr.shape[0], n_out), np.nan, dtype=float)
    name_to_j = {c: j for j, c in enumerate(cols)}
    for c, idx in col_to_idx.items():
        if c in name_to_j and 0 <= idx < n_out:
            out[:, idx] = arr[:, name_to_j[c]]
    return out


def _pack_gsr(cols: list[str], arr: np.ndarray, n_out: int = 40) -> np.ndarray:
    skip = {"onset", "Timestamp_raw", "Timestamp_cal", "System_Timestamp_cal", "_raw"}
    feat_cols = [c for c in cols if c not in skip]
    # priority first, then remaining in file order
    ordered = [c for c in GSR_PRIORITY if c in feat_cols]
    ordered += [c for c in feat_cols if c not in ordered]
    ordered = ordered[:n_out]
    name_to_j = {c: j for j, c in enumerate(cols)}
    out = np.full((arr.shape[0], n_out), np.nan, dtype=float)
    for i, c in enumerate(ordered):
        out[:, i] = arr[:, name_to_j[c]]
    return out


def _pack_cursor(cols: list[str], arr: np.ndarray) -> np.ndarray:
    name_to_j = {c: j for j, c in enumerate(cols)}
    n = arr.shape[0]
    out = np.full((n, 4), np.nan, dtype=float)
    if "CX" in name_to_j:
        out[:, 0] = arr[:, name_to_j["CX"]]
    if "CY" in name_to_j:
        out[:, 1] = arr[:, name_to_j["CY"]]
    # velocity from successive cursor positions
    if np.isfinite(out[:, 0]).any():
        dx = np.diff(out[:, 0], prepend=out[0, 0])
        dy = np.diff(out[:, 1], prepend=out[0, 1])
        out[:, 2] = np.sqrt(dx * dx + dy * dy)
    if "CS" in name_to_j:
        out[:, 3] = arr[:, name_to_j["CS"]]
    return out


def _load_eeg_window_means(edf: Path, rng: np.random.Generator, max_win: int = 6) -> list[np.ndarray]:
    import mne

    raw = mne.io.read_raw_edf(edf, preload=True, verbose=False)
    data = raw.get_data()  # (n_ch, n_times)
    sfreq = float(raw.info["sfreq"])
    n_ch = min(63, data.shape[0])
    data = data[:n_ch]
    win = int(WIN_SEC * sfreq)
    if win < 8 or data.shape[1] < win:
        mu = np.nanmean(data, axis=1)
        vec = np.full(63, np.nan)
        vec[: n_ch] = mu
        return [vec]
    n_possible = data.shape[1] // win
    starts = np.arange(0, n_possible * win, win)
    if len(starts) > max_win:
        starts = rng.choice(starts, size=max_win, replace=False)
    out = []
    for s in starts:
        chunk = data[:, int(s) : int(s) + win]
        vec = np.full(63, np.nan)
        vec[:n_ch] = np.nanmean(chunk, axis=1)
        out.append(vec)
    return out


def collect_windows(rng: np.random.Generator):
    """Build multimodal window table keyed by (subject, run, window_id)."""
    # gather files by (sub, run)
    physio_kinds = {
        "gaze": "eye_tracking",
        "pupil": "pupil",
        "cursor": "cursor",
        "gsr": "gsr_eda",
    }
    records = []  # list of dicts with modality vectors + meta

    # physio
    for kind, mod in physio_kinds.items():
        for tsv in DATA.rglob(f"*recording-{kind}_physio.tsv.gz"):
            m = re.search(r"sub-([^/]+).*run-(\d+)", str(tsv))
            if not m:
                continue
            sub, run = m.group(1), int(m.group(2))
            js = tsv.with_suffix("").with_suffix(".json")
            # path is .tsv.gz → replace
            js = Path(str(tsv).replace(".tsv.gz", ".json"))
            if not js.exists():
                continue
            cols, arr = _read_physio(tsv, js)
            if arr.shape[0] < 10:
                continue
            onset = arr[:, cols.index("onset")] if "onset" in cols else np.arange(arr.shape[0])
            if kind == "gaze":
                feat = _pack_by_map(cols, arr, GAZE_COL_TO_IDX, 16)
            elif kind == "pupil":
                feat = _pack_by_map(cols, arr, PUPIL_COL_TO_IDX, 21)
            elif kind == "cursor":
                feat = _pack_cursor(cols, arr)
            else:
                feat = _pack_gsr(cols, arr, 40)
            for wi, mu in enumerate(
                _window_means(feat, onset, WIN_SEC, rng, max_win=MAX_WIN_PHYSIO)
            ):
                records.append(
                    dict(subject=sub, run=run, window=wi, modality=mod, vec=mu.astype(float))
                )

    # eeg
    for edf in DATA.rglob("*_eeg.edf"):
        m = re.search(r"sub-([^/]+).*run-(\d+)", str(edf))
        if not m:
            continue
        sub, run = m.group(1), int(m.group(2))
        try:
            vecs = _load_eeg_window_means(edf, rng, max_win=MAX_WIN_EEG)
        except Exception as e:
            print(f"skip eeg {edf.name}: {e}", flush=True)
            continue
        for wi, mu in enumerate(vecs):
            records.append(dict(subject=sub, run=run, window=wi, modality="eeg", vec=mu.astype(float)))

    return records


def pivot_multimodal(records: list[dict], rng: np.random.Generator):
    """Align modalities on (subject, run, window); fill missing mods with run-level mean."""
    mods = ["eye_tracking", "pupil", "cursor", "gsr_eda", "eeg"]
    dims = {m: max(CHANNEL_NAMES[m]) + 1 for m in mods}
    by_key: dict[tuple, dict] = {}
    run_pool: dict[tuple, dict] = {}  # (sub, run, mod) -> list of vecs
    for r in records:
        key = (r["subject"], r["run"], r["window"])
        by_key.setdefault(key, {})[r["modality"]] = r["vec"]
        rp = (r["subject"], r["run"], r["modality"])
        run_pool.setdefault(rp, []).append(r["vec"])
    run_mean = {k: np.nanmean(np.vstack(v), axis=0) for k, v in run_pool.items()}

    rows_X, rows_W, rows_Y, rows_meta = [], [], [], []
    block_slices = {}
    start = 0
    for m in mods:
        block_slices[m] = (start, start + dims[m])
        start += dims[m]
    p = start

    keys = list(by_key.keys())
    rng.shuffle(keys)
    for key in keys:
        pack = dict(by_key[key])
        sub, run, wi = key
        # fill missing modalities from same-run mean so blocks stay interpretable
        for m in mods:
            if m not in pack and (sub, run, m) in run_mean:
                pack[m] = run_mean[(sub, run, m)]
        if len(pack) < 2:
            continue
        x = np.full(p, np.nan, dtype=float)
        for m, (a, b) in block_slices.items():
            if m in pack:
                v = pack[m]
                x[a : a + len(v)] = v[: b - a]
        if np.isfinite(x).mean() < 0.2:
            continue
        w = 0 if run <= 1 else 1
        if "pupil" in pack and np.isfinite(pack["pupil"][0]):
            y = float(pack["pupil"][0])
        elif "gsr_eda" in pack and np.isfinite(pack["gsr_eda"][0]):
            y = float(pack["gsr_eda"][0])
        elif "cursor" in pack and np.isfinite(pack["cursor"][2]):
            y = float(pack["cursor"][2])
        else:
            y = float(np.nanmean(x))
        rows_X.append(x)
        rows_W.append(w)
        rows_Y.append(y)
        rows_meta.append({"subject": sub, "run": run, "window": wi, "mods": sorted(pack.keys())})
        if len(rows_X) >= N_SUBSAMPLE * 3:
            break

    X = np.asarray(rows_X, float)
    W = np.asarray(rows_W, int)
    Y = np.asarray(rows_Y, float)
    col_med = np.nanmedian(X, axis=0)
    inds = np.where(~np.isfinite(X))
    X[inds] = np.take(col_med, inds[1])
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    Y = np.nan_to_num(Y, nan=float(np.nanmedian(Y)))

    rng2 = np.random.default_rng(SEED + 7)
    idx0 = np.where(W == 0)[0]
    idx1 = np.where(W == 1)[0]
    n_each = min(len(idx0), len(idx1), N_SUBSAMPLE // 2)
    if n_each < 100:
        raise RuntimeError(f"too few aligned windows: n0={len(idx0)} n1={len(idx1)}")
    # with-replacement if needed to hit ~10k balanced
    replace0 = len(idx0) < N_SUBSAMPLE // 2
    replace1 = len(idx1) < N_SUBSAMPLE // 2
    n_target = N_SUBSAMPLE // 2
    take = np.concatenate(
        [
            rng2.choice(idx0, n_target, replace=replace0 or len(idx0) < n_target),
            rng2.choice(idx1, n_target, replace=replace1 or len(idx1) < n_target),
        ]
    )
    rng2.shuffle(take)
    return X[take], Y[take], W[take], [rows_meta[i] for i in take], block_slices, mods


def _rbf_mmd2(X0: np.ndarray, X1: np.ndarray, *, gamma: float | None = None, max_n: int = 800, rng: np.random.Generator | None = None) -> float:
    """Unbiased squared MMD with RBF kernel (median heuristic)."""
    rng = rng or np.random.default_rng(0)
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    if len(X0) > max_n:
        X0 = X0[rng.choice(len(X0), max_n, replace=False)]
    if len(X1) > max_n:
        X1 = X1[rng.choice(len(X1), max_n, replace=False)]
    Z = np.vstack([X0, X1])
    # median heuristic on a subsample of pairwise distances
    if gamma is None:
        idx = rng.choice(len(Z), size=min(400, len(Z)), replace=False)
        S = Z[idx]
        d2 = np.sum((S[:, None, :] - S[None, :, :]) ** 2, axis=-1)
        med = float(np.median(d2[d2 > 0])) if np.any(d2 > 0) else 1.0
        gamma = 1.0 / max(med, 1e-8)

    def k(A, B):
        d2 = np.sum((A[:, None, :] - B[None, :, :]) ** 2, axis=-1)
        return np.exp(-gamma * d2)

    K00 = k(X0, X0)
    K11 = k(X1, X1)
    K01 = k(X0, X1)
    n0, n1 = len(X0), len(X1)
    # unbiased U-statistic diagonals removed
    mmd2 = (
        (K00.sum() - np.trace(K00)) / max(n0 * (n0 - 1), 1)
        + (K11.sum() - np.trace(K11)) / max(n1 * (n1 - 1), 1)
        - 2.0 * K01.mean()
    )
    return float(mmd2)


def po_risk_learner(X, Y, W, *, seed: int = 0, n_splits: int = 5, clip_e: float = 0.01):
    """Cross-fit PO-risk learner (DR-style): m(X), e(X), τ on pseudo-outcome.

    Returns observed PO-risk, cross-fit nuisances, and per-sample PO scores.
    """
    from sklearn.model_selection import StratifiedKFold

    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    W = np.asarray(W, int)
    n, p = X.shape
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=60, max_depth=6, min_samples_leaf=5, random_state=seed + fold, n_jobs=-1
        )
        e = RandomForestClassifier(
            n_estimators=60, max_depth=6, min_samples_leaf=5, random_state=seed + 100 + fold, n_jobs=-1
        )
        m.fit(X[tr], Y[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, clip_e, 1.0 - clip_e)
    # Pseudo-outcome / PO score (DR residual product)
    po_score = (Y - m_hat) * (W - e_hat)
    # τ-learner on X → PO-risk = E[τ(X)^2]
    tau = RandomForestRegressor(
        n_estimators=80, max_depth=6, min_samples_leaf=5, random_state=seed + 7, n_jobs=-1
    )
    tau.fit(X, po_score)
    tau_hat = tau.predict(X)
    po_risk = float(np.mean(tau_hat ** 2))
    return dict(po_risk=po_risk, m_hat=m_hat, e_hat=e_hat, po_score=po_score, tau_hat=tau_hat)


def mmd_loco_importance(X, W, *, seed: int = 0, max_n: int = MMD_MAX_N) -> np.ndarray:
    """MMD-LOCO: Δ_j = MMD²(P0,P1) - MMD²(P0^{-j}, P1^{-j}).

    Large Δ_j ⇒ feature j drives the between-batch discrepancy.
    """
    rng = np.random.default_rng(seed)
    X = np.asarray(X, float)
    W = np.asarray(W, int)
    X0, X1 = X[W == 0], X[W == 1]
    mmd_full = _rbf_mmd2(X0, X1, max_n=max_n, rng=rng)
    p = X.shape[1]
    imps = np.zeros(p, float)
    for j in range(p):
        cols = [c for c in range(p) if c != j]
        imps[j] = mmd_full - _rbf_mmd2(X0[:, cols], X1[:, cols], max_n=max_n, rng=rng)
        if (j + 1) % 20 == 0:
            print(f"  MMD-LOCO {j+1}/{p}", flush=True)
    return imps


def fsds_porisk_mmd_loco(X, Y, W, *, n_perm: int = N_PERM, seed: int = SEED):
    """FSDS = PO-risk learner + MMD-LOCO (the locked attribution logic).

    1) PO-risk learner: global discrepancy between batches via frozen-style DR PO-risk.
    2) MMD-LOCO on X: which coordinates drive the P0 vs P1 shift.
    3) Null: permute W, recompute MMD-LOCO → feature p-values (CFPerm-style).
    """
    rng = np.random.default_rng(seed)
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    W = np.asarray(W, int)
    n, p = X.shape

    print("  PO-risk learner…", flush=True)
    po = po_risk_learner(X, Y, W, seed=seed)
    print(f"  observed PO-risk={po['po_risk']:.6f}", flush=True)

    print("  MMD-LOCO importance…", flush=True)
    imp_obs = mmd_loco_importance(X, W, seed=seed)
    # mild coupling: upweight features that also move the PO score across batches
    po_s = po["po_score"]
    for j in range(p):
        ms = abs(np.mean(X[W == 1, j]) - np.mean(X[W == 0, j]))
        s = np.std(X[:, j]) + 1e-8
        # keep MMD-LOCO primary; PO-risk path as light coupler
        imp_obs[j] = float(imp_obs[j]) + 0.05 * (ms / s) * abs(po["po_risk"])

    perm = np.zeros((p, n_perm), float)
    for b in range(n_perm):
        Wb = rng.permutation(W)
        perm[:, b] = mmd_loco_importance(X, Wb, seed=seed + 17 + b)
        if (b + 1) % 5 == 0:
            print(f"  perm {b+1}/{n_perm}", flush=True)

    pvals = (1.0 + np.sum(perm >= imp_obs[:, None], axis=1)) / (1.0 + n_perm)
    q_feat = np.quantile(perm, 0.95, axis=1)
    thr = float(np.quantile(q_feat, 0.90))
    rejected = [int(j) for j in range(p) if imp_obs[j] > thr and pvals[j] <= ALPHA]
    if len(rejected) == 0:
        rejected = [int(j) for j in np.argsort(pvals) if pvals[j] <= ALPHA][:30]
    return dict(
        imp=imp_obs,
        pvals=pvals,
        threshold=thr,
        rejected=rejected,
        po_risk=po["po_risk"],
        method="PO-risk learner + MMD-LOCO",
    )


def rf_domain_classifier_vimp(X, W, *, seed: int = SEED) -> tuple[np.ndarray, float]:
    """RF Domain Classifier: predict batch W from X → impurity VIMP.

    Ranking for the business board: np.argsort(-VIMP)[:k] (within each modality).
    """
    X = np.asarray(X, float)
    W = np.asarray(W, int)
    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=10,
        min_samples_leaf=5,
        random_state=seed,
        n_jobs=-1,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    # quick holdout AUC as board sanity
    rng = np.random.default_rng(seed)
    idx = np.arange(len(W))
    rng.shuffle(idx)
    cut = int(0.8 * len(W))
    tr, te = idx[:cut], idx[cut:]
    clf2 = RandomForestClassifier(
        n_estimators=120, max_depth=10, min_samples_leaf=5, random_state=seed + 1, n_jobs=-1
    )
    clf2.fit(X[tr], W[tr])
    proba = clf2.predict_proba(X[te])[:, 1]
    auc = float(roc_auc_score(W[te], proba)) if len(np.unique(W[te])) > 1 else float("nan")
    return vimp, auc


def ranked_indices_by_modality(vimp: np.ndarray, block_slices: dict, mods: list[str], *, k: int = TOP_K):
    """Business board dict: modality → indices ordered by RF Domain Classifier VIMP.

    order = np.argsort(-VIMP_mod)[:k]
    """
    vimp = np.asarray(vimp, float)
    out = {}
    detail = {}
    for m in mods:
        a, b = block_slices[m]
        vm = vimp[a:b]
        order = np.argsort(-vm)[: min(k, len(vm))]
        ranked = [int(i) for i in order]
        out[m] = ranked  # ORDER MATTERS
        detail[m] = [
            {
                "rank": r + 1,
                "index": int(idx),
                "global_index": int(a + idx),
                "name": CHANNEL_NAMES[m].get(idx, (f"feat_{idx}", ""))[0],
                "desc": CHANNEL_NAMES[m].get(idx, ("", ""))[1],
                "rf_vimp": float(vm[idx]),
            }
            for r, idx in enumerate(ranked)
        ]
    return out, detail


def map_rejected_to_modalities(rejected, block_slices, mods, pvals, imp):
    """MMD-LOCO selected set (unordered support); RF ranking is the board order."""
    out = {}
    detail = {}
    for m in mods:
        a, b = block_slices[m]
        local = sorted(j - a for j in rejected if a <= j < b)
        out[m] = local
        detail[m] = [
            {
                "index": int(idx),
                "global_index": int(a + idx),
                "name": CHANNEL_NAMES[m].get(idx, (f"feat_{idx}", ""))[0],
                "desc": CHANNEL_NAMES[m].get(idx, ("", ""))[1],
                "pval": float(pvals[a + idx]),
                "imp": float(imp[a + idx]),
            }
            for idx in local
        ]
    return out, detail


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    if not DATA.exists():
        raise SystemExit(f"missing extracted AFFEC data at {DATA}")

    cache = OUT / "affec_fsds_xyw_cache.npz"
    rng = np.random.default_rng(SEED)
    if cache.exists():
        print(f"loading cache {cache}", flush=True)
        z = np.load(cache, allow_pickle=True)
        X, Y, W = z["X"], z["Y"], z["W"]
        block_slices = {k: tuple(v) for k, v in z["block_slices"].item().items()}
        mods = list(z["mods"])
    else:
        print("collecting multimodal windows…", flush=True)
        records = collect_windows(rng)
        print(f"raw modality window records: {len(records)}", flush=True)
        X, Y, W, meta, block_slices, mods = pivot_multimodal(records, rng)
        np.savez_compressed(
            cache,
            X=X,
            Y=Y,
            W=W,
            block_slices=block_slices,
            mods=np.array(mods, dtype=object),
        )
        print(f"wrote cache {cache}", flush=True)

    print(
        f"aligned subsample: n={len(W)} p={X.shape[1]} "
        f"n0={(W==0).sum()} n1={(W==1).sum()}",
        flush=True,
    )

    # --- RF Domain Classifier: ordered board (primary prototype dict) ---
    print("RF Domain Classifier VIMP + argsort(-VIMP)[:k]…", flush=True)
    rf_vimp, rf_auc = rf_domain_classifier_vimp(X, W, seed=SEED)
    feature_indices, detail = ranked_indices_by_modality(
        rf_vimp, block_slices, mods, k=TOP_K
    )
    print(f"  RF domain AUC={rf_auc:.3f}", flush=True)

    compact = "feature-indices (RF Domain Classifier rank): " + ", ".join(
        f"{m}:{feature_indices[m]}" for m in mods
    )
    print(compact, flush=True)

    payload = {
        "method": "PO-risk learner + MMD-LOCO + RF Domain Classifier",
        "board": "RF Domain Classifier VIMP rank · np.argsort(-VIMP)[:k] per modality",
        "top_k": TOP_K,
        "rf_domain_auc": rf_auc,
        "n": int(len(W)),
        "p": int(X.shape[1]),
        "alpha": ALPHA,
        "batch_def": "W=0: run∈{0,1}; W=1: run∈{2,3}",
        "feature_indices": feature_indices,  # ORDERED by RF VIMP
        "feature_indices_detail": detail,
        "block_slices": {m: list(block_slices[m]) for m in mods},
        "compact": compact,
    }

    if not SKIP_MMD:
        print("running FSDS (PO-risk learner + MMD-LOCO)…", flush=True)
        if len(W) > 4000:
            rng_s = np.random.default_rng(SEED + 3)
            i0 = np.where(W == 0)[0]
            i1 = np.where(W == 1)[0]
            take = np.concatenate(
                [rng_s.choice(i0, 2000, replace=False), rng_s.choice(i1, 2000, replace=False)]
            )
            rng_s.shuffle(take)
            Xs, Ys, Ws = X[take], Y[take], W[take]
            print(f"  MMD-LOCO compute subsample n={len(Ws)} (from {len(W)})", flush=True)
        else:
            Xs, Ys, Ws = X, Y, W
        res = fsds_porisk_mmd_loco(Xs, Ys, Ws, n_perm=N_PERM, seed=SEED)
        mmd_sel, mmd_detail = map_rejected_to_modalities(
            res["rejected"], block_slices, mods, res["pvals"], res["imp"]
        )
        payload["po_risk"] = res["po_risk"]
        payload["n_mmd"] = int(len(Ws))
        payload["n_perm"] = N_PERM
        payload["mmd_loco_selected"] = mmd_sel
        payload["mmd_loco_detail"] = mmd_detail
        payload["threshold"] = res["threshold"]
    else:
        print("SKIP_MMD=1 · board is RF Domain Classifier rank only", flush=True)

    dest = OUT / "affec_fsds_feature_indices.json"
    dest.write_text(json.dumps(payload, indent=2))
    print(f"wrote {dest}", flush=True)

    lines = [
        "# AFFEC FSDS board · RF Domain Classifier rank",
        "",
        compact,
        "",
        f"- board order: `np.argsort(-RF_VIMP)[:{TOP_K}]` within each modality",
        f"- RF domain AUC: {rf_auc:.3f}",
        f"- n={payload['n']}, p={payload['p']}",
        f"- batch: {payload['batch_def']}",
        "",
        "| modality | ranked indices (high → low VIMP) | names |",
        "|---|---|---|",
    ]
    for m in mods:
        names = [d["name"] for d in detail[m]]
        lines.append(f"| {m} | {feature_indices[m]} | {', '.join(names) if names else '—'} |")
    md = OUT / "affec_fsds_feature_indices.md"
    md.write_text("\n".join(lines) + "\n")
    print(f"wrote {md}", flush=True)

    # tiny prototype dump for copy-paste
    proto = OUT / "affec_fsds_feature_indices_prototype.py"
    proto.write_text(
        "# RF Domain Classifier ranked feature-indices (order = importance)\n"
        "feature_indices = "
        + repr(feature_indices)
        + "\n"
    )
    print(f"wrote {proto}", flush=True)


if __name__ == "__main__":
    main()
