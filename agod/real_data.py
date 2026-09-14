"""Load real tables as consecutive-batch streams (no shuffle, no K-fold).

Tries local Affec / Tencent / MSR-VTT / COCO packs first, then public
sklearn / OpenML tables ordered by a slow coordinate so neighbouring
batches are temporally or spatially adjacent.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from agod.po_refit import Stream, stream_from_xy

ROOT = Path(__file__).resolve().parents[1]


def _order(X, y, key):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).ravel()
    key = np.asarray(key).ravel()
    order = np.argsort(key, kind="mergesort")
    return X[order], y[order]


def _pca(X, d, seed=0):
    X = np.asarray(X, dtype=float)
    if X.shape[1] <= d:
        return X
    from sklearn.decomposition import PCA

    return PCA(n_components=int(d), random_state=int(seed)).fit_transform(X)


def _stream(X, y, *, name, task, n_per, n_batches, meta=None):
    n = int(n_per) * int(n_batches)
    if len(X) < n:
        return None
    return stream_from_xy(
        X, y, n_per=n_per, n_batches=n_batches, name=name, task=task, meta=meta
    )


def load_affec(root: Path, max_n: int, seed: int, pca_d: int):
    cache = root / "results/affec_fsds/affec_fsds_xyw_cache.npz"
    if not cache.is_file():
        return None
    z = np.load(cache, allow_pickle=True)
    X, y = z["X"].astype(np.float32), z["Y"].astype(np.float64)
    n = min(max_n, len(X))
    return _pca(X[:n], pca_d, seed), y[:n], "mse"


def load_msrvtt(root: Path, max_n: int, seed: int, pca_d: int):
    d = root / "data/msrvtt/packed"
    if not (d / "video_feat.npy").is_file():
        return None
    X = np.concatenate(
        [
            np.load(d / "video_feat.npy"),
            np.load(d / "audio_feat.npy"),
            np.load(d / "text_feat.npy"),
        ],
        axis=1,
    )
    y = np.load(d / "labelsmsr.npy").astype(np.int64)
    Xp = _pca(X, pca_d, seed)
    order = np.argsort(Xp[:, 0])
    Xp, y = Xp[order], y[order]
    n = min(max_n, len(Xp))
    return Xp[:n], y[:n], "acc"


def load_img_txt(root: Path, name: str, max_n: int, seed: int, pca_d: int):
    d = root / "data/img_txt" / name
    if not (d / "img_feats.npy").is_file():
        return None
    X = np.concatenate(
        [np.load(d / "img_feats.npy"), np.load(d / "txt_feats.npy")], axis=1
    )
    y = np.load(d / "labels.npy")
    y = np.asarray(y).ravel().astype(int)
    u, c = np.unique(y, return_counts=True)
    keep = set(u[np.argsort(-c)[:4]].tolist())
    remap = {lab: i for i, lab in enumerate(sorted(keep))}
    y = np.array([remap[v] if v in keep else 4 for v in y], dtype=int)
    n = min(max_n, len(X))
    return _pca(X[:n], pca_d, seed), y[:n], "acc"


def load_tencent(root: Path, max_n: int, seed: int, pca_d: int):
    p = root / "results/tencent_gr/user_feats_lean.parquet"
    if not p.is_file():
        return None
    import pandas as pd

    df = pd.read_parquet(p)
    y_col = "life_ctcvr" if "life_ctcvr" in df.columns else "arpu_mean"
    drop = {"user_id", y_col}
    num = [c for c in df.columns if c not in drop and np.issubdtype(df[c].dtype, np.number)]
    X = df[num].fillna(0).to_numpy(np.float32)
    y = df[y_col].fillna(0).to_numpy(np.float64)
    n = min(max_n, len(X))
    return _pca(X[:n], pca_d, seed), y[:n], "mse"


def _as_numeric(df):
    import pandas as pd

    parts = []
    for c in df.columns:
        s = df[c]
        if str(s.dtype) in ("object", "string", "category") or s.dtype == object:
            codes, _ = pd.factorize(s.astype(str), sort=False)
            parts.append(codes.astype(np.float64)[:, None])
        else:
            v = pd.to_numeric(s, errors="coerce").to_numpy(dtype=np.float64)
            v = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
            parts.append(v[:, None])
    return np.hstack(parts) if parts else np.zeros((len(df), 1))


def load_interstate(root: Path, max_n: int, seed: int, pca_d: int):
    """UCI Metro Interstate Traffic Volume, hourly 2012–2018. Clock = date_time."""
    del seed, pca_d
    import pandas as pd

    gz = root / "data/public/Metro_Interstate_Traffic_Volume.csv.gz"
    csv = root / "data/public/Metro_Interstate_Traffic_Volume.csv"
    if gz.is_file():
        df = pd.read_csv(gz)
    elif csv.is_file():
        df = pd.read_csv(csv)
    else:
        url = "https://archive.ics.uci.edu/static/public/492/metro+interstate+traffic+volume.zip"
        df = pd.read_csv(url, compression="zip")
    df["date_time"] = pd.to_datetime(df["date_time"], errors="coerce")
    df = df.dropna(subset=["date_time", "traffic_volume"]).sort_values("date_time")
    df["hour"] = df["date_time"].dt.hour
    df["dow"] = df["date_time"].dt.dayofweek
    df["month"] = df["date_time"].dt.month
    df["is_holiday"] = (~df["holiday"].isna()) & (df["holiday"].astype(str) != "None")
    df["is_holiday"] = df["is_holiday"].astype(int)
    df["rain_1h"] = np.clip(pd.to_numeric(df["rain_1h"], errors="coerce").fillna(0), 0, 50)
    df["weather_code"] = pd.factorize(df["weather_main"].astype(str))[0]
    cols = ["temp", "rain_1h", "snow_1h", "clouds_all", "hour", "dow", "month", "is_holiday", "weather_code"]
    X = df[cols].to_numpy(np.float64)
    y = df["traffic_volume"].to_numpy(np.float64)
    n = min(int(max_n), len(X))
    return X[:n], y[:n], "mse"


def load_nyc_taxi(root: Path, max_n: int, seed: int, pca_d: int):
    """OpenML NYC green taxi Dec 2016, ordered by pickup clock. y = tip_amount."""
    del seed, pca_d
    import pandas as pd
    from sklearn.datasets import fetch_openml

    cache = root / "data/public/nyc_taxi_xy.npz"
    if cache.is_file():
        z = np.load(cache)
        n = min(int(max_n), len(z["y"]))
        return z["X"][:n], z["y"][:n], "mse"
    ds = fetch_openml(data_id=42729, as_frame=True, parser="auto")
    df = ds.data.copy()
    y = pd.to_numeric(ds.target, errors="coerce").to_numpy(np.float64)
    keep = ~np.isnan(y)
    df, y = df.loc[keep], y[keep]
    # consecutive trips: day → hour → minute (already a stream clock)
    day = pd.to_numeric(df.get("lpep_pickup_datetime_day"), errors="coerce").fillna(0)
    hour = pd.to_numeric(df.get("lpep_pickup_datetime_hour"), errors="coerce").fillna(0)
    minute = pd.to_numeric(df.get("lpep_pickup_datetime_minute"), errors="coerce").fillna(0)
    key = day.to_numpy() * 24 * 60 + hour.to_numpy() * 60 + minute.to_numpy()
    # total_amount includes the tip — drop it
    drop = [c for c in df.columns if c == "total_amount"]
    X = _as_numeric(df.drop(columns=drop, errors="ignore"))
    order = np.argsort(key, kind="mergesort")
    X, y = X[order], y[order]
    n = min(int(max_n), len(X))
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, X=X[: max(n, 12000)], y=y[: max(n, 12000)])
    return X[:n], y[:n], "mse"


def load_diabetes_readmit(root: Path, max_n: int, seed: int, pca_d: int):
    """Local source→target diabetes readmission (datasets.zip). Real shift stream."""
    del seed, pca_d
    import zipfile

    import pandas as pd

    zpath = root / "datasets/datasets.zip"
    if not zpath.is_file():
        return None
    with zipfile.ZipFile(zpath) as zf:
        src = pd.read_csv(zf.open("source_DiabetesReadmission.csv"))
        tgt = pd.read_csv(zf.open("target_DiabetesReadmission.csv"))
    ycol = "readmitted"
    df = pd.concat([src, tgt], ignore_index=True)
    y = df[ycol].to_numpy(np.int64)
    X = _as_numeric(df.drop(columns=[ycol]))
    n = min(int(max_n), len(X))
    # Keep the source→target cut inside the window (not a source-only prefix).
    n_src = len(src)
    half = n // 2
    start = max(0, n_src - half)
    end = min(len(X), start + n)
    return X[start:end], y[start:end], "acc"


def load_california():
    from sklearn.datasets import fetch_california_housing

    ds = fetch_california_housing()
    # Latitude: spatially consecutive census blocks
    return _order(ds.data, ds.target, ds.data[:, 6]) + ("mse",)


def load_diabetes():
    from sklearn.datasets import load_diabetes

    ds = load_diabetes()
    return _order(ds.data, ds.target, ds.data[:, 2]) + ("mse",)  # BMI


def load_wine_red():
    from sklearn.datasets import fetch_openml

    ds = fetch_openml("wine-quality-red", version=1, as_frame=False, parser="auto")
    X = np.asarray(ds.data, dtype=float)
    y = np.asarray(ds.target, dtype=float)
    return _order(X, y, X[:, -1]) + ("mse",)  # alcohol


def load_digits_pca(seed=0):
    from sklearn.datasets import load_digits

    ds = load_digits()
    Xp = _pca(ds.data, 16, seed)
    y = ds.target.astype(int)
    return _order(Xp, y, Xp[:, 0]) + ("acc",)


def _xy_cache_path(root, name):
    return Path(root) / "data/public" / f"{name}_xy.npz"


def _read_xy_cache(root, name, max_n):
    p = _xy_cache_path(root, name)
    if not p.is_file():
        return None
    z = np.load(p)
    n = min(int(max_n), len(z["y"]))
    return z["X"][:n], z["y"][:n]


def _write_xy_cache(root, name, X, y, cap=20000):
    p = _xy_cache_path(root, name)
    p.parent.mkdir(parents=True, exist_ok=True)
    n = min(int(cap), len(y))
    np.savez_compressed(p, X=np.asarray(X)[:n], y=np.asarray(y)[:n])


def _labels(y):
    import pandas as pd

    y = pd.Series(y).astype(str)
    codes, _ = pd.factorize(y, sort=False)
    return codes.astype(np.int64)


def load_electricity(root: Path, max_n: int, seed: int, pca_d: int):
    """OpenML 151 elec2. Row order is the 30-min clock. y = price UP/DOWN."""
    del seed, pca_d
    cached = _read_xy_cache(root, "electricity", max_n)
    if cached is not None:
        return cached[0], cached[1], "acc"
    from sklearn.datasets import fetch_openml

    ds = fetch_openml(data_id=151, as_frame=True, parser="auto")
    y = _labels(ds.target)
    X = _as_numeric(ds.data)
    _write_xy_cache(root, "electricity", X, y)
    n = min(int(max_n), len(y))
    return X[:n], y[:n], "acc"


def load_airlines(root: Path, max_n: int, seed: int, pca_d: int):
    """OpenML 1169 airlines. Row order is the flight clock. y = delay."""
    del seed, pca_d
    cached = _read_xy_cache(root, "airlines", max_n)
    if cached is not None:
        return cached[0], cached[1], "acc"
    from sklearn.datasets import fetch_openml

    ds = fetch_openml(data_id=1169, as_frame=True, parser="auto")
    y = _labels(ds.target)
    X = _as_numeric(ds.data)
    _write_xy_cache(root, "airlines", X, y)
    n = min(int(max_n), len(y))
    return X[:n], y[:n], "acc"


def _read_uci_zip_csv(url, member, **read_kw):
    import io
    import urllib.request
    import zipfile

    import pandas as pd

    raw = urllib.request.urlopen(url, timeout=90).read()
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()
        path = member if member in names else next(n for n in names if n.endswith(member))
        return pd.read_csv(zf.open(path), **read_kw)


def load_bike_hour(root: Path, max_n: int, seed: int, pca_d: int):
    """UCI Bike Sharing hour.csv, already ordered by day×hour. y = cnt."""
    del seed, pca_d
    cached = _read_xy_cache(root, "bike_hour", max_n)
    if cached is not None:
        return cached[0], cached[1], "mse"
    url = "https://archive.ics.uci.edu/static/public/275/bike+sharing+dataset.zip"
    df = _read_uci_zip_csv(url, "hour.csv")
    y = df["cnt"].to_numpy(np.float64)
    drop = [c for c in ("instant", "dteday", "casual", "registered", "cnt") if c in df.columns]
    X = _as_numeric(df.drop(columns=drop))
    _write_xy_cache(root, "bike_hour", X, y)
    n = min(int(max_n), len(y))
    return X[:n], y[:n], "mse"


def load_beijing_pm25(root: Path, max_n: int, seed: int, pca_d: int):
    """UCI Beijing PM2.5, hourly 2010–2014. Clock = year-month-day-hour."""
    del seed, pca_d
    import pandas as pd

    cached = _read_xy_cache(root, "beijing_pm25", max_n)
    if cached is not None:
        return cached[0], cached[1], "mse"
    url = (
        "https://archive.ics.uci.edu/ml/machine-learning-databases/"
        "00381/PRSA_data_2010.1.1-2014.12.31.csv"
    )
    df = pd.read_csv(url)
    df = df.dropna(subset=["pm2.5"]).copy()
    df = df.sort_values(["year", "month", "day", "hour"])
    y = df["pm2.5"].to_numpy(np.float64)
    X = _as_numeric(df.drop(columns=[c for c in ("No", "pm2.5") if c in df.columns]))
    _write_xy_cache(root, "beijing_pm25", X, y)
    n = min(int(max_n), len(y))
    return X[:n], y[:n], "mse"


def load_occupancy(root: Path, max_n: int, seed: int, pca_d: int):
    """UCI occupancy detection, concatenated and sorted by date. y = Occupancy."""
    del seed, pca_d
    import pandas as pd

    cached = _read_xy_cache(root, "occupancy", max_n)
    if cached is not None:
        return cached[0], cached[1], "acc"
    url = "https://archive.ics.uci.edu/static/public/357/occupancy+detection.zip"
    import io
    import urllib.request
    import zipfile

    raw = urllib.request.urlopen(url, timeout=90).read()
    frames = []
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        for name in zf.namelist():
            if not name.endswith(".txt"):
                continue
            frames.append(pd.read_csv(zf.open(name)))
    df = pd.concat(frames, ignore_index=True)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date", "Occupancy"]).sort_values("date")
    y = df["Occupancy"].to_numpy(np.int64)
    X = _as_numeric(df.drop(columns=["date", "Occupancy"]))
    _write_xy_cache(root, "occupancy", X, y)
    n = min(int(max_n), len(y))
    return X[:n], y[:n], "acc"


# (name, loader, clock, remote)
# remote=True loaders are skipped in CI unless a cache exists or download=True.
STREAMS = (
    ("affec", load_affec, "local", False),
    ("tencent", load_tencent, "local", False),
    ("msrvtt", load_msrvtt, "local", False),
    (
        "coco_time",
        lambda root, max_n, seed, pca_d: load_img_txt(
            root, "coco_time_order", max_n, seed, pca_d
        ),
        "local",
        False,
    ),
    (
        "fashion_iq",
        lambda root, max_n, seed, pca_d: load_img_txt(
            root, "fashion_iq", max_n, seed, pca_d
        ),
        "local",
        False,
    ),
    (
        "indiana_cxr",
        lambda root, max_n, seed, pca_d: load_img_txt(
            root, "indiana_cxr", max_n, seed, pca_d
        ),
        "local",
        False,
    ),
    ("interstate", load_interstate, "time", False),
    ("nyc_taxi", load_nyc_taxi, "time", False),
    ("electricity", load_electricity, "time", True),
    ("airlines", load_airlines, "time", True),
    ("bike_hour", load_bike_hour, "time", True),
    ("beijing_pm25", load_beijing_pm25, "time", True),
    ("occupancy", load_occupancy, "time", True),
    ("diabetes_readmit", load_diabetes_readmit, "shift", False),
    ("california", lambda root, max_n, seed, pca_d: load_california(), "spatial", False),
)

# Back-compat aliases used by older scripts/tests.
LOCAL = tuple((n, fn) for n, fn, clock, _ in STREAMS if clock == "local")
PUBLIC = tuple((n, fn) for n, fn, clock, _ in STREAMS if clock != "local")


def iter_real_streams(
    *,
    root: Path | None = None,
    n_per=200,
    n_batches=12,
    pca_d=32,
    seed=0,
    max_n=8000,
    download=False,
    clocks=None,
):
    """Yield consecutive-batch streams. Skip missing packs / too-short tables.

    ``download=False`` (CI default) will not fetch remote OpenML/UCI tables
    unless a ``data/public/{name}_xy.npz`` cache already exists.
    """
    root = Path(root or ROOT)
    allow = set(clocks) if clocks else {"time", "shift", "spatial", "local"}
    seen = set()

    def emit(name, packed, clock):
        if packed is None or name in seen:
            return None
        if len(packed) == 3:
            X, y, task = packed
        else:
            X, y = packed[:2]
            task = packed[2] if len(packed) > 2 else "mse"
        n_b = int(n_batches)
        n_p = int(n_per)
        if len(X) < n_p * n_b:
            n_b = int(len(X) // n_p)
        if n_b < 6:
            print(f"[skip] {name}: only {len(X)} rows", flush=True)
            return None
        st = _stream(
            X,
            y,
            name=name,
            task=task,
            n_per=n_p,
            n_batches=n_b,
            meta={"source": name, "seed": int(seed), "clock": clock},
        )
        if st is None:
            return None
        seen.add(name)
        return st

    for name, fn, clock, remote in STREAMS:
        if clock not in allow:
            continue
        cache = _xy_cache_path(root, name)
        if remote and not download and not cache.is_file():
            print(f"[skip] {name}: no cache (pass download=True)", flush=True)
            continue
        try:
            packed = fn(root, max_n, seed, pca_d)
        except Exception as exc:
            print(f"[skip] {name}: {exc}", flush=True)
            continue
        if packed is None:
            print(f"[skip] {name}: missing local pack", flush=True)
            continue
        st = emit(name, packed, clock)
        if st is not None:
            yield st
