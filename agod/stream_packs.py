"""Load real streaming packs (Metro Interstate, stocks, Beijing PM2.5, Waymo proxy).

All loaders return (X, y, meta) in **time order** for gradual-shift evaluation.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

import numpy as np
import pandas as pd

ArrayPack = Tuple[np.ndarray, np.ndarray, Dict]


def _onehot(s: pd.Series, prefix: str) -> pd.DataFrame:
    return pd.get_dummies(s.astype(str).fillna("NA"), prefix=prefix, dtype=np.float32)


def load_metro_interstate(root: Path, max_n: int = 20000) -> ArrayPack:
    folder = root / "data/stream_packs/metro_interstate"
    matches = sorted(folder.glob("*.csv.gz")) + sorted(folder.glob("*.csv"))
    if not matches:
        raise FileNotFoundError(f"metro_interstate csv missing under {folder}")
    p = matches[0]
    df = pd.read_csv(p)
    df["date_time"] = pd.to_datetime(df["date_time"])
    df = df.sort_values("date_time").reset_index(drop=True)
    df["hour"] = df["date_time"].dt.hour
    df["dow"] = df["date_time"].dt.dayofweek
    df["month"] = df["date_time"].dt.month
    df["is_holiday"] = (
        ~df["holiday"].isna() & (df["holiday"].astype(str) != "None")
    ).astype(float)
    base = df[
        ["temp", "rain_1h", "snow_1h", "clouds_all", "hour", "dow", "month", "is_holiday"]
    ].astype(np.float32)
    oh = _onehot(df["weather_main"], "wm")
    X = pd.concat([base, oh], axis=1).fillna(0).to_numpy(np.float32)
    y = df["traffic_volume"].to_numpy(np.float64)
    n = min(max_n, len(X))
    meta = {
        "name": "metro_interstate",
        "n": n,
        "d": int(X.shape[1]),
        "target": "traffic_volume",
        "freq": "hourly",
        "proxy": False,
    }
    return X[:n], y[:n], meta


def load_beijing_pm25(root: Path, max_n: int = 20000) -> ArrayPack:
    folder = root / "data/stream_packs/beijing_pm25"
    csvs = sorted(folder.glob("*.csv")) if folder.is_dir() else []
    if not csvs:
        raise FileNotFoundError(f"beijing_pm25 csv missing under {folder}")
    p = csvs[0]
    df = pd.read_csv(p)
    ycol = "pm2.5" if "pm2.5" in df.columns else "pm25"
    df = df.dropna(subset=[ycol]).reset_index(drop=True)
    df["stamp"] = pd.to_datetime(df[["year", "month", "day", "hour"]])
    df = df.sort_values("stamp").reset_index(drop=True)
    base = df[["DEWP", "TEMP", "PRES", "Iws", "Is", "Ir", "hour", "month"]].astype(np.float32)
    lag = df[ycol].shift(1).fillna(df[ycol].median())
    base = base.copy()
    base["pm_lag1"] = lag.astype(np.float32)
    oh = _onehot(df["cbwd"], "wind")
    X = pd.concat([base, oh], axis=1).fillna(0).to_numpy(np.float32)
    y = df[ycol].to_numpy(np.float64)
    n = min(max_n, len(X))
    meta = {
        "name": "beijing_pm25",
        "n": n,
        "d": int(X.shape[1]),
        "target": ycol,
        "freq": "hourly",
        "proxy": False,
    }
    return X[:n], y[:n], meta


def load_stocks(root: Path, ticker: str = "SPY", max_n: int = 20000, horizon: int = 1) -> ArrayPack:
    p = root / "data/stream_packs/stocks" / f"{ticker}.csv"
    df = pd.read_csv(p)
    cols = {c.lower().replace(" ", ""): c for c in df.columns}
    if "timestamp" in cols:
        raw = df[cols["timestamp"]]
        dt = pd.to_datetime(raw, unit="s", errors="coerce")
        if dt.isna().all():
            dt = pd.to_datetime(raw, errors="coerce")
        df["date"] = dt
    elif "date" in cols:
        df["date"] = pd.to_datetime(df[cols["date"]])
    else:
        df["date"] = pd.to_datetime(df.iloc[:, 0], errors="coerce")
    close_c = cols.get("adjclose", cols.get("close"))
    vol_c = cols.get("volume")
    high_c = cols.get("high", close_c)
    low_c = cols.get("low", close_c)
    df = df.assign(
        close=df[close_c].astype(float),
        volume=df[vol_c].astype(float) if vol_c else 1.0,
        high=df[high_c].astype(float),
        low=df[low_c].astype(float),
    ).sort_values("date").reset_index(drop=True)
    ret = df["close"].pct_change().fillna(0.0)
    feats = []
    for lag in (1, 2, 3, 5, 10):
        feats.append(ret.shift(lag).fillna(0.0).rename(f"ret_l{lag}"))
    feats.append(ret.rolling(5).std().fillna(0.0).rename("vol5"))
    feats.append(ret.rolling(20).std().fillna(0.0).rename("vol20"))
    feats.append(((df["high"] - df["low"]) / (df["close"] + 1e-8)).rename("hl_range"))
    feats.append(np.log1p(df["volume"]).rename("log_vol"))
    feats.append(df["date"].dt.dayofweek.astype(float).rename("dow"))
    feats.append(df["date"].dt.month.astype(float).rename("month"))
    X = pd.concat(feats, axis=1).fillna(0).to_numpy(np.float32)
    y = ret.shift(-horizon).fillna(0.0).to_numpy(np.float64)
    n = min(max_n, len(X) - horizon)
    meta = {
        "name": f"stocks_{ticker}",
        "n": n,
        "d": int(X.shape[1]),
        "target": f"ret_t+{horizon}",
        "freq": "daily",
        "proxy": False,
    }
    return X[:n], y[:n], meta


def ensure_waymo_proxy(root: Path, n: int = 12000, seed: int = 0) -> Path:
    """Kinematics proxy if real Waymo dump is absent."""
    out = root / "data/stream_packs/waymo_proxy/waymo_proxy_xy.npz"
    if out.is_file():
        import numpy as _np
        try:
            z = _np.load(out)
            if len(z["X"]) >= n:
                return out
        except Exception:
            pass
    out.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    regime = t / max(n - 1, 1)
    speed = 8.0 + 22.0 * regime + rng.normal(0, 1.5, n)
    yaw_rate = (0.25 - 0.22 * regime) * np.sin(0.03 * t) + rng.normal(0, 0.05, n)
    accel = np.gradient(speed) + rng.normal(0, 0.3, n)
    lat_acc = speed * yaw_rate + rng.normal(0, 0.2, n)
    density = 0.7 - 0.45 * regime + rng.normal(0, 0.05, n)
    lead_gap = 8.0 + 35.0 * regime + rng.normal(0, 3.0, n)
    dt = 0.1
    y = speed * dt + 0.5 * accel * dt**2 + 0.15 * lat_acc * dt + rng.normal(0, 0.08, n)
    y = y * (1.0 + 0.25 * regime)
    X = np.column_stack(
        [
            speed,
            yaw_rate,
            accel,
            lat_acc,
            density,
            lead_gap,
            np.sin(2 * np.pi * (t % 100) / 100),
            np.cos(2 * np.pi * (t % 100) / 100),
            regime,
        ]
    ).astype(np.float32)
    np.savez_compressed(out, X=X, y=y.astype(np.float64), note="waymo_kinematics_proxy_gradual")
    (out.parent / "README.md").write_text(
        "# Waymo proxy\n\n"
        "Real Waymo Open Motion files were not present. "
        "This NPZ is a **kinematics proxy** with gradual urban→highway drift "
        "for streaming PO-OOD vs DRE tests.\n",
        encoding="utf-8",
    )
    return out


def load_waymo_proxy(root: Path, max_n: int = 12000, seed: int = 0) -> ArrayPack:
    p = ensure_waymo_proxy(root, n=max_n, seed=seed)
    z = np.load(p)
    X, y = z["X"].astype(np.float32), z["y"].astype(np.float64)
    n = min(max_n, len(X))
    meta = {
        "name": "waymo_proxy",
        "n": n,
        "d": int(X.shape[1]),
        "target": "next_disp",
        "freq": "0.1s",
        "proxy": True,
    }
    return X[:n], y[:n], meta


def load_affec(root: Path, max_n: int = 10000) -> Optional[ArrayPack]:
    candidates = [
        root / "results/affec_fsds/affec_fsds_xyw_cache.npz",
        root / "results/affec_fsds/affec_fsds_xyw_cache.npz",
    ]
    cache = next((c for c in candidates if c.is_file()), None)
    if cache is None:
        return None
    z = np.load(cache, allow_pickle=True)
    ykey = "Y" if "Y" in z.files else "y"
    X, y = z["X"].astype(np.float32), z[ykey].astype(np.float64)
    n = min(max_n, len(X))
    meta = {
        "name": "affec",
        "n": n,
        "d": int(X.shape[1]),
        "target": "affect",
        "freq": "block",
        "proxy": False,
    }
    return X[:n], y[:n], meta


def _stocks(ticker: str):
    def _fn(root: Path, max_n: int = 20000) -> ArrayPack:
        return load_stocks(root, ticker=ticker, max_n=max_n)

    return _fn


LOADERS: Dict[str, Callable[..., Optional[ArrayPack]]] = {
    "metro_interstate": load_metro_interstate,
    "beijing_pm25": load_beijing_pm25,
    "stocks_SPY": _stocks("SPY"),
    "stocks_QQQ": _stocks("QQQ"),
    "stocks_AAPL": _stocks("AAPL"),
    "stocks_MSFT": _stocks("MSFT"),
    "stocks_IWM": _stocks("IWM"),
    "waymo_proxy": load_waymo_proxy,
    "affec": load_affec,
}
