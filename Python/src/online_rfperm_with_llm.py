"""OnlineRFPerm with flexible LLM / RF routing and a detector benchmark.

Fit the blend once on D_ref. Each new batch is T_t = MSE_t − E_ref.
ADDIS is the primary mark (same as the OnlineRFPerm board). SAFFRON,
fix-α, Page-Hinkley, EWMA, DDM, ADWIN-lite, CUSUM, and a betting
martingale sit on the same MSE stream as the benchmark suite.

Y is never a feature. llm_cols / dl_cols are integer positions in X.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional, Sequence

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from online_fdr import addis, saffron
from online_rfperm import hop_fires

BuildFn = Callable[[], Any]
FitFn = Callable[..., Any]
PredictFn = Callable[..., np.ndarray]


@dataclass(frozen=True)
class DLModelAdapter:
    name: str
    build: BuildFn
    fit: FitFn
    predict: PredictFn
    save: Optional[Callable] = None
    load: Optional[Callable] = None


# --------------------------------------------------------------------------- #
# TabPFN-style LLM (context kNN). Not a hosted LLM.
# --------------------------------------------------------------------------- #
class TabPFNStyleLLM:
    def __init__(
        self,
        n_context: int = 10,
        window: int = 2000,
        *,
        kernel: str = "gaussian",
        standardize: bool = True,
        eps: float = 1e-6,
        seed: int = 0,
    ):
        self.n_context = int(n_context)
        self.window = int(window)
        self.kernel = str(kernel)
        self.standardize = bool(standardize)
        self.eps = float(eps)
        self.rng = np.random.default_rng(int(seed))
        self.X_ctx = None
        self.Y_ctx = None
        self._center = None
        self._scale = None
        self._y2 = None

    def fit_context(self, X, Y):
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float).ravel()
        if X.ndim == 1:
            X = X.reshape(1, -1)
        if self.standardize:
            self._center = X.mean(axis=0)
            self._scale = np.where(X.std(axis=0) < self.eps, 1.0, X.std(axis=0))
            X = (X - self._center) / self._scale
        self.X_ctx = X[-self.window :]
        self.Y_ctx = Y[-self.window :]
        self._y2 = (self.X_ctx**2).sum(axis=1)
        return self

    def _dist_matrix(self, X_q):
        x2 = (X_q**2).sum(axis=1, keepdims=True)
        xy = X_q @ self.X_ctx.T
        d2 = x2 + self._y2[None, :] - 2.0 * xy
        return np.sqrt(np.maximum(d2, 0.0))

    def _weight(self, d_k):
        if self.kernel == "gaussian":
            sigma = d_k[:, -1:].mean(axis=1, keepdims=True) + self.eps
            w = np.exp(-(d_k**2) / (2.0 * sigma**2))
        elif self.kernel == "inverse":
            w = 1.0 / (d_k + self.eps)
        else:
            w = np.ones_like(d_k)
        return w / w.sum(axis=1, keepdims=True)

    def predict(self, X_q):
        X_q = np.asarray(X_q, dtype=float)
        if X_q.ndim == 1:
            X_q = X_q.reshape(1, -1)
        if self.standardize:
            X_q = (X_q - self._center) / self._scale
        n_q, n_ctx = len(X_q), len(self.X_ctx)
        k = min(self.n_context, n_ctx)
        D = self._dist_matrix(X_q)
        if k < n_ctx:
            idx = np.argpartition(D, k, axis=1)[:, :k]
        else:
            idx = np.tile(np.arange(n_ctx), (n_q, 1))
        d_k = np.take_along_axis(D, idx, axis=1)
        y_k = self.Y_ctx[idx]
        w = self._weight(d_k)
        return (w * y_k).sum(axis=1)


class TabPFNRegistry:
    def __init__(self):
        self._registry = {}

    def register(self, name, llm):
        self._registry[name] = llm
        return self

    def get(self, name):
        if name not in self._registry:
            raise KeyError(f"{name!r} is not in the TabPFN registry")
        return self._registry[name]

    def __contains__(self, name):
        return name in self._registry


def make_rf_adapter(*, n_estimators: int = 80, max_depth: int = 6) -> DLModelAdapter:
    def fit(X, Y, seed=None, config=None):
        model = RandomForestRegressor(
            n_estimators=int(n_estimators),
            max_depth=int(max_depth),
            n_jobs=1,
            random_state=int(seed or 0),
        )
        model.fit(np.asarray(X, dtype=float), np.asarray(Y, dtype=float).ravel())
        return model

    def predict(fitted, X_new):
        return np.asarray(fitted.predict(np.asarray(X_new, dtype=float)), dtype=float)

    return DLModelAdapter(name="rf", build=lambda: None, fit=fit, predict=predict)


def make_tabpfn_adapter(
    *, tabpfn_registry: TabPFNRegistry, tabpfn_name: str, n_context: int = 10
) -> DLModelAdapter:
    if tabpfn_name not in tabpfn_registry:
        raise KeyError(f"{tabpfn_name!r} is not in the TabPFN registry")

    def fit(X, Y, seed=None, config=None):
        tabpfn = tabpfn_registry.get(tabpfn_name)
        tabpfn.fit_context(np.asarray(X, dtype=float), np.asarray(Y, dtype=float).ravel())
        return {"tabpfn_name": tabpfn_name, "n_context": int(n_context)}

    def predict(fitted, X_new):
        tabpfn = tabpfn_registry.get(fitted["tabpfn_name"])
        return np.asarray(tabpfn.predict(np.asarray(X_new, dtype=float)), dtype=float)

    return DLModelAdapter(
        name=f"tabpfn_{tabpfn_name}", build=lambda: None, fit=fit, predict=predict
    )


def compute_vimp(X, Y, *, n_estimators: int = 60, seed: int = 0) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    rf = RandomForestRegressor(
        n_estimators=int(n_estimators), max_depth=6, n_jobs=1, random_state=int(seed)
    ).fit(X, Y)
    baseline = float(np.mean((Y - rf.predict(X)) ** 2))
    rng = np.random.default_rng(int(seed))
    vimp = np.zeros(X.shape[1], dtype=float)
    for j in range(X.shape[1]):
        Xp = X.copy()
        rng.shuffle(Xp[:, j])
        vimp[j] = float(np.mean((Y - rf.predict(Xp)) ** 2) - baseline)
    return vimp


def split_by_vimp(vimp, *, n_high: int = 4) -> tuple[tuple[int, ...], tuple[int, ...]]:
    order = np.argsort(-np.asarray(vimp, dtype=float))
    n_high = int(n_high)
    return tuple(int(i) for i in order[:n_high]), tuple(int(i) for i in order[n_high:])


# --------------------------------------------------------------------------- #
# Detector suite on an MSE / T stream. No river / frouros / skmultiflow.
# --------------------------------------------------------------------------- #
def first_k_consecutive_rej(det, k: int = 1) -> int | None:
    det = np.asarray(det, dtype=bool).ravel()
    k = int(k)
    if len(det) < k:
        return None
    for i in range(len(det) - k + 1):
        if det[i : i + k].all():
            return int(i)
    return None


def empirical_pval_large(stream, burnin: int = 5) -> np.ndarray:
    """p_t = P(history ≥ current). Large MSE / T → small p. Burn-in is p=1."""
    stream = np.asarray(stream, dtype=float).ravel()
    pvals = np.ones(len(stream), dtype=float)
    for i, v in enumerate(stream):
        if i < int(burnin):
            continue
        hist = stream[:i]
        pvals[i] = (1.0 + float(np.sum(hist >= v))) / (1.0 + float(len(hist)))
    return pvals


def page_hinkley(x, *, delta: float = 0.005, threshold: float = 5.0, min_instances: int = 20) -> np.ndarray:
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    mean = 0.0
    ph = 0.0
    ph_min = 0.0
    for i, v in enumerate(x):
        mean = mean + (v - mean) / float(i + 1)
        ph = ph + (v - mean - float(delta))
        ph_min = min(ph_min, ph)
        if i + 1 >= int(min_instances) and (ph - ph_min) > float(threshold):
            det[i] = True
    return det


def ewma_shift(x, *, r: float = 0.2, burnin: int = 8, z: float = 3.0) -> np.ndarray:
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    if len(x) == 0:
        return det
    mu0 = float(np.mean(x[: max(int(burnin), 1)]))
    s0 = float(max(np.std(x[: max(int(burnin), 1)]), 1e-8))
    zbar = mu0
    for i, v in enumerate(x):
        zbar = float(r) * float(v) + (1.0 - float(r)) * zbar
        if i >= int(burnin) and abs(zbar - mu0) > float(z) * s0:
            det[i] = True
    return det


def cusum_shift(x, *, drift: float = 0.5, threshold: float = 4.0, burnin: int = 8) -> np.ndarray:
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    if len(x) == 0:
        return det
    mu0 = float(np.mean(x[: max(int(burnin), 1)]))
    s0 = float(max(np.std(x[: max(int(burnin), 1)]), 1e-8))
    gp = 0.0
    for i, v in enumerate(x):
        gp = max(0.0, gp + (float(v) - mu0) / s0 - float(drift))
        if i >= int(burnin) and gp > float(threshold):
            det[i] = True
    return det


def ddm_on_flags(flags, *, warning: float = 2.0, drift: float = 3.0, min_n: int = 20) -> np.ndarray:
    flags = np.asarray(flags, dtype=float).ravel()
    det = np.zeros(len(flags), dtype=bool)
    n = 0
    p = 0.0
    p_min, s_min = 1.0, 1.0
    for i, e in enumerate(flags):
        n += 1
        p = p + (float(e) - p) / float(n)
        s = float(np.sqrt(max(p * (1.0 - p) / max(n, 1), 0.0)))
        if n >= int(min_n) and p + s < p_min + s_min:
            p_min, s_min = p, s
        if n >= int(min_n) and p + s > p_min + float(drift) * s_min:
            det[i] = True
    return det


def adwin_lite(x, *, min_len: int = 16, z: float = 2.8) -> np.ndarray:
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    w: list[float] = []
    for i, v in enumerate(x):
        w.append(float(v))
        if len(w) < int(min_len):
            continue
        arr = np.asarray(w, dtype=float)
        n = len(arr)
        for cut in (n // 3, n // 2, (2 * n) // 3):
            if cut < 6 or n - cut < 6:
                continue
            m1, m2 = float(arr[:cut].mean()), float(arr[cut:].mean())
            v1 = float(arr[:cut].var())
            v2 = float(arr[cut:].var())
            se = float(np.sqrt(v1 / cut + v2 / (n - cut) + 1e-12))
            if abs(m1 - m2) > float(z) * se:
                det[i] = True
                w = w[cut:]
                break
    return det


def martingale_reject(pvals, *, alpha: float = 0.05, epsilon: float = 0.7) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float).ravel()
    det = np.zeros(len(pvals), dtype=bool)
    m = 1.0
    thr = 1.0 / max(float(alpha), 1e-12)
    for i, p in enumerate(pvals):
        p = float(np.clip(p, 1e-12, 1.0))
        m *= float(epsilon) * (p ** (float(epsilon) - 1.0))
        det[i] = bool(m >= thr)
    return det


def run_detector_benchmark(mse, T, p_rank, *, alpha: float = 0.05, burnin: int = 5) -> dict:
    """Same MSE / T stream, every detector. first = first batch index or None."""
    mse = np.asarray(mse, dtype=float).ravel()
    T = np.asarray(T, dtype=float).ravel()
    p_rank = np.asarray(p_rank, dtype=float).ravel()
    p_emp = empirical_pval_large(T, burnin=burnin)
    p_addis = np.where(T > 0.0, p_rank, 1.0)
    scale = float(max(np.std(mse[: max(len(mse) // 4, 1)]), 1e-8))
    q90 = float(np.quantile(mse[: max(int(burnin), 1)], 0.9)) if len(mse) else 0.0
    flags = (mse > q90).astype(float)

    packs = {
        "addis": addis(p_addis, alpha=alpha)["reject"],
        "saffron": saffron(p_addis, alpha=alpha)["reject"],
        "fix_alpha": p_addis <= float(alpha),
        "page_hinkley": page_hinkley(
            mse, delta=0.005 * scale, threshold=5.0 * scale, min_instances=max(int(burnin) * 2, 8)
        ),
        "ewma": ewma_shift(mse, burnin=burnin),
        "cusum": cusum_shift(mse, burnin=burnin),
        "ddm": ddm_on_flags(flags, min_n=max(int(burnin) * 2, 8)),
        "adwin": adwin_lite(mse),
        "martingale": martingale_reject(p_addis, alpha=alpha),
        "emp_fix": p_emp <= float(alpha),
    }
    out = {"p_emp": p_emp, "p_addis": p_addis}
    for name, det in packs.items():
        det = np.asarray(det, dtype=bool).ravel()
        out[name] = {
            "reject": det,
            "n_reject": int(det.sum()),
            "first": first_k_consecutive_rej(det, 1),
            "first_3": first_k_consecutive_rej(det, 3),
        }
    hop = np.array(
        [False]
        + [hop_fires(mse[i], mse[i - 1]) for i in range(1, len(mse))],
        dtype=bool,
    )
    out["hop"] = {
        "reject": hop,
        "n_reject": int(hop.sum()),
        "first": first_k_consecutive_rej(hop, 1),
        "first_3": first_k_consecutive_rej(hop, 3),
    }
    return out


# --------------------------------------------------------------------------- #
# Main: OnlineRFPerm with LLM / RF routing + detector benchmark
# --------------------------------------------------------------------------- #
def _as_cols(cols, p: int) -> np.ndarray:
    cols = np.asarray(list(cols), dtype=int).ravel()
    if cols.size == 0:
        raise ValueError("column set is empty")
    if np.any(cols < 0) or np.any(cols >= p):
        raise ValueError(f"column index out of range for p={p}: {cols}")
    return cols


def _blend_predict(llm_adapter, fit_llm, dl_adapter, fit_dl, X, llm_cols, dl_cols, w_llm, w_dl):
    y_llm = np.asarray(llm_adapter.predict(fit_llm, X[:, llm_cols]), dtype=float).ravel()
    y_dl = np.asarray(dl_adapter.predict(fit_dl, X[:, dl_cols]), dtype=float).ravel()
    return float(w_llm) * y_llm + float(w_dl) * y_dl


def online_rfperm_with_llm(
    X,
    Y,
    *,
    llm_adapter: DLModelAdapter,
    dl_adapter: DLModelAdapter,
    llm_cols: Sequence[int],
    dl_cols: Sequence[int],
    w_llm: float = 0.5,
    w_dl: float = 0.5,
    n_ref: int = 500,
    n_new: int = 40,
    seed: int = 2026,
    alpha: float = 0.05,
    burnin: int = 5,
    config_name: str = "custom",
) -> dict:
    """Frozen blend on D_ref, then OnlineRFPerm T_t plus the detector suite."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    if X.ndim != 2 or len(Y) != len(X):
        raise ValueError("X must be 2d and aligned with Y")
    p = X.shape[1]
    llm_cols = _as_cols(llm_cols, p)
    dl_cols = _as_cols(dl_cols, p)
    n_ref = int(n_ref)
    n_new = int(n_new)
    if n_ref + n_new > len(Y):
        raise ValueError("n_ref + n_new exceeds stream length")

    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    n_cal = min(max(int(n_new) * 3, int(n_new)), max(int(n_ref) // 4, int(n_new)))
    n_fit = int(n_ref) - int(n_cal)
    if n_fit < 32:
        n_fit, n_cal = int(n_ref), 0
    X_fit, Y_fit = X_ref[:n_fit], Y_ref[:n_fit]
    fit_llm = llm_adapter.fit(X_fit[:, llm_cols], Y_fit, seed=seed)
    fit_dl = dl_adapter.fit(X_fit[:, dl_cols], Y_fit, seed=seed)

    def _mse(Xa, Ya):
        y_hat = _blend_predict(
            llm_adapter, fit_llm, dl_adapter, fit_dl, Xa, llm_cols, dl_cols, w_llm, w_dl
        )
        return float(np.mean((Ya - y_hat) ** 2))

    # E_ref and rank-p null are the blend on a held-out slice of D_ref, same
    # batch size as the stream. Do not use the RF-only FrozenRFPerm pool.
    rng = np.random.default_rng(int(seed) + 17)
    X_null = X_ref[n_fit:] if n_cal else X_ref
    Y_null = Y_ref[n_fit:] if n_cal else Y_ref
    m = min(int(n_new), len(Y_null))
    null_mse = np.empty(80, dtype=float)
    for i in range(len(null_mse)):
        sl = rng.choice(len(Y_null), size=m, replace=False)
        null_mse[i] = _mse(X_null[sl], Y_null[sl])
    mse_ref = float(np.mean(null_mse))
    null_T = null_mse - mse_ref

    def _pval(T: float) -> float:
        return float(np.sum(float(T) <= null_T) + 1.0) / float(len(null_T) + 1.0)

    mse_list = []
    T_list = []
    p_rank = []
    hops = []
    e_prev = None
    n_batches = (len(Y) - n_ref) // n_new
    for t in range(n_batches):
        lo = n_ref + t * n_new
        hi = lo + n_new
        Xb, Yb = X[lo:hi], Y[lo:hi]
        mse = _mse(Xb, Yb)
        T = float(mse - mse_ref)
        p = _pval(T)
        hop = hop_fires(mse, e_prev)
        mse_list.append(mse)
        T_list.append(T)
        p_rank.append(p)
        hops.append(hop)
        e_prev = mse

    mse_arr = np.asarray(mse_list, dtype=float)
    T_arr = np.asarray(T_list, dtype=float)
    p_arr = np.asarray(p_rank, dtype=float)
    bench = run_detector_benchmark(mse_arr, T_arr, p_arr, alpha=alpha, burnin=burnin)
    return {
        "config_name": str(config_name),
        "n_llm_cols": int(len(llm_cols)),
        "n_dl_cols": int(len(dl_cols)),
        "w_llm": float(w_llm),
        "w_dl": float(w_dl),
        "mse_ref": mse_ref,
        "MSE_list": mse_arr,
        "T_list": T_arr,
        "p_rank": p_arr,
        "n_batches": int(n_batches),
        "detectors": bench,
        "first_rejection_batch": bench["addis"]["first"],
        "mean_mse_pre": float(mse_arr[: max(len(mse_arr) // 2, 1)].mean()) if len(mse_arr) else None,
        "mean_mse_post": float(mse_arr[len(mse_arr) // 2 :].mean()) if len(mse_arr) else None,
    }


def default_adapters(*, seed: int = 2026) -> tuple[DLModelAdapter, DLModelAdapter]:
    reg = TabPFNRegistry()
    reg.register("default", TabPFNStyleLLM(n_context=10, window=2000, kernel="gaussian", seed=seed))
    llm = make_tabpfn_adapter(tabpfn_registry=reg, tabpfn_name="default", n_context=10)
    dl = make_rf_adapter(n_estimators=80, max_depth=6)
    return llm, dl


def default_routing_configs(p: int, high_idx, low_idx) -> list[dict]:
    all_cols = list(range(int(p)))
    mid = int(p) // 2
    return [
        {
            "name": "highVIMP_RF + lowVIMP_TabPFN",
            "llm_cols": list(low_idx),
            "dl_cols": list(high_idx),
            "w_llm": 0.5,
            "w_dl": 0.5,
        },
        {
            "name": "highVIMP_TabPFN + lowVIMP_RF",
            "llm_cols": list(high_idx),
            "dl_cols": list(low_idx),
            "w_llm": 0.5,
            "w_dl": 0.5,
        },
        {
            "name": "first_half_TabPFN + last_half_RF",
            "llm_cols": list(range(mid)),
            "dl_cols": list(range(mid, int(p))),
            "w_llm": 0.5,
            "w_dl": 0.5,
        },
        {"name": "all_RF", "llm_cols": all_cols, "dl_cols": all_cols, "w_llm": 0.0, "w_dl": 1.0},
        {"name": "all_TabPFN", "llm_cols": all_cols, "dl_cols": all_cols, "w_llm": 1.0, "w_dl": 0.0},
    ]


def benchmark_routing(X, Y, *, routing_configs, n_ref=500, n_new=40, seed=2026, alpha=0.05) -> list[dict]:
    llm, dl = default_adapters(seed=seed)
    out = []
    for cfg in routing_configs:
        rec = online_rfperm_with_llm(
            X,
            Y,
            llm_adapter=llm,
            dl_adapter=dl,
            llm_cols=cfg["llm_cols"],
            dl_cols=cfg["dl_cols"],
            w_llm=cfg.get("w_llm", 0.5),
            w_dl=cfg.get("w_dl", 0.5),
            n_ref=n_ref,
            n_new=n_new,
            seed=seed,
            alpha=alpha,
            config_name=cfg["name"],
        )
        out.append(rec)
    return out


# --------------------------------------------------------------------------- #
# Stationary / random-noise DGPs (no labeled onset)
# --------------------------------------------------------------------------- #
def make_stationary(n_ref: int, n_new: int, n_batches: int, p: int = 10, seed: int = 0):
    rng = np.random.default_rng(int(seed))
    n = int(n_ref) + int(n_new) * int(n_batches)
    X = rng.normal(size=(n, int(p)))
    w = np.zeros(int(p))
    k = min(4, int(p))
    w[:k] = np.array([2.0, -1.5, 1.2, 0.9][:k])
    Y = X @ w + rng.normal(0.0, 0.3, size=n)
    if p >= 6:
        Y = Y + 0.3 * X[:, 4] * X[:, 5]
    meta = {
        "kind": "stationary",
        "onset_batch": None,
        "n_ref": int(n_ref),
        "n_new": int(n_new),
        "n_batches": int(n_batches),
    }
    return X, Y, meta


def make_random_noise(n_ref: int, n_new: int, n_batches: int, p: int = 10, seed: int = 0):
    rng = np.random.default_rng(int(seed))
    n = int(n_ref) + int(n_new) * int(n_batches)
    X = rng.normal(size=(n, int(p)))
    Y = rng.normal(size=n)
    meta = {
        "kind": "random_noise",
        "onset_batch": None,
        "n_ref": int(n_ref),
        "n_new": int(n_new),
        "n_batches": int(n_batches),
    }
    return X, Y, meta


DETECTOR_ORDER = (
    "addis",
    "saffron",
    "fix_alpha",
    "hop",
    "page_hinkley",
    "ewma",
    "cusum",
    "ddm",
    "adwin",
    "martingale",
)


def far_from_runs(runs: list[dict]) -> dict[str, float]:
    """Stream FAR = share of replications with at least one rejection. No onset."""
    out = {}
    n = max(len(runs), 1)
    for name in DETECTOR_ORDER:
        fired = 0
        for rec in runs:
            first = rec["detectors"][name]["first"]
            fired += int(first is not None)
        out[name] = float(fired) / float(n)
    return out


def run_far_study(
    *,
    n_reps: int = 25,
    n_ref: int = 400,
    n_new: int = 40,
    n_batches: int = 20,
    p: int = 10,
    n_high: int = 4,
    seed0: int = 2026,
) -> dict:
    llm, dl = default_adapters(seed=seed0)
    by_kind = {}
    makers = {"stationary": make_stationary, "random_noise": make_random_noise}
    for kind, maker in makers.items():
        runs = []
        for r in range(int(n_reps)):
            seed = int(seed0) + r
            X, Y, _ = maker(n_ref, n_new, n_batches, p=p, seed=seed)
            vimp = compute_vimp(X[:n_ref], Y[:n_ref], seed=seed)
            high_idx, low_idx = split_by_vimp(vimp, n_high=n_high)
            rec = online_rfperm_with_llm(
                X,
                Y,
                llm_adapter=llm,
                dl_adapter=dl,
                llm_cols=low_idx,
                dl_cols=high_idx,
                w_llm=0.5,
                w_dl=0.5,
                n_ref=n_ref,
                n_new=n_new,
                seed=seed,
                config_name=f"{kind}/highRF+lowTabPFN",
            )
            runs.append(rec)
        by_kind[kind] = {"far": far_from_runs(runs), "n_reps": int(n_reps), "n_batches": int(n_batches)}
    return by_kind


def far_table_latex(study: dict) -> str:
    rows = []
    labels = {
        "addis": r"OnlineRFPerm + ADDIS",
        "saffron": r"OnlineRFPerm + SAFFRON",
        "fix_alpha": r"OnlineRFPerm + fix-$\alpha$",
        "hop": r"last-two hop ($1.5\times$)",
        "page_hinkley": r"Page--Hinkley",
        "ewma": r"EWMA",
        "cusum": r"CUSUM",
        "ddm": r"DDM",
        "adwin": r"ADWIN-lite",
        "martingale": r"betting martingale",
    }
    for name in DETECTOR_ORDER:
        s = 100.0 * float(study["stationary"]["far"][name])
        n = 100.0 * float(study["random_noise"]["far"][name])
        rows.append(f"{labels[name]} & {s:.1f} & {n:.1f} \\\\")
    n_reps = study["stationary"]["n_reps"]
    n_batches = study["stationary"]["n_batches"]
    body = "\n".join(rows)
    return rf"""\begin{{table}}[ht]
\centering
\caption{{False-alarm rate on streams with no labeled onset. Routing is high-VIMP RF + low-VIMP TabPFN-style kNN. FAR is the share of {n_reps} replications that reject at least once over {n_batches} batches. ADDIS is the primary OnlineRFPerm mark.}}
\label{{tab:far}}
\begin{{tabular}}{{lcc}}
\toprule
method & stationary FAR (\%) & random-noise FAR (\%) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def main() -> int:
    n_ref, n_new, n_batches, p = 400, 40, 16, 10
    X, Y, _ = make_stationary(n_ref, n_new, n_batches, p=p, seed=0)
    vimp = compute_vimp(X[:n_ref], Y[:n_ref], seed=0)
    high_idx, low_idx = split_by_vimp(vimp, n_high=4)
    configs = default_routing_configs(p, high_idx, low_idx)
    routed = benchmark_routing(X, Y, routing_configs=configs, n_ref=n_ref, n_new=n_new, seed=2026)
    print("routing (stationary, one seed)")
    print(f"{'config':<40} {'1st ADDIS':>10} {'n_ADDIS':>8}")
    for rec in routed:
        print(
            f"{rec['config_name']:<40} {str(rec['detectors']['addis']['first']):>10} "
            f"{rec['detectors']['addis']['n_reject']:>8}"
        )
    study = run_far_study(n_reps=20, n_ref=n_ref, n_new=n_new, n_batches=20, p=p, seed0=2026)
    print(far_table_latex(study))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
