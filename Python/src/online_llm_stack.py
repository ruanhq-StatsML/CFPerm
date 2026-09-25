"""Online stacking on frozen OnlineRFPerm: RF + residual LLM.

Practical kit (nothing else):
  1. Hold-out stack weight  w_stack ∈ [0, 1]  (OLS on ref hold-out; 0 = skip LLM)
  2. Local shrink            w_i = w_stack * clip(σ_med / σ_i)
  3. Inference               posterior MSE vs RF **OOB** ref MSE  → FDR p

MSE drops when the stacked residual is real. FAR drops because that lower MSE
is scored against a frozen OOB null, not re-ranked inside the trail.

The LLM's only input is the sliding-window prompt. The stack weight stays
frozen. Each trail batch renders a prompt from the past window, predicts the
residual, then appends that batch's output once Y is observed.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestRegressor

GATE = 1.5
PROMPT_WINDOW = 8


class TabPFNStyleLLM:
    def __init__(
        self,
        n_context=10,
        window=2000,
        kernel="gaussian",
        standardize=True,
        eps=1e-6,
        seed=0,
        temperature=1.0,
        auto_calibrate=False,
    ):
        self.n_context = int(n_context)
        self.window = int(window)
        self.kernel = kernel
        self.standardize = bool(standardize)
        self.eps = float(eps)
        self.temperature = max(float(temperature), self.eps)
        self.auto_calibrate = bool(auto_calibrate)
        self.X_ctx = None
        self.Y_ctx = None
        self._center = None
        self._scale = None
        self._y2 = None
        self._sigma_ref = None
        self._y_std_ref = None

    def fit_context(self, X, Y):
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float).ravel()
        if X.ndim == 1:
            X = X.reshape(1, -1)
        if X.shape[1] == 0:
            self.X_ctx = np.zeros((len(Y), 1))
            self.Y_ctx = Y[-self.window :]
            self._center = np.zeros(1)
            self._scale = np.ones(1)
            self._y2 = (self.X_ctx ** 2).sum(axis=1)
            self._sigma_ref = 1.0
            self._y_std_ref = float(np.std(Y) + self.eps)
            return self
        if self.standardize:
            self._center = X.mean(axis=0)
            self._scale = np.where(X.std(axis=0) < self.eps, 1.0, X.std(axis=0))
            X = (X - self._center) / self._scale
        self.X_ctx = X[-self.window :]
        self.Y_ctx = Y[-self.window :]
        self._y2 = (self.X_ctx ** 2).sum(axis=1)
        self._y_std_ref = float(np.std(self.Y_ctx) + self.eps)
        self._sigma_ref = self._estimate_bandwidth(self.X_ctx)
        return self

    def _estimate_bandwidth(self, X_std: np.ndarray) -> float:
        n_ctx = len(X_std)
        k = min(self.n_context, n_ctx)
        if k < 2:
            return 1.0 * self.temperature
        x2 = (X_std ** 2).sum(axis=1, keepdims=True)
        D = np.sqrt(np.maximum(x2 + x2.T - 2.0 * (X_std @ X_std.T), 0.0))
        np.fill_diagonal(D, np.inf)
        kth = np.partition(D, kth=k - 1, axis=1)[:, k - 1]
        med = float(np.median(kth[np.isfinite(kth)]))
        if not np.isfinite(med) or med <= self.eps:
            med = float(np.median(D[np.isfinite(D)]))
        return max(med * self.temperature, self.eps)

    def _kernel_weights(self, d_k: np.ndarray) -> np.ndarray:
        if self.kernel == "gaussian":
            if self.auto_calibrate and self._sigma_ref is not None:
                local = d_k.max(axis=1, keepdims=True) + self.eps
                sigma = np.maximum(self._sigma_ref, local)
            else:
                sigma = (d_k.max(axis=1, keepdims=True) + self.eps) * self.temperature
            w = np.exp(-(d_k ** 2) / (2.0 * sigma ** 2))
        else:
            w = 1.0 / (d_k * self.temperature + self.eps)
        return w / np.maximum(w.sum(axis=1, keepdims=True), self.eps)

    def predict_with_uncertainty(self, X_q):
        """Weighted kNN mean and weighted Y-std (uncertainty proxy)."""
        X_q = np.asarray(X_q, dtype=float)
        if X_q.ndim == 1:
            X_q = X_q.reshape(1, -1)
        if self.X_ctx is None or len(self.X_ctx) == 0:
            z = np.zeros(len(X_q))
            return z, np.full(len(X_q), self._y_std_ref or 1.0)
        if X_q.shape[1] == 0:
            mu = float(np.mean(self.Y_ctx))
            return np.full(len(X_q), mu), np.full(len(X_q), self._y_std_ref or 1.0)
        if self.standardize:
            X_q = (X_q - self._center) / self._scale
        n_q, n_ctx = len(X_q), len(self.X_ctx)
        k = min(self.n_context, n_ctx)
        x2 = (X_q ** 2).sum(axis=1, keepdims=True)
        D = np.sqrt(np.maximum(x2 + self._y2[None, :] - 2.0 * (X_q @ self.X_ctx.T), 0.0))
        if k < n_ctx:
            idx = np.argpartition(D, kth=k - 1, axis=1)[:, :k]
        else:
            idx = np.tile(np.arange(n_ctx), (n_q, 1))
        d_k = np.take_along_axis(D, idx, axis=1)
        y_k = self.Y_ctx[idx]
        w = self._kernel_weights(d_k)
        mu = (w * y_k).sum(axis=1)
        var = (w * (y_k - mu[:, None]) ** 2).sum(axis=1)
        sig = np.sqrt(np.maximum(var, self.eps))
        floor = (self._y_std_ref or 1.0) * 0.05
        sig = np.maximum(sig, floor)
        return mu, sig

    def predict(self, X_q):
        mu, _ = self.predict_with_uncertainty(X_q)
        return mu


def make_rf_adapter(n_estimators=80, max_depth=6):
    def fit(X, Y, seed=None, config=None):
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float).ravel()
        if X.shape[1] == 0:
            mu = float(np.mean(Y))
            return ("const", mu)
        model = RandomForestRegressor(
            n_estimators=int(n_estimators),
            max_depth=int(max_depth),
            n_jobs=1,
            random_state=int(seed or 0),
            oob_score=True,
            bootstrap=True,
        )
        model.fit(X, Y)
        return ("rf", model)

    def predict(fitted, X_new):
        kind, obj = fitted
        X_new = np.asarray(X_new, dtype=float)
        if kind == "const":
            return np.full(len(X_new), float(obj))
        return np.asarray(obj.predict(X_new), dtype=float)

    return {"fit": fit, "predict": predict, "name": "rf"}


def make_tabpfn_adapter(
    n_context=10,
    window=2000,
    seed=0,
    temperature=4.0,
    auto_calibrate=True,
):
    def fit(X, Y, seed=None, config=None):
        llm = TabPFNStyleLLM(
            n_context=n_context,
            window=window,
            seed=seed or 0,
            temperature=temperature,
            auto_calibrate=auto_calibrate,
        )
        llm.fit_context(X, Y)
        return llm

    def predict(fitted, X_new):
        return np.asarray(fitted.predict(X_new), dtype=float)

    def predict_uncertainty(fitted, X_new):
        mu, sig = fitted.predict_with_uncertainty(X_new)
        return np.asarray(mu, float), np.asarray(sig, float)

    return {
        "fit": fit,
        "predict": predict,
        "predict_uncertainty": predict_uncertainty,
        "name": "tabpfn",
    }


def compute_vimp(X, Y, n_estimators=60, seed=0):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    if X.shape[1] == 0:
        return np.zeros(0)
    rf = RandomForestRegressor(
        n_estimators=int(n_estimators),
        max_depth=6,
        n_jobs=1,
        random_state=int(seed),
    ).fit(X, Y)
    baseline = float(np.mean((Y - rf.predict(X)) ** 2))
    rng = np.random.default_rng(int(seed))
    vimp = np.zeros(X.shape[1], dtype=float)
    for j in range(X.shape[1]):
        Xp = X.copy()
        rng.shuffle(Xp[:, j])
        vimp[j] = float(np.mean((Y - rf.predict(Xp)) ** 2) - baseline)
    return vimp


def split_by_vimp(vimp, n_high=4):
    vimp = np.asarray(vimp, dtype=float)
    order = np.argsort(-vimp)
    n_high = int(min(max(int(n_high), 1), max(len(order) - 1, 1))) if len(order) else 0
    return [int(i) for i in order[:n_high]], [int(i) for i in order[n_high:]]


def hop_fires(e_now, e_prev, gate: float = GATE, e_floor: float = 0.0) -> bool:
    if e_prev is None:
        return False
    denom = max(float(e_prev), float(e_floor or 0.0), 1e-8)
    ratio = float(e_now) / denom
    return bool(np.isfinite(ratio) and ratio >= float(gate))


def _llm_uncertainty(llm_adapter, fit_llm, Xb, llm_cols):
    if "predict_uncertainty" in llm_adapter:
        return llm_adapter["predict_uncertainty"](fit_llm, Xb[:, llm_cols])
    mu = np.asarray(llm_adapter["predict"](fit_llm, Xb[:, llm_cols]), float).ravel()
    return mu, np.full(len(mu), np.nan)


def holdout_stack_weight(e_hold, y_corr_hold) -> float:
    """Layer-1 OLS: ŷ = RF + w·corr, fit on hold-out RF residuals.

    w = <e, corr> / ||corr||² clipped to [0, 1]. Negative / useless LLM → 0.
    """
    e = np.asarray(e_hold, float).ravel()
    corr = np.asarray(y_corr_hold, float).ravel()
    denom = float(np.dot(corr, corr))
    if denom <= 1e-12:
        return 0.0
    return float(np.clip(np.dot(e, corr) / denom, 0.0, 1.0))


def online_stack_weights(sig_llm, *, w_stack: float, sigma_med: float) -> np.ndarray:
    """Online application: freeze w_stack, shrink by local LLM uncertainty."""
    sig = np.maximum(np.asarray(sig_llm, float).ravel(), 1e-8)
    rel = np.clip(max(float(sigma_med), 1e-8) / sig, 0.0, 1.0)
    return np.clip(float(w_stack) * rel, 0.0, 1.0)


def posterior_weights(sig_llm, *, w_stack: float, sigma_med: float) -> np.ndarray:
    """Alias of ``online_stack_weights``."""
    return online_stack_weights(sig_llm, w_stack=w_stack, sigma_med=sigma_med)


def blend_variance_matched(
    y_dl: np.ndarray,
    y_llm: np.ndarray,
    sig_llm: np.ndarray,
    *,
    w_stack: float,
    sigma_med: float,
) -> np.ndarray:
    """Stacked posterior mean: (1-w) RF + w LLM."""
    y_dl = np.asarray(y_dl, float).ravel()
    y_llm = np.asarray(y_llm, float).ravel()
    w = online_stack_weights(sig_llm, w_stack=w_stack, sigma_med=sigma_med)
    return (1.0 - w) * y_dl + w * y_llm


def pval_vs_ref(trail: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Right-tail p of each trail score vs a frozen ref score law."""
    trail = np.asarray(trail, float).ravel()
    ref = np.asarray(ref, float).ravel()
    ref = ref[np.isfinite(ref)]
    n = max(len(ref), 1)
    if len(ref) == 0:
        return np.ones(len(trail), dtype=float)
    return np.array([(1.0 + float(np.sum(ref >= s))) / (n + 1.0) for s in trail], dtype=float)


def render_sliding_prompt(window_rows, query_x, *, w_stack: float, sigma_med: float) -> str:
    """Prompt entry. The model sees only this text.

    The window is past batches only. The query is the current x, without Y.
    The requested output is the RF residual, because the stack adds it on top
    of the frozen forest.
    """
    lines = [
        "Task: predict the residual e = Y - RF(X). Do not predict Y.",
        "The forest is frozen. w_stack and sigma_med are frozen.",
        f"w_stack={float(w_stack):.4f} sigma_med={float(sigma_med):.4f}",
        "Sliding window, oldest to newest:",
    ]
    if not window_rows:
        lines.append("(empty)")
    for row in window_rows:
        lines.append(
            "t={t} y={y:.4f} rf={rf:.4f} residual={e:.4f} w={w:.3f} x={x}".format(**row)
        )
    q = np.asarray(query_x, dtype=float).ravel()
    lines.append("query_x=" + np.array2string(q, precision=4, separator=","))
    lines.append("output: residual")
    return "\n".join(lines)


def _batch_row(t, Xb, Yb, y_rf, residual, w):
    return {
        "t": int(t),
        "y": float(np.mean(Yb)),
        "rf": float(np.mean(y_rf)),
        "e": float(np.mean(residual)),
        "w": float(np.mean(w)),
        "x": np.round(np.mean(np.asarray(Xb, float), axis=0), 4).tolist(),
    }


def frozen_llm_mse_stream(
    X,
    Y,
    *,
    llm_adapter=None,
    dl_adapter=None,
    llm_cols=None,
    dl_cols=None,
    w_llm=0.5,
    w_dl=0.5,
    ref_batch_size=1000,
    batch_size=50,
    seed=2026,
    route="prototype",
    n_high=4,
    variance_match=True,
    prompt_window=PROMPT_WINDOW,
):
    """Freeze on ref; trail uses the stacked posterior.

    Routes:
      rf         — RF only (no stack)
      prototype  — online stacking (RF + residual LLM, all columns)
      vimp       — high-VIMP → RF, low-VIMP → residual LLM, then same stack

    Each trail batch also emits a prompt built from the sliding window of
    past batch outputs, then appends this batch once Y is in hand.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n_ref = int(ref_batch_size)
    batch = int(batch_size)
    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    p_x = X_ref.shape[1]
    if dl_adapter is None:
        dl_adapter = make_rf_adapter()
    if route == "vimp" and (llm_cols is None or dl_cols is None):
        high_idx, low_idx = split_by_vimp(compute_vimp(X_ref, Y_ref, seed=seed), n_high=n_high)
        dl_cols = high_idx if dl_cols is None else dl_cols
        llm_cols = low_idx if llm_cols is None else llm_cols
    if llm_cols is None:
        llm_cols = list(range(p_x))
    if dl_cols is None:
        dl_cols = list(range(p_x))
    llm_cols = np.asarray(list(llm_cols), dtype=int)
    dl_cols = np.asarray(list(dl_cols), dtype=int)

    fit_dl = dl_adapter["fit"](X_ref[:, dl_cols], Y_ref, seed=seed)
    y_dl_ref = np.asarray(dl_adapter["predict"](fit_dl, X_ref[:, dl_cols]), float).ravel()
    y_null_ref = y_dl_ref
    kind, obj = fit_dl if isinstance(fit_dl, tuple) and len(fit_dl) == 2 else (None, None)
    if kind == "rf" and hasattr(obj, "oob_prediction_"):
        oob = np.asarray(obj.oob_prediction_, float).ravel()
        if len(oob) == len(Y_ref) and np.isfinite(oob).mean() > 0.5:
            y_null_ref = np.where(np.isfinite(oob), oob, y_dl_ref)
    sigma_rf = float(np.std(Y_ref - y_null_ref) + 1e-8)

    use_llm = route != "rf" and float(w_llm) > 0
    w_stack = 0.0
    sigma_med = sigma_rf
    if use_llm:
        if llm_adapter is None:
            llm_adapter = make_tabpfn_adapter(seed=seed)
        e_ref = Y_ref - y_dl_ref
        n_hold = max(int(0.2 * n_ref), 50)
        n_fit = max(n_ref - n_hold, 1)
        fit_llm = llm_adapter["fit"](X_ref[:n_fit, llm_cols], e_ref[:n_fit], seed=seed)
        y_corr_hold, sig_hold = _llm_uncertainty(
            llm_adapter, fit_llm, X_ref[n_fit:], llm_cols
        )
        w_stack = holdout_stack_weight(e_ref[n_fit:], y_corr_hold)
        finite = sig_hold[np.isfinite(sig_hold)]
        sigma_med = float(np.median(finite)) if len(finite) else sigma_rf
        if not np.isfinite(sigma_med) or sigma_med <= 0:
            sigma_med = sigma_rf
    else:
        fit_llm = None
        llm_adapter = llm_adapter or {}

    def _posterior(Xb):
        """RF prior + LLM residual posterior mean and per-index weight."""
        y_dl = np.asarray(dl_adapter["predict"](fit_dl, Xb[:, dl_cols]), dtype=float).ravel()
        if not use_llm:
            return y_dl, np.zeros(len(y_dl), dtype=float), y_dl
        y_corr, sig = _llm_uncertainty(llm_adapter, fit_llm, Xb, llm_cols)
        if variance_match:
            w = online_stack_weights(sig, w_stack=w_stack, sigma_med=sigma_med)
            return y_dl + w * y_corr, w, y_dl
        w = float(w_llm) / max(float(w_llm) + float(w_dl), 1e-8)
        return y_dl + w * y_corr, np.full(len(y_dl), w), y_dl

    mse_ref = float(np.mean((Y_ref - y_null_ref) ** 2))
    ref_mse_batches = []
    window = []
    for i in range(0, n_ref, batch):
        sl = slice(i, i + batch)
        if len(Y_ref[sl]) == 0:
            break
        y_hat, w, y_rf = _posterior(X_ref[sl])
        ref_mse_batches.append(float(np.mean((Y_ref[sl] - y_null_ref[sl]) ** 2)))
        window.append(_batch_row(i // batch, X_ref[sl], Y_ref[sl], y_rf, Y_ref[sl] - y_hat, w))
    ref_mse_batches = np.asarray(ref_mse_batches, dtype=float)
    window = window[-int(prompt_window) :]

    mse_list = []
    w_list = []
    hop_det = []
    prompts = []
    outputs = []
    e_prev = None
    X_trail, Y_trail = X[n_ref:], Y[n_ref:]
    for i in range(0, len(Y_trail), batch):
        sl = slice(i, i + batch)
        if len(Y_trail[sl]) == 0:
            break
        Xb, Yb = X_trail[sl], Y_trail[sl]
        query_x = np.mean(Xb, axis=0)
        prompt = render_sliding_prompt(window, query_x, w_stack=w_stack, sigma_med=sigma_med)
        y_hat, w, y_rf = _posterior(Xb)
        mse = float(np.mean((Yb - y_hat) ** 2))
        mse_list.append(mse)
        w_list.append(float(np.mean(w)))
        hop_det.append(hop_fires(mse, e_prev, e_floor=1.0 / max(len(Yb), 1)))
        e_prev = mse
        prompts.append(prompt)
        outputs.append(
            {
                "t": int(n_ref // batch + i // batch),
                "y_hat": float(np.mean(y_hat)),
                "residual_hat": float(np.mean(y_hat - y_rf)),
                "mse": mse,
                "w": float(np.mean(w)),
            }
        )
        window.append(_batch_row(outputs[-1]["t"], Xb, Yb, y_rf, Yb - y_hat, w))
        window = window[-int(prompt_window) :]

    mse_arr = np.asarray(mse_list, dtype=float)
    pval = pval_vs_ref(mse_arr, ref_mse_batches)
    return {
        "MSE_list": mse_arr,
        "T_list": mse_arr - mse_ref,
        "pval": pval,
        "hop_det": np.asarray(hop_det, dtype=bool),
        "mse_ref": mse_ref,
        "mse_trail": float(np.mean(mse_arr)) if len(mse_arr) else np.nan,
        "ref_mse_batches": ref_mse_batches,
        "w_batch": np.asarray(w_list, dtype=float),
        "n_llm_cols": int(len(llm_cols)) if use_llm else 0,
        "n_dl_cols": int(len(dl_cols)),
        "llm_cols": llm_cols.tolist() if use_llm else [],
        "dl_cols": dl_cols.tolist(),
        "w_llm": float(w_llm) if use_llm else 0.0,
        "w_dl": float(w_dl),
        "n_trail_batches": int(len(mse_arr)),
        "route": str(route),
        "sigma_rf": float(sigma_rf),
        "w_stack": float(w_stack),
        "variance_match": bool(variance_match and use_llm),
        "stacking": "holdout_w * local_σ  vs  RF_OOB_null",
        "prompts": prompts,
        "outputs": outputs,
        "prompt_window": int(prompt_window),
    }


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n, p = 400, 4
    X = rng.normal(size=(n, p))
    Y = X[:, 0] - 0.5 * X[:, 1] + rng.normal(scale=0.3, size=n)
    out = frozen_llm_mse_stream(X, Y, ref_batch_size=200, batch_size=20, seed=0, prompt_window=4)
    print("w_stack", round(out["w_stack"], 4), "batches", out["n_trail_batches"])
    print(out["prompts"][0])
    print("output", out["outputs"][0])
