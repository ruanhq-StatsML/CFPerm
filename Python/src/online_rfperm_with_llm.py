"""OnlineRFPerm with LLM / RF routing, plus the detector board.

Same contract as onlinePermOOB_wholedf: last column of df is Y, never a
feature. Fit once on the reference window. Trail batches score T = MSE − E_ref.
ADDIS is the primary mark. The rest of the board sits on the same MSE stream.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from online_fdr import addis, saffron
from online_rfperm import hop_fires


class TabPFNStyleLLM:
    def __init__(self, n_context=10, window=2000, kernel="gaussian",
                 standardize=True, eps=1e-6, seed=0):
        self.n_context = int(n_context)
        self.window = int(window)
        self.kernel = kernel
        self.standardize = bool(standardize)
        self.eps = float(eps)
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
        self.X_ctx = X[-self.window:]
        self.Y_ctx = Y[-self.window:]
        self._y2 = (self.X_ctx ** 2).sum(axis=1)
        return self

    def predict(self, X_q):
        X_q = np.asarray(X_q, dtype=float)
        if X_q.ndim == 1:
            X_q = X_q.reshape(1, -1)
        if self.standardize:
            X_q = (X_q - self._center) / self._scale
        n_q, n_ctx = len(X_q), len(self.X_ctx)
        k = min(self.n_context, n_ctx)
        x2 = (X_q ** 2).sum(axis=1, keepdims=True)
        D = np.sqrt(np.maximum(x2 + self._y2[None, :] - 2.0 * (X_q @ self.X_ctx.T), 0.0))
        if k < n_ctx:
            idx = np.argpartition(D, k, axis=1)[:, :k]
        else:
            idx = np.tile(np.arange(n_ctx), (n_q, 1))
        d_k = np.take_along_axis(D, idx, axis=1)
        y_k = self.Y_ctx[idx]
        if self.kernel == "gaussian":
            sigma = d_k[:, -1:].mean(axis=1, keepdims=True) + self.eps
            w = np.exp(-(d_k ** 2) / (2.0 * sigma ** 2))
        else:
            w = 1.0 / (d_k + self.eps)
        w = w / w.sum(axis=1, keepdims=True)
        return (w * y_k).sum(axis=1)


def make_rf_adapter(n_estimators=80, max_depth=6):
    def fit(X, Y, seed=None, config=None):
        model = RandomForestRegressor(
            n_estimators=int(n_estimators), max_depth=int(max_depth),
            n_jobs=1, random_state=int(seed or 0),
        )
        model.fit(np.asarray(X, dtype=float), np.asarray(Y, dtype=float).ravel())
        return model

    def predict(fitted, X_new):
        return np.asarray(fitted.predict(np.asarray(X_new, dtype=float)), dtype=float)

    return {"fit": fit, "predict": predict, "name": "rf"}


def make_tabpfn_adapter(n_context=10, window=2000, seed=0):
    llm = TabPFNStyleLLM(n_context=n_context, window=window, seed=seed)

    def fit(X, Y, seed=None, config=None):
        llm.fit_context(X, Y)
        return llm

    def predict(fitted, X_new):
        return np.asarray(fitted.predict(X_new), dtype=float)

    return {"fit": fit, "predict": predict, "name": "tabpfn"}


def compute_vimp(X, Y, n_estimators=60, seed=0):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    rf = RandomForestRegressor(
        n_estimators=int(n_estimators), max_depth=6, n_jobs=1, random_state=int(seed),
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
    order = np.argsort(-np.asarray(vimp, dtype=float))
    n_high = int(n_high)
    return [int(i) for i in order[:n_high]], [int(i) for i in order[n_high:]]


def first_k_consecutive_rej(det, k):
    det = np.asarray(det, dtype=bool)
    if len(det) < k:
        return -1
    for i in range(len(det) - k + 1):
        if det[i:i + k].all():
            return int(i)
    return -1


def empirical_pval(stream, burnin=5):
    """Large MSE / T → small p. Burn-in is p=1."""
    stream = np.asarray(stream, dtype=float).ravel()
    pvals = np.ones(len(stream), dtype=float)
    for i, v in enumerate(stream):
        if i < int(burnin):
            continue
        hist = stream[:i]
        pvals[i] = (1.0 + float(np.sum(hist >= v))) / (1.0 + float(len(hist)))
    return pvals


def page_hinkley(x, delta=0.005, threshold=5.0, min_instances=20):
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


def ewma_shift(x, r=0.2, burnin=8, z=3.0):
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    if len(x) == 0:
        return det
    mu0 = float(np.mean(x[:max(int(burnin), 1)]))
    s0 = float(max(np.std(x[:max(int(burnin), 1)]), 1e-8))
    zbar = mu0
    for i, v in enumerate(x):
        zbar = float(r) * float(v) + (1.0 - float(r)) * zbar
        if i >= int(burnin) and abs(zbar - mu0) > float(z) * s0:
            det[i] = True
    return det


def cusum_shift(x, drift=0.5, threshold=4.0, burnin=8):
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    if len(x) == 0:
        return det
    mu0 = float(np.mean(x[:max(int(burnin), 1)]))
    s0 = float(max(np.std(x[:max(int(burnin), 1)]), 1e-8))
    gp = 0.0
    for i, v in enumerate(x):
        gp = max(0.0, gp + (float(v) - mu0) / s0 - float(drift))
        if i >= int(burnin) and gp > float(threshold):
            det[i] = True
    return det


def ddm_on_flags(flags, warning=2.0, drift=3.0, min_n=20):
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


def adwin_lite(x, min_len=16, z=2.8):
    x = np.asarray(x, dtype=float).ravel()
    det = np.zeros(len(x), dtype=bool)
    w = []
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
            se = float(np.sqrt(arr[:cut].var() / cut + arr[cut:].var() / (n - cut) + 1e-12))
            if abs(m1 - m2) > float(z) * se:
                det[i] = True
                w = w[cut:]
                break
    return det


def martingale_reject(pvals, alpha=0.05, epsilon=0.7):
    pvals = np.asarray(pvals, dtype=float).ravel()
    det = np.zeros(len(pvals), dtype=bool)
    m = 1.0
    thr = 1.0 / max(float(alpha), 1e-12)
    for i, p in enumerate(pvals):
        p = float(np.clip(p, 1e-12, 1.0))
        m *= float(epsilon) * (p ** (float(epsilon) - 1.0))
        det[i] = bool(m >= thr)
    return det


def _pack_det(results, name, det):
    det = np.asarray(det, dtype=bool).ravel()
    results[f"{name}_SUM"] = int(det.sum())
    results[f"{name}_1"] = first_k_consecutive_rej(det, 1)
    results[f"{name}_2"] = first_k_consecutive_rej(det, 2)
    results[f"{name}_3"] = first_k_consecutive_rej(det, 3)


def onlinePermOOB_with_LLM(
    df,
    llm_adapter=None,
    dl_adapter=None,
    llm_cols=None,
    dl_cols=None,
    w_llm=0.5,
    w_dl=0.5,
    ref_batch_size=1000,
    batch_size=50,
    burnin=5,
    seed=2026,
    alpha=0.05,
):
    """Frozen blend on df[:ref], then the OnlineRFPerm board on the trail.

    df: last column is Y. llm_cols / dl_cols index the X block only.
    """
    results = {}
    df = np.asarray(df, dtype=float)
    p = df.shape[1]
    ref_batch_size = int(ref_batch_size)
    batch_size = int(batch_size)
    df_ref = df[:ref_batch_size, :]
    df_ref_X = np.asarray(df_ref[:, :(p - 1)]).astype(float)
    df_ref_Y = np.asarray(df_ref[:, (p - 1)]).astype(float)
    p_x = df_ref_X.shape[1]
    if llm_cols is None:
        llm_cols = list(range(p_x))
    if dl_cols is None:
        dl_cols = list(range(p_x))
    llm_cols = np.asarray(list(llm_cols), dtype=int)
    dl_cols = np.asarray(list(dl_cols), dtype=int)
    if llm_adapter is None:
        llm_adapter = make_tabpfn_adapter(seed=seed)
    if dl_adapter is None:
        dl_adapter = make_rf_adapter()

    fit_llm = llm_adapter["fit"](df_ref_X[:, llm_cols], df_ref_Y, seed=seed)
    fit_dl = dl_adapter["fit"](df_ref_X[:, dl_cols], df_ref_Y, seed=seed)

    def _predict(X):
        y_llm = np.asarray(llm_adapter["predict"](fit_llm, X[:, llm_cols]), dtype=float).ravel()
        y_dl = np.asarray(dl_adapter["predict"](fit_dl, X[:, dl_cols]), dtype=float).ravel()
        return float(w_llm) * y_llm + float(w_dl) * y_dl

    y_hat_ref = _predict(df_ref_X)
    mse_ref = float(np.mean((df_ref_Y - y_hat_ref) ** 2))
    results["mse_ref"] = mse_ref

    df_trail = df[ref_batch_size:, :]
    MSE_list = []
    hop_det = []
    e_prev = None
    for i in range(0, len(df_trail), batch_size):
        df_batch = df_trail[i:i + batch_size, :]
        if len(df_batch) == 0:
            break
        Xb = np.asarray(df_batch[:, :(p - 1)]).astype(float)
        Yb = np.asarray(df_batch[:, (p - 1)]).astype(float)
        y_hat = _predict(Xb)
        mse = float(np.mean((Yb - y_hat) ** 2))
        MSE_list.append(mse)
        hop_det.append(hop_fires(mse, e_prev))
        e_prev = mse

    MSE_list = np.asarray(MSE_list, dtype=float)
    T_list = MSE_list - mse_ref
    results["MSE_list"] = MSE_list
    results["T_list"] = T_list
    pvals = empirical_pval(T_list, burnin=burnin)
    p_addis = np.where(T_list > 0.0, pvals, 1.0)
    results["pvals"] = p_addis

    rej_addis = addis(p_addis, alpha=alpha)["reject"]
    rej_saffron = saffron(p_addis, alpha=alpha)["reject"]
    rej_fix = p_addis <= float(alpha)
    _pack_det(results, "addis", rej_addis)
    _pack_det(results, "saffron", rej_saffron)
    _pack_det(results, "fix", rej_fix)
    _pack_det(results, "hop", hop_det)

    scale = float(max(np.std(MSE_list[:max(len(MSE_list) // 4, 1)]), 1e-8)) if len(MSE_list) else 1e-8
    _pack_det(
        results, "PH",
        page_hinkley(MSE_list, delta=0.005 * scale, threshold=5.0 * scale,
                     min_instances=max(int(burnin) * 2, 8)),
    )
    _pack_det(results, "EWMA", ewma_shift(MSE_list, burnin=burnin))
    _pack_det(results, "CUSUM", cusum_shift(MSE_list, burnin=burnin))
    q90 = float(np.quantile(MSE_list[:max(int(burnin), 1)], 0.9)) if len(MSE_list) else 0.0
    _pack_det(results, "DDM", ddm_on_flags((MSE_list > q90).astype(float), min_n=max(int(burnin) * 2, 8)))
    _pack_det(results, "ADWIN", adwin_lite(MSE_list))
    _pack_det(results, "martingale", martingale_reject(p_addis, alpha=alpha))

    results["n_llm_cols"] = int(len(llm_cols))
    results["n_dl_cols"] = int(len(dl_cols))
    results["w_llm"] = float(w_llm)
    results["w_dl"] = float(w_dl)
    results["first_rejection_batch"] = results["addis_1"]
    return results


def make_stationary_df(n=3000, p=10, seed=0):
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(int(n), int(p)))
    w = np.zeros(int(p))
    k = min(4, int(p))
    w[:k] = np.array([2.0, -1.5, 1.2, 0.9][:k])
    Y = X @ w + rng.normal(0.0, 0.3, size=int(n))
    if p >= 6:
        Y = Y + 0.3 * X[:, 4] * X[:, 5]
    return np.column_stack([X, Y])


def make_random_noise_df(n=3000, p=10, seed=0):
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(int(n), int(p)))
    Y = rng.normal(size=int(n))
    return np.column_stack([X, Y])


DETECTORS = (
    "addis", "saffron", "fix", "hop", "PH", "EWMA", "CUSUM", "DDM", "ADWIN", "martingale",
)


def far_from_runs(runs):
    out = {}
    n = max(len(runs), 1)
    for name in DETECTORS:
        fired = sum(int(rec.get(f"{name}_1", -1) not in (-1, None)) for rec in runs)
        out[name] = float(fired) / float(n)
    return out


def run_far_study(n_reps=20, n=2400, p=10, ref_batch_size=400, batch_size=40, seed0=2026):
    llm = make_tabpfn_adapter(seed=seed0)
    dl = make_rf_adapter()
    by_kind = {}
    makers = {"stationary": make_stationary_df, "random_noise": make_random_noise_df}
    for kind, maker in makers.items():
        runs = []
        for r in range(int(n_reps)):
            seed = int(seed0) + r
            df = maker(n=n, p=p, seed=seed)
            X_ref = df[:ref_batch_size, :-1]
            Y_ref = df[:ref_batch_size, -1]
            high_idx, low_idx = split_by_vimp(compute_vimp(X_ref, Y_ref, seed=seed), n_high=4)
            rec = onlinePermOOB_with_LLM(
                df, llm_adapter=llm, dl_adapter=dl,
                llm_cols=low_idx, dl_cols=high_idx,
                w_llm=0.5, w_dl=0.5,
                ref_batch_size=ref_batch_size, batch_size=batch_size,
                burnin=5, seed=seed,
            )
            runs.append(rec)
        n_batches = (int(n) - int(ref_batch_size)) // int(batch_size)
        by_kind[kind] = {"far": far_from_runs(runs), "n_reps": int(n_reps), "n_batches": int(n_batches)}
    return by_kind


def far_table_latex(study):
    labels = {
        "addis": r"OnlineRFPerm + ADDIS",
        "saffron": r"OnlineRFPerm + SAFFRON",
        "fix": r"OnlineRFPerm + fix-$\alpha$",
        "hop": r"last-two hop ($1.5\times$)",
        "PH": r"Page--Hinkley",
        "EWMA": r"EWMA",
        "CUSUM": r"CUSUM",
        "DDM": r"DDM",
        "ADWIN": r"ADWIN-lite",
        "martingale": r"betting martingale",
    }
    rows = []
    for name in DETECTORS:
        s = 100.0 * float(study["stationary"]["far"][name])
        n = 100.0 * float(study["random_noise"]["far"][name])
        rows.append(f"{labels[name]} & {s:.1f} & {n:.1f} \\\\")
    n_reps = study["stationary"]["n_reps"]
    n_batches = study["stationary"]["n_batches"]
    body = "\n".join(rows)
    return (
        r"\begin{table}[ht]" "\n"
        r"\centering" "\n"
        rf"\caption{{False-alarm rate, no labeled onset. high-VIMP RF + low-VIMP TabPFN-style $k$NN. "
        rf"Share of {n_reps} replications with at least one rejection over {n_batches} batches. "
        r"ADDIS is the primary OnlineRFPerm mark.}" "\n"
        r"\label{tab:far}" "\n"
        r"\begin{tabular}{lcc}" "\n"
        r"\toprule" "\n"
        r"method & stationary FAR (\%) & random-noise FAR (\%) \\" "\n"
        r"\midrule" "\n"
        f"{body}\n"
        r"\bottomrule" "\n"
        r"\end{tabular}" "\n"
        r"\end{table}" "\n"
    )


if __name__ == "__main__":
    df = make_stationary_df(n=2000, p=10, seed=0)
    X_ref = df[:500, :-1]
    Y_ref = df[:500, -1]
    high_idx, low_idx = split_by_vimp(compute_vimp(X_ref, Y_ref, seed=0), n_high=4)
    rec = onlinePermOOB_with_LLM(
        df,
        llm_cols=low_idx,
        dl_cols=high_idx,
        ref_batch_size=500,
        batch_size=50,
        seed=2026,
    )
    print("addis_1", rec["addis_1"], "saffron_1", rec["saffron_1"], "PH_1", rec["PH_1"])

