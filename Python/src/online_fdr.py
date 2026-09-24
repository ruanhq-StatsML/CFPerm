"""Online FDR: SAFFRON and ADDIS on a stream of p-values.

SAFFRON (Ramdas, Zrnic, Wainwright, Jordan 2018): alpha-investing with
candidate threshold λ. p > λ are not candidates.

ADDIS (Tian & Ramdas 2019): discard p > τ, then SAFFRON-like on the
selected subsequence with candidate threshold λ ≤ τ. Better when most
nulls are conservative — the MSE-quiet batches look like that.

Default for the board: ADDIS. MSE-broken batches feed OnlineRFPerm
rank-p; quiet batches are discarded (p=1).
"""
from __future__ import annotations

import numpy as np

ALPHA = 0.05
SAFFRON_LAMBDA = 0.5
ADDIS_LAMBDA = 0.25
ADDIS_TAU = 0.5


def gamma_infinite(n: int) -> np.ndarray:
    """γ_t = 1/(t(t+1)), Σ_t γ_t = 1. Does not peek at a future horizon."""
    t = np.arange(1, max(int(n), 1) + 1, dtype=float)
    return 1.0 / (t * (t + 1.0))


def _g(gamma: np.ndarray, idx: int) -> float:
    if idx < 0:
        return float(gamma[0])
    if idx >= len(gamma):
        return float(gamma[-1])
    return float(gamma[idx])


def saffron(pvals, *, alpha: float = ALPHA, lambda_: float = SAFFRON_LAMBDA, w0: float | None = None) -> dict:
    """SAFFRON. α_t = min(λ, (1−λ) · (W0 γ_{t−τ*} + (α−W0) Σ_j γ_{t−τ_j}))."""
    p = np.asarray(pvals, dtype=float).ravel()
    n = len(p)
    w0 = float(alpha) / 2.0 if w0 is None else float(w0)
    gamma = gamma_infinite(max(n, 2))
    reject = np.zeros(n, dtype=bool)
    alph = np.zeros(n)
    cand = np.zeros(n, dtype=bool)
    last_cand = None
    for t in range(n):
        g_idx = t if last_cand is None else t - int(last_cand) - 1
        wealth = w0 * _g(gamma, g_idx)
        for rt in np.flatnonzero(reject[:t]):
            wealth += (float(alpha) - w0) * _g(gamma, t - int(rt) - 1)
        alph[t] = min(float(lambda_), (1.0 - float(lambda_)) * wealth)
        reject[t] = bool(p[t] <= alph[t])
        cand[t] = bool(p[t] <= float(lambda_))
        if cand[t]:
            last_cand = t
    return {
        "method": "saffron",
        "reject": reject,
        "alpha_t": alph,
        "candidate": cand,
        "discarded": np.zeros(n, dtype=bool),
        "p": p,
    }


def addis(
    pvals,
    *,
    alpha: float = ALPHA,
    lambda_: float = ADDIS_LAMBDA,
    tau: float = ADDIS_TAU,
    w0: float | None = None,
) -> dict:
    """ADDIS. Discard p > τ. On the rest, candidate if p ≤ λ, factor (τ−λ)/τ."""
    p = np.asarray(pvals, dtype=float).ravel()
    n = len(p)
    w0 = float(alpha) / 2.0 if w0 is None else float(w0)
    gamma = gamma_infinite(max(n, 2))
    reject = np.zeros(n, dtype=bool)
    alph = np.zeros(n)
    cand = np.zeros(n, dtype=bool)
    discarded = p > float(tau)
    selected: list[int] = []
    last_cand_k = None
    factor = (float(tau) - float(lambda_)) / max(float(tau), 1e-12)
    for t in range(n):
        if discarded[t]:
            alph[t] = 0.0
            continue
        selected.append(t)
        k = len(selected)  # 1-based position on the selected stream
        g_idx = k - 1 if last_cand_k is None else k - int(last_cand_k) - 1
        wealth = w0 * _g(gamma, g_idx)
        for rt in np.flatnonzero(reject[:t]):
            k_j = selected.index(int(rt)) + 1
            wealth += (float(alpha) - w0) * _g(gamma, k - k_j - 1)
        alph[t] = min(float(lambda_), factor * wealth)
        reject[t] = bool(p[t] <= alph[t])
        cand[t] = bool(p[t] <= float(lambda_))
        if cand[t]:
            last_cand_k = k
    return {
        "method": "addis",
        "reject": reject,
        "alpha_t": alph,
        "candidate": cand,
        "discarded": discarded,
        "p": p,
    }


def pvals_mse_gated(rows) -> np.ndarray:
    """OnlineRFPerm rank-p only when serving MSE is broken. Else p=1 (discard)."""
    out = []
    for r in rows:
        mse_b = bool(r.get("mse_broken") or r.get("mse_large"))
        if mse_b and r.get("rfperm_p") is not None:
            out.append(float(r["rfperm_p"]))
        else:
            out.append(1.0)
    return np.asarray(out, dtype=float)


def annotate_mse_rfperm_fdr(rows, *, method: str = "addis", alpha: float = ALPHA) -> dict:
    """MSE collapse → OnlineRFPerm p → SAFFRON or ADDIS. Mutates rows."""
    rows = list(rows)
    p = pvals_mse_gated(rows)
    name = str(method or "addis").lower()
    out = saffron(p, alpha=alpha) if name == "saffron" else addis(p, alpha=alpha)
    tested = [bool(r.get("mse_broken") or r.get("mse_large")) for r in rows]
    for r, pv, rej, a, cand, disc, te in zip(
        rows, out["p"], out["reject"], out["alpha_t"], out["candidate"], out["discarded"], tested
    ):
        r["fdr_method"] = out["method"]
        r["fdr_p"] = float(pv)
        r["fdr_alpha"] = float(a)
        r["fdr_reject"] = bool(rej)
        r["fdr_candidate"] = bool(cand)
        r["fdr_discarded"] = bool(disc)
        r["fdr_tested"] = bool(te)
    first = next((int(r["t"]) for r in rows if r.get("fdr_reject")), None)
    return {
        "fdr_method": out["method"],
        "fdr_level": float(alpha),
        "n_fdr_tested": int(sum(tested)),
        "n_fdr_reject": int(np.sum(out["reject"])),
        "n_fdr_candidate": int(np.sum(out["candidate"])),
        "n_fdr_discarded": int(np.sum(out["discarded"])),
        "onset_fdr": first,
    }


def both_fdr(rows, *, alpha: float = ALPHA) -> dict:
    """ADDIS primary, SAFFRON as the contrast. Restores ADDIS marks on rows."""
    saff = annotate_mse_rfperm_fdr(rows, method="saffron", alpha=alpha)
    saff_rej = [bool(r.get("fdr_reject")) for r in rows]
    add = annotate_mse_rfperm_fdr(rows, method="addis", alpha=alpha)
    for r, s in zip(rows, saff_rej):
        r["fdr_reject_saffron"] = bool(s)
        r["fdr_reject_addis"] = bool(r.get("fdr_reject"))
    return {"addis": add, "saffron": saff}
