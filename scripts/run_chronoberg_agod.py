#!/usr/bin/env python3
"""ChronoBerg prototype for Attribution-Guided Online Distillation (AGOD).

ChronoBerg is a temporally ordered literary corpus (1750–2000). This script
chunks the public *test* splits into passages, builds three native feature
blocks (text hashing + period-calibrated valence / arousal), and runs the
online AGOD loop from the method note:

  MSG_m = Normalize(AUC_m * VIMP_m + gamma * PO-risk_m)
  alpha  = Softmax(MSG / tau)
  L      = L_global + sum_m alpha_m L_m + lambda * Omega(alpha)

PO-risk uses dummy labels ``Y = np.arange(n)`` on the stacked (D0, Dt) batch
as requested. Teacher is a frozen random projection; the student is a linear
map whose per-modality gradient is routed by alpha.

  python3 scripts/run_chronoberg_agod.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "chronoberg"
RAW = DATA / "raw"
LEX = DATA / "lexicons"
OUT = ROOT / "results" / "chronoberg_agod"
DOCS = ROOT / "docs" / "method"
NPZ = DATA / "chronoberg_agod_windows.npz"

HF_BASE = "https://huggingface.co/datasets/spaul25/Chronoberg/resolve/main"
ERAS = (1750, 1800, 1850, 1900, 1950)
MODALITIES = ("text", "valence", "arousal")
SEED = 2026
N_PER = 180
N_WORDS = 80
D_MOD = 24
D_OUT = 24
TOP_K = 1
GAMMA = 1.0
TAU = 0.35
LAMBDA_SMOOTH = 0.15
THETA_LOW = 0.55
THETA_HIGH = 0.80
OVERLAP_FALLBACK = 0.30
ALPHA_CLIP = 0.01
N_RF = 80
N_PO = 40
N_SPLITS = 3
STEP_BUDGET = 9  # routed head-updates per window (capacity)
GLOBAL_STEPS = 1
LR = 0.25
INIT_NOISE = 0.90
INJECT_ALPHA = 0.85
MOD_SEED = {"text": 11, "valence": 23, "arousal": 37}

WORD_RE = re.compile(r"[a-zA-Z']+")


def _rng(seed: int = SEED) -> np.random.Generator:
    return np.random.default_rng(seed)


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    import urllib.request

    print(f"download {url} -> {dest}", flush=True)
    urllib.request.urlretrieve(url, dest)


def ensure_raw() -> None:
    for year in ERAS:
        _download(
            f"{HF_BASE}/train_test_split/test/{year}_test.json",
            RAW / f"{year}_test.json",
        )
        for kind in ("Valence", "Arousal", "Dominance"):
            _download(
                f"{HF_BASE}/Lexicons/{kind}_lexicon_{year}.csv",
                LEX / f"{kind}_lexicon_{year}.csv",
            )


def load_lexicon(year: int, kind: str) -> dict[str, float]:
    path = LEX / f"{kind}_lexicon_{year}.csv"
    df = pd.read_csv(path)
    word_col = "words" if "words" in df.columns else df.columns[1]
    score_col = [c for c in df.columns if c.lower() == kind.lower()]
    score_col = score_col[0] if score_col else df.columns[2]
    out = {}
    for w, s in zip(df[word_col], df[score_col]):
        if w is None or (isinstance(w, float) and np.isnan(w)):
            continue
        try:
            out[str(w).lower()] = float(s)
        except (TypeError, ValueError):
            continue
    return out


def passages_from_books(year: int) -> list[str]:
    with (RAW / f"{year}_test.json").open() as f:
        books = json.load(f)
    chunks: list[str] = []
    for book in books:
        words = WORD_RE.findall(book["text"].lower())
        if len(words) < N_WORDS:
            continue
        for i in range(0, len(words) - N_WORDS + 1, N_WORDS):
            chunks.append(" ".join(words[i : i + N_WORDS]))
            if len(chunks) >= 8000:
                return chunks
    return chunks


def sample_passages(chunks: list[str], n: int, seed: int) -> list[str]:
    rng = _rng(seed)
    if len(chunks) <= n:
        return list(chunks)
    idx = rng.choice(len(chunks), size=n, replace=False)
    return [chunks[i] for i in idx]


def hash_text(texts: list[str], n_features: int = D_MOD) -> np.ndarray:
    vec = HashingVectorizer(
        n_features=n_features,
        alternate_sign=True,
        norm="l2",
        ngram_range=(1, 2),
        lowercase=True,
    )
    return np.asarray(vec.transform(texts).toarray(), dtype=np.float64)


def _stable_bucket(token: str, n_features: int) -> int:
    digest = hashlib.md5(token.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little") % n_features


def lex_hash(texts: list[str], lex: dict[str, float], n_features: int = D_MOD) -> np.ndarray:
    X = np.zeros((len(texts), n_features), dtype=np.float64)
    for i, text in enumerate(texts):
        for tok in text.split():
            score = lex.get(tok)
            if score is None:
                continue
            X[i, _stable_bucket(tok, n_features)] += score
        nrm = np.linalg.norm(X[i])
        if nrm > 0:
            X[i] /= nrm
    return X


def build_window_features() -> dict:
    ensure_raw()
    blocks = {}
    passages = {}
    for year in ERAS:
        raw_chunks = passages_from_books(year)
        texts = sample_passages(raw_chunks, N_PER, seed=SEED + year)
        valence = load_lexicon(year, "Valence")
        arousal = load_lexicon(year, "Arousal")
        Xt = hash_text(texts)
        Xv = lex_hash(texts, valence)
        Xa = lex_hash(texts, arousal)
        X = np.hstack([Xt, Xv, Xa])
        blocks[str(year)] = X
        passages[str(year)] = np.array(texts)
        print(
            f"era {year}: n_passages={len(texts)} from {len(raw_chunks)} chunks, "
            f"X={X.shape}",
            flush=True,
        )
    return {
        "eras": np.array(ERAS, dtype=int),
        "modalities": np.array(MODALITIES),
        "d_mod": np.array([D_MOD], dtype=int),
        **{f"X_{year}": blocks[str(year)] for year in ERAS},
    }


def save_npz(payload: dict) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(NPZ, **payload)
    print(f"wrote {NPZ} ({NPZ.stat().st_size} bytes)", flush=True)


def load_windows() -> dict[int, np.ndarray]:
    if not NPZ.exists():
        payload = build_window_features()
        save_npz(payload)
    else:
        payload = dict(np.load(NPZ, allow_pickle=True))
    return {int(year): np.asarray(payload[f"X_{year}"], dtype=np.float64) for year in ERAS}


def modality_slices(d_mod: int = D_MOD) -> dict[str, tuple[int, int]]:
    return {
        "text": (0, d_mod),
        "valence": (d_mod, 2 * d_mod),
        "arousal": (2 * d_mod, 3 * d_mod),
    }


def slice_mod(X: np.ndarray, name: str) -> np.ndarray:
    a, b = modality_slices()[name]
    return X[:, a:b]


def rf_domain_auc_vimp(X0: np.ndarray, X1: np.ndarray, *, seed: int = SEED):
    X = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    p = X.shape[1]
    depth = max(2, int(round(np.sqrt(p))))
    leaf = max(1, int(round(np.sqrt(len(X)) // 2)))
    clf = RandomForestClassifier(
        n_estimators=N_RF,
        max_depth=depth,
        min_samples_leaf=leaf,
        oob_score=True,
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(np.float64)
    try:
        Xtr, Xte, Wtr, Wte = train_test_split(
            X, W, test_size=0.25, random_state=seed, stratify=W
        )
        clf2 = RandomForestClassifier(
            n_estimators=max(40, N_RF // 2),
            max_depth=depth,
            min_samples_leaf=leaf,
            n_jobs=-1,
            random_state=seed + 1,
        )
        clf2.fit(Xtr, Wtr)
        auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    except ValueError:
        auc = 0.5
    oob = float(getattr(clf, "oob_score_", 0.5))
    return vimp, auc, oob


def po_risk_block(X0: np.ndarray, X1: np.ndarray, *, seed: int = SEED) -> tuple[np.ndarray, float]:
    """CFPerm-style PO-risk with Y = np.arange(n) on the stacked window."""
    X = np.vstack([X0, X1])
    n = X.shape[0]
    Y = np.arange(n, dtype=np.float64)  # requested dummy labels
    Y = (Y - Y.mean()) / (Y.std() + 1e-8)
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=N_PO,
            max_depth=8,
            min_samples_leaf=3,
            random_state=seed + fold,
            n_jobs=-1,
        )
        e = RandomForestClassifier(
            n_estimators=N_PO,
            max_depth=8,
            min_samples_leaf=3,
            random_state=seed + 40 + fold,
            n_jobs=-1,
        )
        m.fit(X[tr], Y[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1.0 - ALPHA_CLIP)
    po = (Y - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=N_PO,
        max_depth=8,
        min_samples_leaf=3,
        random_state=seed + 7,
        n_jobs=-1,
    )
    tau.fit(X, po)
    vimp = tau.feature_importances_.astype(np.float64)
    risk = float(np.mean(tau.predict(X) ** 2))
    return vimp, risk


def _safe_div(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    s = float(np.sum(np.abs(v)))
    if s <= 1e-12:
        return np.ones_like(v) / len(v)
    return v / s


def softmax(z: np.ndarray, tau: float = TAU) -> np.ndarray:
    z = np.asarray(z, dtype=np.float64) / max(tau, 1e-8)
    z = z - np.max(z)
    e = np.exp(z)
    return e / e.sum()


def allocate_steps(alpha: np.ndarray, budget: int = STEP_BUDGET) -> np.ndarray:
    """Largest-remainder allocation so uniform alpha gets equal steps."""
    p = np.clip(np.asarray(alpha, dtype=np.float64), 1e-8, None)
    p = p / p.sum()
    raw = p * budget
    counts = np.floor(raw).astype(int)
    leftover = int(budget - counts.sum())
    order = np.argsort(-(raw - counts), kind="stable")
    for i in order[:leftover]:
        counts[i] += 1
    return counts


def msg_from_components(
    auc: dict[str, float],
    vimp: dict[str, float],
    po: dict[str, float],
    *,
    gamma: float = GAMMA,
) -> np.ndarray:
    raw = np.array(
        [auc[m] * vimp[m] + gamma * po[m] for m in MODALITIES],
        dtype=np.float64,
    )
    return _safe_div(np.maximum(raw, 0.0))


def rank_overlap(a: np.ndarray, b: np.ndarray, k: int = 8) -> float:
    sa = set(np.argsort(-np.asarray(a))[:k].tolist())
    sb = set(np.argsort(-np.asarray(b))[:k].tolist())
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / max(1, len(sa | sb))


def l2_normalize(Z: np.ndarray) -> np.ndarray:
    nrm = np.linalg.norm(Z, axis=1, keepdims=True)
    nrm = np.maximum(nrm, 1e-8)
    return Z / nrm


def recall_at_k(query: np.ndarray, gallery: np.ndarray, k: int = TOP_K) -> float:
    q = l2_normalize(query)
    g = l2_normalize(gallery)
    sim = q @ g.T
    n = len(q)
    k = min(k, n)
    hits = 0
    for i in range(n):
        top = np.argpartition(-sim[i], kth=k - 1)[:k]
        if i in top:
            hits += 1
    return hits / n


def mse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean((a - b) ** 2))


class RoutedHeads:
    """Per-modality linear student heads. Distillation *capacity* is a discrete
    step budget allocated by alpha (plus a tiny unweighted global trickle)."""

    def __init__(self, d_mod: int, d_out: int, rng: np.random.Generator):
        self.d_mod = d_mod
        self.Wt = {m: rng.normal(0.0, 0.08, size=(d_mod, d_out)) for m in MODALITIES}
        self.Ws = {
            m: self.Wt[m] + rng.normal(0.0, INIT_NOISE, size=(d_mod, d_out))
            for m in MODALITIES
        }

    def embed_mod(self, X: np.ndarray, m: str, teacher: bool = False) -> np.ndarray:
        Xm = slice_mod(X, m)
        W = self.Wt[m] if teacher else self.Ws[m]
        return Xm @ W

    def embed(self, X: np.ndarray, teacher: bool = False) -> np.ndarray:
        return np.hstack([self.embed_mod(X, m, teacher=teacher) for m in MODALITIES])

    def _step_mod(self, X: np.ndarray, m: str, lr: float) -> float:
        Xm = slice_mod(X, m)
        err = Xm @ self.Ws[m] - Xm @ self.Wt[m]
        grad = (Xm.T @ err) / max(len(Xm), 1) + 1e-4 * self.Ws[m]
        self.Ws[m] = self.Ws[m] - lr * grad
        return mse(Xm @ self.Ws[m], Xm @ self.Wt[m])

    def losses(self, X: np.ndarray) -> dict[str, float]:
        out = {m: mse(self.embed_mod(X, m, False), self.embed_mod(X, m, True)) for m in MODALITIES}
        out["global"] = mse(self.embed(X, False), self.embed(X, True))
        return out

    def distill(self, X: np.ndarray, alpha: np.ndarray, rng: np.random.Generator) -> dict[str, float]:
        for _ in range(GLOBAL_STEPS):
            for m in MODALITIES:
                self._step_mod(X, m, lr=0.25 * LR)
        counts = allocate_steps(alpha, STEP_BUDGET)
        for i, m in enumerate(MODALITIES):
            for _ in range(int(counts[i])):
                self._step_mod(X, m, lr=LR)
        out = self.losses(X)
        out["routed"] = float(np.sum(alpha * np.array([out[m] for m in MODALITIES])))
        out["step_counts"] = counts.tolist()
        return out


def compute_msg(X0: np.ndarray, X1: np.ndarray, *, seed: int = SEED) -> dict:
    """Per-modality RF-domain AUC / VIMP + PO-risk, then MSG state vector."""
    auc, vimp_share, po_risk, po_vimp = {}, {}, {}, {}
    joint_vimp, joint_auc, _ = rf_domain_auc_vimp(X0, X1, seed=seed)
    slices = modality_slices()
    joint_mass = {m: float(joint_vimp[a:b].sum()) for m, (a, b) in slices.items()}
    joint_mass = {m: v / (sum(joint_mass.values()) + 1e-12) for m, v in joint_mass.items()}
    for m in MODALITIES:
        Xm0, Xm1 = slice_mod(X0, m), slice_mod(X1, m)
        vimp_m, auc_m, _oob = rf_domain_auc_vimp(Xm0, Xm1, seed=seed + MOD_SEED[m])
        _, risk_m = po_risk_block(Xm0, Xm1, seed=seed + 17 + MOD_SEED[m])
        auc[m] = float(auc_m)
        # OOB-VIMP scalar: mass of this block inside the joint domain RF
        vimp_share[m] = float(joint_mass[m])
        po_risk[m] = float(risk_m)
        po_vimp[m] = float(vimp_m.sum())  # normalized RF always ~1; kept for logs
    po_norm = _safe_div(np.array([po_risk[m] for m in MODALITIES]))
    po_share = {m: float(po_norm[i]) for i, m in enumerate(MODALITIES)}
    g = msg_from_components(auc, vimp_share, po_share, gamma=GAMMA)
    overlap = rank_overlap(
        np.array([vimp_share[m] for m in MODALITIES]),
        np.array([po_share[m] for m in MODALITIES]),
        k=2,
    )
    used_fallback = overlap < OVERLAP_FALLBACK
    if used_fallback:
        g = _safe_div(np.array([auc[m] * vimp_share[m] for m in MODALITIES]))
    return {
        "auc": auc,
        "vimp": vimp_share,
        "po_risk": po_risk,
        "po_share": po_share,
        "joint_auc": float(joint_auc),
        "g": g.tolist(),
        "overlap": float(overlap),
        "fallback_rf_only": bool(used_fallback),
    }


def gate_alpha(alpha: np.ndarray, auc: dict[str, float], alpha_prev: np.ndarray | None) -> np.ndarray:
    a = alpha.copy()
    for i, m in enumerate(MODALITIES):
        if auc[m] < THETA_LOW:
            a[i] *= 0.35
        elif auc[m] > THETA_HIGH:
            a[i] *= 1.15
    a = _safe_div(np.maximum(a, 1e-8))
    if alpha_prev is not None:
        a = (1.0 - LAMBDA_SMOOTH) * a + LAMBDA_SMOOTH * alpha_prev
        a = _safe_div(a)
    return a


def inject_valence(windows: dict[int, np.ndarray], rng: np.random.Generator) -> dict[int, np.ndarray]:
    out = {k: v.copy() for k, v in windows.items()}
    a, b = modality_slices()["valence"]
    u = rng.normal(0.0, 1.0, size=(b - a,))
    u = u / (np.linalg.norm(u) + 1e-8)
    for year in (1900, 1950):
        out[year][:, a:b] += INJECT_ALPHA * u
        nrm = np.linalg.norm(out[year][:, a:b], axis=1, keepdims=True)
        out[year][:, a:b] /= np.maximum(nrm, 1e-8)
    return out


def run_protocol(windows: dict[int, np.ndarray], *, name: str, seed: int = SEED) -> dict:
    rng = _rng(seed)
    students = {
        "B1_static": RoutedHeads(D_MOD, D_OUT, rng),
        "B2_covariate": RoutedHeads(D_MOD, D_OUT, rng),
        "B3_agod": RoutedHeads(D_MOD, D_OUT, rng),
    }
    # share the same frozen teacher so baselines are comparable
    Wt = {m: students["B1_static"].Wt[m].copy() for m in MODALITIES}
    for s in students.values():
        s.Wt = {m: Wt[m].copy() for m in MODALITIES}
        s.Ws = {
            m: Wt[m] + rng.normal(0.0, INIT_NOISE, size=Wt[m].shape) for m in MODALITIES
        }

    X_ref = windows[ERAS[0]]
    alpha_prev = None
    rows = []
    t0 = time.perf_counter()
    for year in ERAS[1:]:
        X_t = windows[year]
        msg = compute_msg(X_ref, X_t, seed=seed + year)
        g = np.asarray(msg["g"], dtype=np.float64)
        auc_vec = np.array([msg["auc"][m] for m in MODALITIES])
        alpha_b1 = np.ones(3) / 3.0
        alpha_b2 = softmax(auc_vec, tau=TAU)
        alpha_b3 = softmax(g, tau=TAU)
        alpha_b3 = gate_alpha(alpha_b3, msg["auc"], alpha_prev)
        alpha_prev = alpha_b3.copy()
        alphas = {"B1_static": alpha_b1, "B2_covariate": alpha_b2, "B3_agod": alpha_b3}

        metrics = {"year": int(year), "msg": msg, "alpha": {}}
        for key, student in students.items():
            alpha = alphas[key]
            last_loss = student.distill(X_t, alpha, rng)
            Zs = student.embed(X_t, teacher=False)
            Zt = student.embed(X_t, teacher=True)
            Zs_ref = student.embed(X_ref, teacher=False)
            Zt_ref = student.embed(X_ref, teacher=True)
            rec = recall_at_k(Zs, Zt, k=TOP_K)
            rec_ref = recall_at_k(Zs_ref, Zt_ref, k=TOP_K)
            metrics["alpha"][key] = alpha.tolist()
            metrics[key] = {
                "recall_at_k_current": rec,
                "recall_at_k_reference": rec_ref,
                "forgetting": float(max(0.0, 1.0 - rec_ref)),
                "mse_global": last_loss.get("global"),
                "mse_routed": last_loss.get("routed"),
                "mse_mod": {m: last_loss.get(m) for m in MODALITIES},
                "step_counts": last_loss.get("step_counts"),
                "flops_share": alpha.tolist(),
            }
            print(
                f"[{name}] {year} {key}: alpha={np.round(alpha, 3).tolist()} "
                f"steps={last_loss.get('step_counts')} "
                f"MSE={last_loss.get('global'):.3f} val={last_loss.get('valence'):.3f} "
                f"R@{TOP_K}={rec:.3f}",
                flush=True,
            )
        rows.append(metrics)

    elapsed = time.perf_counter() - t0
    summary = {}
    for key in students:
        rec = [r[key]["recall_at_k_current"] for r in rows]
        forget = [r[key]["forgetting"] for r in rows]
        mse_g = [r[key]["mse_global"] for r in rows]
        mse_v = [r[key]["mse_mod"]["valence"] for r in rows]
        late = [r for r in rows if r["year"] in (1900, 1950)]
        summary[key] = {
            "mean_recall_at_k": float(np.mean(rec)),
            "mean_forgetting": float(np.mean(forget)),
            "mean_mse_global": float(np.mean(mse_g)),
            "mean_mse_valence": float(np.mean(mse_v)),
            "late_mse_valence": float(np.mean([r[key]["mse_mod"]["valence"] for r in late])),
            "last_recall_at_k": float(rec[-1]),
        }
    return {
        "protocol": name,
        "n_per": int(len(X_ref)),
        "d_mod": D_MOD,
        "d_out": D_OUT,
        "y_rule": "np.arange(n)",
        "modalities": list(MODALITIES),
        "elapsed_seconds": elapsed,
        "windows": rows,
        "summary": summary,
    }


def write_latex(results: dict) -> str:
    obs = results["observational"]
    inj = results["valence_inject"]
    lines = [
        "% ChronoBerg AGOD prototype · RF-domain MSG + dummy-label PO-risk\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{ChronoBerg AGOD prototype. Passages from era test splits; "
        "modalities are hashed text plus period-calibrated valence/arousal. "
        "PO-risk labels are $Y=\\texttt{np.arange}(n)$ (standardized). "
        "Distillation uses a fixed step budget routed by $\\alpha$. "
        "Late valence MSE is the 1900/1950 mean.}\n",
        "\\label{tab:chronoberg-agod-summary}\n\\small\n",
        "\\begin{tabular}{@{}l cc cc@{}}\n\\toprule\n",
        "Method & \\multicolumn{2}{c}{Observational} & "
        "\\multicolumn{2}{c}{Valence inject (1900/1950)} \\\\\n",
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n",
        " & Mean MSE & Late val. MSE & Mean MSE & Late val. MSE \\\\\n\\midrule\n",
    ]
    names = [
        ("B1\\_static", "B1_static"),
        ("B2\\_covariate", "B2_covariate"),
        ("B3\\_AGOD", "B3_agod"),
    ]
    for lab, key in names:
        s0, s1 = obs["summary"][key], inj["summary"][key]
        lines.append(
            f"{lab} & ${s0['mean_mse_global']:.3f}$ & ${s0['late_mse_valence']:.3f}$ & "
            f"${s1['mean_mse_global']:.3f}$ & ${s1['late_mse_valence']:.3f}$ \\\\\n"
        )
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    lines.append("\\begin{table}[ht]\n\\centering\n")
    lines.append(
        "\\caption{MSG state and AGOD routing weights on observational ChronoBerg "
        "(reference = 1750). $g_m=\\mathrm{Normalize}(\\mathrm{AUC}_m\\cdot"
        "\\mathrm{VIMP}_m+\\gamma\\cdot\\mathrm{PO}$-$\\mathrm{risk}_m)$.}\n"
    )
    lines.append("\\label{tab:chronoberg-agod-msg}\n\\small\n")
    lines.append("\\begin{tabular}{@{}l l ccc ccc@{}}\n\\toprule\n")
    lines.append(
        "Year & $W$ vs 1750 & "
        "AUC$_t$ & AUC$_v$ & AUC$_a$ & "
        "$\\alpha_t$ & $\\alpha_v$ & $\\alpha_a$ \\\\\n\\midrule\n"
    )
    for row in obs["windows"]:
        auc = row["msg"]["auc"]
        a = row["alpha"]["B3_agod"]
        lines.append(
            f"{row['year']} & era window & "
            f"${auc['text']:.3f}$ & ${auc['valence']:.3f}$ & ${auc['arousal']:.3f}$ & "
            f"${a[0]:.2f}$ & ${a[1]:.2f}$ & ${a[2]:.2f}$ \\\\\n"
        )
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
    return "".join(lines)


def write_prototype_py(results: dict) -> str:
    lines = [
        "# ChronoBerg AGOD prototype · MSG clever covariate, Y=np.arange(n)\n",
        f"modalities = {list(MODALITIES)!r}\n",
        "y_rule = 'np.arange(n)'\n",
    ]
    for proto in ("observational", "valence_inject"):
        r = results[proto]
        lines.append(f"\n# --- {proto} ---\n")
        lines.append(f"summary_{proto} = {repr(r['summary'])}\n")
        for row in r["windows"]:
            y = row["year"]
            lines.append(
                f"msg_{proto}_{y} = {repr(row['msg']['g'])}\n"
            )
            lines.append(
                f"alpha_agod_{proto}_{y} = {repr(row['alpha']['B3_agod'])}\n"
            )
            lines.append(
                f"auc_{proto}_{y} = {repr(row['msg']['auc'])}\n"
            )
    return "".join(lines)


def write_plot(results: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.4))
    years = [r["year"] for r in results["observational"]["windows"]]
    for ax, proto, title in [
        (axes[0], "observational", "Observational ChronoBerg"),
        (axes[1], "valence_inject", "Valence inject 1900/1950"),
    ]:
        rows = results[proto]["windows"]
        A = np.array([r["alpha"]["B3_agod"] for r in rows])
        for i, m in enumerate(MODALITIES):
            ax.plot(years, A[:, i], marker="o", label=m)
        ax.set_ylim(0, 1)
        ax.set_title(title)
        ax.set_xlabel("era")
        ax.set_ylabel(r"AGOD $\alpha$")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.3)
    ax = axes[2]
    keys = ["B1_static", "B2_covariate", "B3_agod"]
    labs = ["B1 static", "B2 cov.", "B3 AGOD"]
    x = np.arange(len(keys))
    w = 0.35
    obs = [results["observational"]["summary"][k]["late_mse_valence"] for k in keys]
    inj = [results["valence_inject"]["summary"][k]["late_mse_valence"] for k in keys]
    ax.bar(x - w / 2, obs, w, label="observational")
    ax.bar(x + w / 2, inj, w, label="valence inject")
    ax.set_xticks(x, labs)
    ax.set_ylabel("late valence MSE")
    ax.set_title("Routed distillation error")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-features", action="store_true")
    args = parser.parse_args()
    if args.rebuild_features or not NPZ.exists():
        save_npz(build_window_features())
    windows = load_windows()
    rng = _rng(SEED)
    print("=== observational ChronoBerg ===", flush=True)
    obs = run_protocol(windows, name="observational", seed=SEED)
    print("=== valence inject ===", flush=True)
    inj = run_protocol(inject_valence(windows, rng), name="valence_inject", seed=SEED + 1)
    payload = {
        "dataset": "spaul25/Chronoberg",
        "source": "era test splits chunked to 80-token passages",
        "modalities": list(MODALITIES),
        "y_rule": "np.arange(n) on stacked (D0, Dt) for PO-risk",
        "hyperparameters": {
            "n_per": N_PER,
            "d_mod": D_MOD,
            "gamma": GAMMA,
            "tau": TAU,
            "top_k": TOP_K,
            "step_budget": STEP_BUDGET,
            "inject_alpha": INJECT_ALPHA,
        },
        "observational": obs,
        "valence_inject": inj,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    (OUT / "chronoberg_agod_summary.json").write_text(json.dumps(payload, indent=2))
    (OUT / "chronoberg_agod_prototype.py").write_text(write_prototype_py(payload))
    tex = write_latex(payload)
    (OUT / "ChronoBerg_AGOD_tables_only.tex").write_text(tex)
    (DOCS / "ChronoBerg_AGOD_tables_only.tex").write_text(tex)
    write_plot(payload, OUT / "chronoberg_agod_routing.png")
    md = [
        "# ChronoBerg AGOD prototype\n\n",
        "Labels for PO-risk: `Y = np.arange(n)` on stacked reference + current window.\n\n",
        "| Protocol | B1 MSE | B2 MSE | B3 AGOD MSE | B3 late val. MSE |\n",
        "|----------|--------|--------|-------------|------------------|\n",
        f"| observational | {obs['summary']['B1_static']['mean_mse_global']:.3f} | "
        f"{obs['summary']['B2_covariate']['mean_mse_global']:.3f} | "
        f"{obs['summary']['B3_agod']['mean_mse_global']:.3f} | "
        f"{obs['summary']['B3_agod']['late_mse_valence']:.3f} |\n",
        f"| valence inject | {inj['summary']['B1_static']['mean_mse_global']:.3f} | "
        f"{inj['summary']['B2_covariate']['mean_mse_global']:.3f} | "
        f"{inj['summary']['B3_agod']['mean_mse_global']:.3f} | "
        f"{inj['summary']['B3_agod']['late_mse_valence']:.3f} |\n",
    ]
    (OUT / "README.md").write_text("".join(md))
    print(json.dumps(payload["observational"]["summary"], indent=2))
    print(json.dumps(payload["valence_inject"]["summary"], indent=2))
    print(f"elapsed obs={obs['elapsed_seconds']:.1f}s inject={inj['elapsed_seconds']:.1f}s")


if __name__ == "__main__":
    main()
