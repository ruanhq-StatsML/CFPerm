"""Chronoberg prototype: modality-gap clever covariates for temporal shift.

Downloads the 1750 vs 1950 test books + period VAD lexicons from Hugging Face
(`spaul25/Chronoberg`), scores sentences on text / valence / arousal / dominance,
decomposes each modality's relative contribution to the batch shift, and feeds
those contributions in as TMLE-style clever covariates.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
SRC = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC))

from clever_covariate_gap import (  # noqa: E402
    ModalitySpec,
    compare_raw_vs_clever,
    decompose_modality_gap,
    inject_block_shift,
    make_synthetic_shift,
    sample_efficiency_curve,
    shuffle_block,
)

TOKEN_RE = re.compile(r"[A-Za-z]{2,}")
SENT_SPLIT = re.compile(r"(?<=[.!?])\s+")


def _hf_download(filename: str) -> Path:
    from huggingface_hub import hf_hub_download

    p = hf_hub_download(repo_id="spaul25/Chronoberg", filename=filename, repo_type="dataset")
    return Path(p)


def load_lexicon(path: Path) -> dict[str, float]:
    df = pd.read_csv(path)
    word_col = "words" if "words" in df.columns else df.columns[1]
    score_cols = [
        c
        for c in df.columns
        if str(c).lower() not in {"unnamed: 0", "words", "word"} and not str(c).startswith("Unnamed")
    ]
    score_col = score_cols[0]
    words = df[word_col].astype(str).str.lower().str.strip()
    scores = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
    out: dict[str, float] = {}
    for w, s in zip(words.tolist(), scores.astype(float).tolist()):
        if not w or w == "nan":
            continue
        out[w] = float(s)
    return out


def pool_lexicons(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    out: dict[str, float] = {}
    for w, sa in a.items():
        out[str(w)] = float(sa)
    for w, sb in b.items():
        w = str(w)
        out[w] = 0.5 * (out[w] + float(sb)) if w in out else float(sb)
    return out


def sentence_split(text: str) -> list[str]:
    text = re.sub(r"\s+", " ", str(text)).strip()
    parts = SENT_SPLIT.split(text)
    out = []
    for s in parts:
        s = s.strip()
        n = len(TOKEN_RE.findall(s))
        if 8 <= n <= 60:
            # drop Gutenberg boilerplate
            low = s.lower()
            if "project gutenberg" in low or "end of the project" in low:
                continue
            out.append(s)
    return out


def load_period_sentences(era: str, *, max_docs: int | None = None) -> list[str]:
    path = _hf_download(f"train_test_split/test/{era}_test.json")
    docs = json.loads(path.read_text(encoding="utf-8"))
    if max_docs is not None:
        docs = docs[:max_docs]
    sents: list[str] = []
    for d in docs:
        sents.extend(sentence_split(d.get("text", "")))
    return sents


def vad_stats(tokens: list[str], lex: dict[str, float]) -> np.ndarray:
    vals = [lex[t] for t in tokens if t in lex]
    if not vals:
        return np.zeros(5, dtype=float)
    arr = np.asarray(vals, dtype=float)
    cov = len(vals) / max(len(tokens), 1)
    return np.array(
        [float(arr.mean()), float(arr.std()), float(arr.min()), float(arr.max()), float(cov)],
        dtype=float,
    )


def build_chronoberg_features(
    *,
    n_per_batch: int = 1800,
    d_text: int = 32,
    seed: int = 2026,
) -> dict:
    """X = [text(d) | valence(5) | arousal(5) | dominance(5)], W = 1750 vs 1950.

    VAD uses a *static pooled* 1750+1950 lexicon so the scores reflect usage,
    not a circular period-specific recoding. Y is the sentence valence mean
    (held out of the CD feature set separately).
    """
    s0 = load_period_sentences("1750")
    s1 = load_period_sentences("1950")
    rng = np.random.default_rng(seed)
    n0 = min(len(s0), n_per_batch)
    n1 = min(len(s1), n_per_batch)
    s0 = [s0[i] for i in rng.choice(len(s0), n0, replace=False)]
    s1 = [s1[i] for i in rng.choice(len(s1), n1, replace=False)]
    texts = s0 + s1
    W = np.concatenate([np.zeros(n0, dtype=int), np.ones(n1, dtype=int)])

    vec = HashingVectorizer(
        n_features=d_text,
        alternate_sign=False,
        norm="l2",
        ngram_range=(1, 2),
        lowercase=True,
        token_pattern=r"(?u)\b[a-zA-Z]{2,}\b",
    )
    X_text = vec.transform(texts).toarray().astype(float)

    v1750 = load_lexicon(_hf_download("Lexicons/Valence_lexicon_1750.csv"))
    v1950 = load_lexicon(_hf_download("Lexicons/Valence_lexicon_1950.csv"))
    a1750 = load_lexicon(_hf_download("Lexicons/Arousal_lexicon_1750.csv"))
    a1950 = load_lexicon(_hf_download("Lexicons/Arousal_lexicon_1950.csv"))
    d1750 = load_lexicon(_hf_download("Lexicons/Dominance_lexicon_1750.csv"))
    d1950 = load_lexicon(_hf_download("Lexicons/Dominance_lexicon_1950.csv"))
    lex_v = pool_lexicons(v1750, v1950)
    lex_a = pool_lexicons(a1750, a1950)
    lex_d = pool_lexicons(d1750, d1950)

    V = np.zeros((len(texts), 5))
    A = np.zeros((len(texts), 5))
    D = np.zeros((len(texts), 5))
    for i, s in enumerate(texts):
        tok = [t.lower() for t in TOKEN_RE.findall(s)]
        V[i] = vad_stats(tok, lex_v)
        A[i] = vad_stats(tok, lex_a)
        D[i] = vad_stats(tok, lex_d)

    X = np.hstack([X_text, V, A, D])
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    Y = V[:, 0].copy()  # valence mean (pre-scale stored above)
    names = ["text", "valence", "arousal", "dominance"]
    slices = [
        slice(0, d_text),
        slice(d_text, d_text + 5),
        slice(d_text + 5, d_text + 10),
        slice(d_text + 10, d_text + 15),
    ]
    spec = ModalitySpec(names=names, slices=slices)
    return {
        "X": X,
        "W": W,
        "Y": Y,
        "spec": spec,
        "n0": n0,
        "n1": n1,
        "d_text": d_text,
        "texts": texts,
    }


def _cd_design(X: np.ndarray, spec: ModalitySpec) -> tuple[np.ndarray, ModalitySpec]:
    """Drop valence from X when Y is valence (no leakage)."""
    keep_names = [n for n in spec.names if n != "valence"]
    keep_slices = [spec.slices[spec.names.index(n)] for n in keep_names]
    cols = []
    new_slices = []
    start = 0
    for sl in keep_slices:
        block = X[:, sl]
        cols.append(block)
        new_slices.append(slice(start, start + block.shape[1]))
        start += block.shape[1]
    return np.hstack(cols), ModalitySpec(names=keep_names, slices=new_slices)


def plot_results(summary: dict, out_dir: Path) -> dict[str, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    names = summary["gap"]["pi_consensus"].keys()
    names = list(names)

    def _vals(key):
        return [summary["gap"][key][n] for n in names]

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    x = np.arange(len(names))
    width = 0.16
    series = [
        ("pi_vimp", "RF VIMP share"),
        ("pi_auc", "block AUC−0.5"),
        ("pi_mmd", "block MMD"),
        ("pi_lomo", "LOMO AUC drop"),
        ("pi_consensus", "consensus π"),
    ]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#111111"]
    for i, ((key, lab), c) in enumerate(zip(series, colors)):
        ax.bar(x + (i - 2) * width, _vals(key), width, label=lab, color=c, alpha=0.9 if key == "pi_consensus" else 0.75)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("relative contribution")
    ax.set_ylim(0, 1.02)
    ax.set_title("Chronoberg 1750 vs 1950 · modality-specific gap")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / "chronoberg_modality_gap_shares.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths["gap_shares"] = p

    # GT board: domain AUC with fold/seed error bars
    gt_rows = [r for r in summary["comparisons"] if r.get("gt")]
    if gt_rows:
        fig, ax = plt.subplots(figsize=(8.4, 4.4))
        labs = [r["name"] for r in gt_rows]
        x = np.arange(len(labs))
        raw = [r["domain_auc_raw"] for r in gt_rows]
        clever = [r["domain_auc_clever"] for r in gt_rows]
        stack = [r.get("domain_auc_stack_xz", r["domain_auc_clever"]) for r in gt_rows]
        e_raw = [r.get("domain_auc_raw_sd", 0.0) for r in gt_rows]
        e_z = [r.get("domain_auc_clever_sd", 0.0) for r in gt_rows]
        e_s = [r.get("domain_auc_stack_sd", 0.0) for r in gt_rows]
        ax.bar(x - 0.25, raw, 0.24, yerr=e_raw, capsize=3, label="raw X", color="#9aa0a6")
        ax.bar(x, clever, 0.24, yerr=e_z, capsize=3, label="clever-Z = π·logit ê", color="#1f77b4")
        ax.bar(x + 0.25, stack, 0.24, yerr=e_s, capsize=3, label="X + clever-Z", color="#2ca02c")
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=12, ha="right")
        ax.set_ylabel("5-fold domain AUC")
        ax.set_ylim(0.45, 1.02)
        ax.axhline(0.5, ls="--", lw=0.8, color="0.5")
        ax.set_title("Tuned GT board · concentrated shift recovery")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / "chronoberg_clever_vs_raw_auc.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["auc"] = p

    se = summary.get("sample_efficiency") or []
    if se:
        fig, ax = plt.subplots(figsize=(7.4, 4.3))
        ns = [r["n"] for r in se]
        ax.errorbar(ns, [r["auc_raw"] for r in se], yerr=[r.get("auc_raw_sd", 0) for r in se],
                    fmt="o-", color="#9aa0a6", label="raw X", capsize=3)
        ax.errorbar(ns, [r["auc_clever"] for r in se], yerr=[r.get("auc_clever_sd", 0) for r in se],
                    fmt="s-", color="#1f77b4", label="clever-Z", capsize=3)
        ax.errorbar(ns, [r.get("auc_stack", r["auc_clever"]) for r in se],
                    yerr=[r.get("auc_stack_sd", 0) for r in se],
                    fmt="^-", color="#2ca02c", label="X + clever-Z", capsize=3)
        ax.set_xlabel("n (balanced 1750+1950)")
        ax.set_ylabel("5-fold domain AUC (mean ± sd)")
        ax.set_title("Sample efficiency · multi-seed Chronoberg")
        ax.legend(frameon=False)
        ax.set_ylim(0.52, 0.85)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / "chronoberg_sample_efficiency.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["efficiency"] = p

    inj = [r for r in summary["comparisons"] if r.get("gt")]
    if inj:
        fig, ax = plt.subplots(figsize=(7.6, 4.3))
        labs = [r["name"] for r in inj]
        x = np.arange(len(labs))
        ax.bar(x - 0.22, [r.get("mass_on_gt_raw", 0) for r in inj], 0.22,
               yerr=[r.get("mass_on_gt_raw_sd", 0) for r in inj], capsize=3,
               label="raw VIMP mass on GT", color="#9aa0a6")
        ax.bar(x, [r.get("pi_consensus_on_gt", 0) for r in inj], 0.22,
               yerr=[r.get("pi_consensus_on_gt_sd", 0) for r in inj], capsize=3,
               label="consensus π on GT", color="#1f77b4")
        ax.bar(x + 0.22, [r.get("z_logit_on_gt", 0) for r in inj], 0.22,
               yerr=[r.get("z_logit_on_gt_sd", 0) for r in inj], capsize=3,
               label="clever-Z VIMP on GT", color="#2ca02c")
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=12, ha="right")
        ax.set_ylabel("mass on ground-truth modality")
        ax.set_ylim(0, 1.02)
        ax.set_title("Controlled inject recovery · GT modality mass")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / "chronoberg_inject_recovery.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["inject"] = p
    return paths


def _tex_escape(s: str) -> str:
    return str(s).replace("_", r"\_")


def write_latex(summary: dict, path: Path) -> None:
    gap = summary["gap"]
    names = list(gap["pi_consensus"].keys())
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        r"\caption{Chronoberg 1750 vs.\ 1950: modality-specific gap shares and clever-covariate FSDS.}",
        r"\begin{tabular}{l" + "c" * len(names) + "}",
        r"\toprule",
        "Decomposition & " + " & ".join(names) + r" \\",
        r"\midrule",
    ]
    for key, lab in [
        ("pi_vimp", "RF VIMP"),
        ("pi_auc", "block AUC"),
        ("pi_mmd", "block MMD"),
        ("pi_lomo", "LOMO drop"),
        ("pi_consensus", "consensus $\\pi$"),
    ]:
        vals = " & ".join(f"{gap[key][n]:.3f}" for n in names)
        lines.append(f"{lab} & {vals} \\\\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        r"\caption{Raw $X$ vs.\ clever covariates $Z_m=\pi_m\log\hat e_m$ (5-fold CV AUC).}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Setting & AUC raw $X$ & AUC clever-$Z$ & AUC $X{+}Z$ & $\Delta_{X+Z}$ \\",
        r"\midrule",
    ]
    for row in summary["comparisons"]:
        dstack = row.get("domain_auc_stack_xz", row["domain_auc_clever"]) - row["domain_auc_raw"]
        lines.append(
            f"{_tex_escape(row['name'])} & {row['domain_auc_raw']:.3f} & "
            f"{row['domain_auc_clever']:.3f} & {row.get('domain_auc_stack_xz', row['domain_auc_clever']):.3f} & "
            f"{dstack:+.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _mean_numeric_rows(rows: list[dict]) -> dict:
    """Average numeric fields across seeds; keep nested dicts from the first row."""
    out = dict(rows[0])
    keys = []
    for k, v in rows[0].items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            keys.append(k)
    skip = {"n_seeds"}
    for k in keys:
        if k.endswith("_sd") or k in skip:
            continue
        vals = [float(r[k]) for r in rows if k in r and r[k] is not None]
        out[k] = round(float(np.mean(vals)), 4)
        if len(vals) > 1:
            out[k + "_sd"] = round(float(np.std(vals, ddof=1)), 4)
    out["n_seeds"] = len(rows)
    return out


def _eval_setting(X, W, spec, *, seed, n_estimators, gt=None, light=True):
    gap = decompose_modality_gap(
        X, W, spec, seed=seed, n_estimators=n_estimators, light=light,
    )
    row = compare_raw_vs_clever(
        X, W, None, spec, gap, seed=seed, gt=gt, with_po=False,
        n_estimators=n_estimators, n_splits=5,
    )
    return row, gap


def run_prototype(
    *,
    n_per_batch: int = 1000,
    d_text: int = 48,
    seed: int = 2026,
    n_estimators: int = 100,
    out_dir: Path | None = None,
    include_synthetic: bool = True,
    efficiency_ns: tuple[int, ...] = (200, 400, 800, 1200),
    n_repeats: int = 4,
    n_gt_seeds: int = 3,
) -> dict:
    out_dir = Path(out_dir or (ROOT / "results" / "chronoberg_clever_cov"))
    out_dir.mkdir(parents=True, exist_ok=True)

    pack = build_chronoberg_features(n_per_batch=n_per_batch, d_text=d_text, seed=seed)
    X, W, Y, spec = pack["X"], pack["W"], pack["Y"], pack["spec"]

    # Full LOMO gap on the observational split (scientific board).
    gap = decompose_modality_gap(X, W, spec, seed=seed, n_estimators=n_estimators, light=False)
    obs, _ = _eval_setting(X, W, spec, seed=seed, n_estimators=n_estimators, light=False)
    obs["name"] = "observational 1750/1950"
    obs["gt"] = None

    comparisons = []

    def _subsample(Xin, Win, s, frac=0.85):
        rng = np.random.default_rng(s)
        i0 = np.where(Win == 0)[0]
        i1 = np.where(Win == 1)[0]
        n0 = max(120, int(frac * len(i0)))
        n1 = max(120, int(frac * len(i1)))
        sel = np.concatenate([
            rng.choice(i0, n0, replace=False),
            rng.choice(i1, n1, replace=False),
        ])
        return Xin[sel], Win[sel]

    def _repeat_gt(transform, name, gt):
        rows = []
        for s in range(n_gt_seeds):
            sid = seed + 31 * (s + 1)
            Xs, Ws = _subsample(X, W, sid)
            Xt = transform(Xs, Ws, sid)
            row, _g = _eval_setting(
                Xt, Ws, spec, seed=seed + 7 * s, n_estimators=n_estimators, gt=gt, light=True,
            )
            rows.append(row)
        agg = _mean_numeric_rows(rows)
        agg["name"] = name
        agg["gt"] = gt
        return agg

    comparisons.append(_repeat_gt(
        lambda Xs, Ws, s: inject_block_shift(Xs, Ws, spec, "valence", alpha=0.70),
        "inject valence (α=0.70)",
        "valence",
    ))
    comparisons.append(_repeat_gt(
        lambda Xs, Ws, s: inject_block_shift(
            shuffle_block(Xs, spec, "text", seed=s), Ws, spec, "valence", alpha=0.80
        ),
        "text⊥ + inject valence",
        "valence",
    ))
    comparisons.append(_repeat_gt(
        lambda Xs, Ws, s: inject_block_shift(Xs, Ws, spec, "text", alpha=0.35),
        "inject text (α=0.35)",
        "text",
    ))

    se = sample_efficiency_curve(
        X, W, spec, ns=efficiency_ns, seed=seed,
        n_estimators=max(80, n_estimators - 20), n_repeats=n_repeats,
    )

    synthetic = None
    if include_synthetic:
        syn_rows = []
        for s in range(n_gt_seeds):
            Xs, Ws, _Ys, specs = make_synthetic_shift(
                n=280, d_text=48, d_vad=4, gt="valence", mean_shift=0.90, seed=seed + 100 * s,
            )
            row, g = _eval_setting(Xs, Ws, specs, seed=seed + s, n_estimators=n_estimators, gt="valence", light=True)
            syn_rows.append(row)
            last_g = g
        synthetic = _mean_numeric_rows(syn_rows)
        synthetic["name"] = "synthetic GT=valence"
        synthetic["gt"] = "valence"
        synthetic["gap"] = last_g.as_dict()
        comparisons.append(synthetic)

    summary = {
        "dataset": "spaul25/Chronoberg test 1750 vs 1950",
        "protocol": {
            "d_text": d_text,
            "n_per_batch": n_per_batch,
            "rf": "max_depth=8, min_samples_leaf=5, n_estimators=%d" % n_estimators,
            "auc": "5-fold stratified CV",
            "clever_Z": "Z_m = π_m * logit ê_m(X_m)",
            "gt_seeds": n_gt_seeds,
            "efficiency_repeats": n_repeats,
        },
        "n0": pack["n0"],
        "n1": pack["n1"],
        "d_text": d_text,
        "modalities": spec.names,
        "gap": gap.as_dict(),
        "observational_auc": obs,
        "comparisons": comparisons,
        "sample_efficiency": se,
        "synthetic": synthetic,
        "notes": {
            "clever_Z": "Z_m = π_m · logit ê_m(X_m); OOF propensity, 5-fold CV AUC",
            "clever_H": "TMLE H_m = π_m (W-ê_m)/(ê_m(1-ê_m)); PO targeting only",
            "vad_lexicon": "static pool of 1750+1950 Chronoberg VAD lexicons",
            "observational": "diffuse multi-modality shift; GT board is inject/synthetic",
        },
    }
    (out_dir / "chronoberg_clever_cov.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_latex(summary, out_dir / "Chronoberg_CleverCov_tables_only.tex")
    docs = ROOT / "docs" / "method"
    docs.mkdir(parents=True, exist_ok=True)
    write_latex(summary, docs / "Chronoberg_CleverCov_tables_only.tex")
    plots = plot_results(summary, out_dir)
    summary["plots"] = {k: str(v) for k, v in plots.items()}
    (out_dir / "chronoberg_clever_cov.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _write_readme(summary, out_dir)
    old_po = out_dir / "chronoberg_porisk_clever.png"
    if old_po.exists():
        old_po.unlink()
    return summary


def _write_readme(summary: dict, out_dir: Path) -> None:
    gap = summary["gap"]
    lines = [
        "# Chronoberg clever-covariate prototype",
        "",
        "Temporal batch **W = 1750 vs 1950** on Hugging Face `spaul25/Chronoberg` test books.",
        "Modalities: hashed text n-grams, valence / arousal / dominance (pooled lexicons).",
        "",
        "## Consensus relative contribution π",
        "",
        "| modality | π |",
        "|---|---:|",
    ]
    for k, v in gap["pi_consensus"].items():
        lines.append(f"| {k} | {v:.3f} |")
    lines += [
        "",
        "## GT board (mean over seeds, 5-fold CV AUC)",
        "",
        "| setting | AUC raw | AUC Z | AUC X+Z | Δ stack | π on GT |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary["comparisons"]:
        st = row.get("domain_auc_stack_xz", row["domain_auc_clever"])
        pi_gt = row.get("pi_consensus_on_gt", float("nan"))
        lines.append(
            f"| {row['name']} | {row['domain_auc_raw']:.3f} | {row['domain_auc_clever']:.3f} | "
            f"{st:.3f} | {st - row['domain_auc_raw']:+.3f} | {pi_gt:.3f} |"
        )
    obs = summary.get("observational_auc") or {}
    if obs:
        lines += [
            "",
            f"Observational 1750/1950 (no GT, diffuse shift): raw AUC "
            f"{obs.get('domain_auc_raw', float('nan')):.3f} vs clever-Z "
            f"{obs.get('domain_auc_clever', float('nan')):.3f}.",
        ]
    lines += [
        "",
        "Clever-Z is `Z_m = π_m · logit ê_m(X_m)` (n_modalities columns), estimated OOF.",
        "Detection uses 5-fold stratified CV. GT rows average 3 subsample seeds.",
        "",
        "```bash",
        "python3 scripts/run_chronoberg_clever_cov.py",
        "```",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    run_prototype()
