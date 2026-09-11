"""Labeled Hugging Face datasets → human feature-blocks → π-guided Stage-2.

Three labeled boards (Y is the dataset label; W is the batch/domain):

  * scikit-learn/adult-census-income  — Y = income, W = sex
      blocks: demography / work / capital
  * stanfordnlp/imdb vs rotten_tomatoes — Y = sentiment, W = source
      blocks: text / style / polarity
  * nyu-mll/multi_nli (fiction vs telephone) — Y = entailment, W = genre
      blocks: premise / hypothesis / overlap

Plots + a LaTeX dashboard drop into results/hf_block_aware/ and docs/method/.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Callable
from urllib.parse import quote

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
    shuffle_block,
)

TOKEN_RE = re.compile(r"[A-Za-z]{2,}")
POS_LEX = {
    "good", "great", "excellent", "love", "amazing", "wonderful", "best", "nice",
    "enjoy", "perfect", "happy", "positive", "beautiful", "brilliant", "fantastic",
    "superb", "delightful", "enjoyable", "impressive", "recommend",
}
NEG_LEX = {
    "bad", "terrible", "poor", "hate", "awful", "worst", "boring", "waste",
    "disappointing", "horrible", "stupid", "negative", "dull", "mediocre",
    "fails", "worse", "ugly", "predictable", "mess", "unfortunately",
}


def _parquet_url(dataset: str, config: str, split: str, shard: str = "0000.parquet") -> str:
    rev = quote("refs/convert/parquet", safe="")
    return f"https://huggingface.co/datasets/{dataset}/resolve/{rev}/{config}/{split}/{shard}"


def _read_parquet(dataset: str, config: str, split: str, columns=None) -> pd.DataFrame:
    url = _parquet_url(dataset, config, split)
    return pd.read_parquet(url, columns=columns)


def _balanced_take(mask: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return idx
    n = min(int(n), len(idx))
    return rng.choice(idx, n, replace=False)


def _onehot(series: pd.Series, prefix: str, max_levels: int = 10) -> tuple[np.ndarray, list[str]]:
    s = series.astype(str).str.strip().replace({"?": "UNK", "": "UNK"})
    keep = set(s.value_counts().head(max_levels).index)
    s = s.where(s.isin(keep), other="OTHER")
    dummies = pd.get_dummies(s, prefix=prefix, drop_first=False)
    return dummies.to_numpy(dtype=float), list(dummies.columns)


def encode_adult_blocks(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, ModalitySpec]:
    """X without sex/relationship (those leak W=sex); Y = income>50K."""
    workclass, _ = _onehot(df["workclass"], "wc", 8)
    occupation, _ = _onehot(df["occupation"], "occ", 10)
    marital, _ = _onehot(df["marital.status"], "mar", 6)
    race, _ = _onehot(df["race"], "race", 5)
    native = (df["native.country"].astype(str).str.strip() == "United-States").to_numpy(dtype=float).reshape(-1, 1)
    demo = np.hstack([
        df[["age"]].to_numpy(dtype=float),
        race,
        marital,
        native,
    ])
    work = np.hstack([
        df[["education.num", "hours.per.week"]].to_numpy(dtype=float),
        workclass,
        occupation,
    ])
    capital = df[["capital.gain", "capital.loss", "fnlwgt"]].to_numpy(dtype=float)
    scaler = StandardScaler()
    demo = scaler.fit_transform(demo)
    work = StandardScaler().fit_transform(work)
    capital = StandardScaler().fit_transform(capital)
    X = np.hstack([demo, work, capital])
    slices = [
        slice(0, demo.shape[1]),
        slice(demo.shape[1], demo.shape[1] + work.shape[1]),
        slice(demo.shape[1] + work.shape[1], X.shape[1]),
    ]
    spec = ModalitySpec(names=["demography", "work", "capital"], slices=slices)
    W = (df["sex"].astype(str).str.strip() == "Male").to_numpy(dtype=int)
    Y = df["income"].astype(str).str.contains(">50K", regex=False).to_numpy(dtype=float)
    return X, W, Y, spec


def load_adult(*, n_per_batch: int = 1000, seed: int = 2026) -> dict:
    df = _read_parquet("scikit-learn/adult-census-income", "default", "train")
    rng = np.random.default_rng(seed)
    sex = df["sex"].astype(str).str.strip()
    i0 = _balanced_take((sex == "Female").to_numpy(), n_per_batch, rng)
    i1 = _balanced_take((sex == "Male").to_numpy(), n_per_batch, rng)
    sub = df.iloc[np.concatenate([i0, i1])].reset_index(drop=True)
    X, W, Y, spec = encode_adult_blocks(sub)
    return {
        "X": X, "W": W, "Y": Y, "spec": spec,
        "n0": int((W == 0).sum()), "n1": int((W == 1).sum()),
        "dataset": "scikit-learn/adult-census-income",
        "label": "income >50K",
        "batch": "W = sex (Female vs Male); relationship dropped (Husband/Wife leak)",
        "slug": "adult",
        "gt_block": "capital",
        "title": "Adult census · income label · sex as batch",
    }


def _tokens(text: str) -> list[str]:
    return TOKEN_RE.findall(str(text).lower())


def style_features(texts: list[str]) -> np.ndarray:
    rows = []
    for t in texts:
        toks = _tokens(t)
        n = max(len(toks), 1)
        raw = str(t)
        nchar = max(len(raw), 1)
        punct = sum(ch in ".,;:!?" for ch in raw) / nchar
        upper = sum(ch.isupper() for ch in raw) / nchar
        digit = sum(ch.isdigit() for ch in raw) / nchar
        excl = raw.count("!") / nchar
        mean_len = float(np.mean([len(w) for w in toks])) if toks else 0.0
        rows.append([np.log1p(len(toks)), mean_len, punct, upper, digit, excl])
    return np.asarray(rows, dtype=float)


def polarity_features(texts: list[str]) -> np.ndarray:
    rows = []
    for t in texts:
        toks = _tokens(t)
        n = max(len(toks), 1)
        pos = sum(w in POS_LEX for w in toks)
        neg = sum(w in NEG_LEX for w in toks)
        rows.append([pos / n, neg / n, (pos - neg) / n, np.log1p(pos + neg)])
    return np.asarray(rows, dtype=float)


def hashed_text(texts: list[str], *, d_text: int, seed: int) -> np.ndarray:
    vec = HashingVectorizer(
        n_features=d_text, alternate_sign=False, ngram_range=(1, 2),
        lowercase=True, token_pattern=r"(?u)\b[a-zA-Z]{2,}\b",
    )
    return np.asarray(vec.fit_transform(texts).toarray(), dtype=float)


def overlap_features(prem: list[str], hyp: list[str]) -> np.ndarray:
    rows = []
    for a, b in zip(prem, hyp):
        ta, tb = set(_tokens(a)), set(_tokens(b))
        union = max(len(ta | tb), 1)
        inter = len(ta & tb) / union
        rows.append([
            np.log1p(len(ta)), np.log1p(len(tb)),
            np.log1p(abs(len(ta) - len(tb))),
            inter,
        ])
    return np.asarray(rows, dtype=float)


def _pack_text_blocks(
    texts0: list[str],
    texts1: list[str],
    y0: np.ndarray,
    y1: np.ndarray,
    *,
    d_text: int,
    seed: int,
) -> dict:
    texts = list(texts0) + list(texts1)
    text = hashed_text(texts, d_text=d_text, seed=seed)
    style = StandardScaler().fit_transform(style_features(texts))
    polar = StandardScaler().fit_transform(polarity_features(texts))
    X = np.hstack([text, style, polar])
    spec = ModalitySpec(
        names=["text", "style", "polarity"],
        slices=[slice(0, d_text), slice(d_text, d_text + 6), slice(d_text + 6, d_text + 10)],
    )
    W = np.array([0] * len(texts0) + [1] * len(texts1), dtype=int)
    Y = np.concatenate([y0, y1]).astype(float)
    return {"X": X, "W": W, "Y": Y, "spec": spec, "n0": len(texts0), "n1": len(texts1)}


def load_imdb_rt(*, n_per_batch: int = 800, d_text: int = 32, seed: int = 2026) -> dict:
    rng = np.random.default_rng(seed)
    imdb = _read_parquet("stanfordnlp/imdb", "plain_text", "test", columns=["text", "label"])
    rt = _read_parquet(
        "cornell-movie-review-data/rotten_tomatoes", "default", "train", columns=["text", "label"]
    )
    i_imdb = rng.choice(len(imdb), min(n_per_batch, len(imdb)), replace=False)
    i_rt = rng.choice(len(rt), min(n_per_batch, len(rt)), replace=False)
    a = imdb.iloc[i_imdb]
    b = rt.iloc[i_rt]
    pack = _pack_text_blocks(
        a["text"].astype(str).tolist(),
        b["text"].astype(str).tolist(),
        a["label"].to_numpy(dtype=float),
        b["label"].to_numpy(dtype=float),
        d_text=d_text, seed=seed,
    )
    pack.update({
        "dataset": "stanfordnlp/imdb vs cornell-movie-review-data/rotten_tomatoes",
        "label": "binary sentiment",
        "batch": "W = source (IMDB vs Rotten Tomatoes)",
        "slug": "imdb_rt",
        "gt_block": "polarity",
        "title": "IMDB vs Rotten Tomatoes · sentiment label · source as batch",
    })
    return pack


def load_mnli_genre(*, n_per_batch: int = 800, d_text: int = 24, seed: int = 2026) -> dict:
    rng = np.random.default_rng(seed)
    df = _read_parquet(
        "nyu-mll/multi_nli", "default", "validation_matched",
        columns=["premise", "hypothesis", "genre", "label"],
    )
    g = df["genre"].astype(str)
    fic = df[g == "fiction"]
    tel = df[g == "telephone"]
    i0 = rng.choice(len(fic), min(n_per_batch, len(fic)), replace=False)
    i1 = rng.choice(len(tel), min(n_per_batch, len(tel)), replace=False)
    a, b = fic.iloc[i0], tel.iloc[i1]
    prem = a["premise"].astype(str).tolist() + b["premise"].astype(str).tolist()
    hyp = a["hypothesis"].astype(str).tolist() + b["hypothesis"].astype(str).tolist()
    prem_x = hashed_text(prem, d_text=d_text, seed=seed)
    hyp_x = hashed_text(hyp, d_text=d_text, seed=seed + 1)
    ov = StandardScaler().fit_transform(overlap_features(prem, hyp))
    X = np.hstack([prem_x, hyp_x, ov])
    spec = ModalitySpec(
        names=["premise", "hypothesis", "overlap"],
        slices=[
            slice(0, d_text),
            slice(d_text, 2 * d_text),
            slice(2 * d_text, 2 * d_text + ov.shape[1]),
        ],
    )
    W = np.array([0] * len(a) + [1] * len(b), dtype=int)
    Y = np.concatenate([
        (a["label"].to_numpy() == 0).astype(float),
        (b["label"].to_numpy() == 0).astype(float),
    ])
    return {
        "X": X, "W": W, "Y": Y, "spec": spec,
        "n0": int((W == 0).sum()), "n1": int((W == 1).sum()),
        "dataset": "nyu-mll/multi_nli validation_matched",
        "label": "entailment vs rest",
        "batch": "W = genre (fiction vs telephone)",
        "slug": "mnli",
        "gt_block": "overlap",
        "title": "MultiNLI · entailment label · fiction vs telephone",
    }


def _mean_numeric_rows(rows: list[dict]) -> dict:
    out = dict(rows[0])
    keys = [k for k, v in rows[0].items() if isinstance(v, (int, float)) and not isinstance(v, bool)]
    for k in keys:
        if k.endswith("_sd") or k == "n_seeds":
            continue
        vals = [float(r[k]) for r in rows if k in r and r[k] is not None]
        out[k] = round(float(np.mean(vals)), 4)
        if len(vals) > 1:
            out[k + "_sd"] = round(float(np.std(vals, ddof=1)), 4)
    out["n_seeds"] = len(rows)
    return out


def _eval_setting(X, W, spec, *, seed, n_estimators, gt=None, light=True, with_subspace=False):
    gap = decompose_modality_gap(
        X, W, spec, seed=seed, n_estimators=n_estimators, light=light,
    )
    row = compare_raw_vs_clever(
        X, W, None, spec, gap, seed=seed, gt=gt, with_po=False,
        n_estimators=n_estimators, n_splits=5,
        with_block_train=True, with_subspace=with_subspace,
    )
    return row, gap


def _subsample(X, W, seed, frac=0.85):
    rng = np.random.default_rng(seed)
    i0 = np.where(W == 0)[0]
    i1 = np.where(W == 1)[0]
    n0 = max(80, int(frac * len(i0)))
    n1 = max(80, int(frac * len(i1)))
    sel = np.concatenate([
        rng.choice(i0, min(n0, len(i0)), replace=False),
        rng.choice(i1, min(n1, len(i1)), replace=False),
    ])
    return X[sel], W[sel]


def run_labeled_board(
    pack: dict,
    *,
    seed: int = 2026,
    n_estimators: int = 80,
    n_gt_seeds: int = 2,
    out_dir: Path | None = None,
    inject_alpha: float = 0.80,
) -> dict:
    slug = pack["slug"]
    out_dir = Path(out_dir or (ROOT / "results" / "hf_block_aware" / slug))
    out_dir.mkdir(parents=True, exist_ok=True)
    X, W, spec = pack["X"], pack["W"], pack["spec"]
    gt = pack["gt_block"]

    gap = decompose_modality_gap(X, W, spec, seed=seed, n_estimators=n_estimators, light=False)
    obs, _ = _eval_setting(X, W, spec, seed=seed, n_estimators=n_estimators, light=False)
    obs["name"] = "observational"
    obs["gt"] = None

    def _repeat(transform: Callable, name: str):
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

    comparisons = [
        _repeat(
            lambda Xs, Ws, s: inject_block_shift(Xs, Ws, spec, gt, alpha=inject_alpha),
            f"inject {gt} (α={inject_alpha:.2f})",
        ),
        _repeat(
            lambda Xs, Ws, s: inject_block_shift(
                shuffle_block(Xs, spec, spec.names[0], seed=s), Ws, spec, gt, alpha=inject_alpha + 0.05
            ),
            f"{spec.names[0]}⊥ + inject {gt}",
        ),
    ]

    summary = {
        "dataset": pack["dataset"],
        "title": pack["title"],
        "label": pack["label"],
        "batch": pack["batch"],
        "slug": slug,
        "gt_block": gt,
        "n0": pack["n0"],
        "n1": pack["n1"],
        "modalities": spec.names,
        "gap": gap.as_dict(),
        "observational_auc": obs,
        "comparisons": comparisons,
        "protocol": {
            "rf": f"max_depth=8, min_samples_leaf=5, n_estimators={n_estimators}",
            "auc": "5-fold stratified CV",
            "stage2": "π-weighted RF; opinion pool; X+Z; adaptive logit",
            "gt_seeds": n_gt_seeds,
        },
    }
    plots = plot_board(summary, out_dir)
    summary["plots"] = {k: str(v) for k, v in plots.items()}
    (out_dir / f"{slug}_board.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_dataset_latex(summary, out_dir / f"{slug}_tables.tex")
    return summary


def _tex_escape(s: str) -> str:
    return str(s).replace("_", r"\_").replace("%", r"\%")


def plot_board(summary: dict, out_dir: Path) -> dict[str, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    slug = summary["slug"]
    title = summary.get("title", slug)
    names = list(summary["gap"]["pi_consensus"].keys())

    def _vals(key):
        return [summary["gap"][key][n] for n in names]

    fig, ax = plt.subplots(figsize=(7.6, 4.0))
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
    ax.set_title(f"{title}\nStage-1 block gap")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / f"{slug}_gap_shares.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths["gap_shares"] = p

    gt_rows = [
        r for r in summary["comparisons"]
        if r.get("gt") and r.get("domain_auc_raw", 0) < 0.98
    ]
    if gt_rows:
        fig, ax = plt.subplots(figsize=(8.4, 4.2))
        labs = [r["name"] for r in gt_rows]
        x = np.arange(len(labs))
        raw = [r["domain_auc_raw"] for r in gt_rows]
        pirf = [r.get("domain_auc_pi_rf", r["domain_auc_raw"]) for r in gt_rows]
        pool = [r.get("domain_auc_pool", r["domain_auc_clever"]) for r in gt_rows]
        stack = [r.get("domain_auc_stack_xz", r["domain_auc_clever"]) for r in gt_rows]
        e_raw = [r.get("domain_auc_raw_sd", 0.0) for r in gt_rows]
        e_pi = [r.get("domain_auc_pi_rf_sd", 0.0) for r in gt_rows]
        e_s = [r.get("domain_auc_stack_sd", 0.0) for r in gt_rows]
        w = 0.18
        ax.bar(x - 1.5 * w, raw, w, yerr=e_raw, capsize=3, label="raw RF on X", color="#9aa0a6")
        ax.bar(x - 0.5 * w, pirf, w, yerr=e_pi, capsize=3, label="π-weighted RF", color="#ff7f0e")
        ax.bar(x + 0.5 * w, pool, w, yerr=0, capsize=3, label="π-opinion pool", color="#1f77b4")
        ax.bar(x + 1.5 * w, stack, w, yerr=e_s, capsize=3, label="X + clever-Z", color="#2ca02c")
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=10, ha="right")
        ax.set_ylabel("5-fold domain AUC")
        ax.set_ylim(0.45, 1.02)
        ax.axhline(0.5, ls="--", lw=0.8, color="0.5")
        ax.set_title(f"{title}\nπ-guided Stage-2")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / f"{slug}_stage2_auc.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["auc"] = p

        fig, ax = plt.subplots(figsize=(7.4, 4.0))
        x = np.arange(len(labs))
        ax.bar(x - 0.22, [r.get("mass_on_gt_raw", 0) for r in gt_rows], 0.22,
               yerr=[r.get("mass_on_gt_raw_sd", 0) for r in gt_rows], capsize=3,
               label="raw VIMP on GT", color="#9aa0a6")
        ax.bar(x, [r.get("pi_consensus_on_gt", 0) for r in gt_rows], 0.22,
               yerr=[r.get("pi_consensus_on_gt_sd", 0) for r in gt_rows], capsize=3,
               label="consensus π on GT", color="#1f77b4")
        ax.bar(x + 0.22, [r.get("pi_rf_on_gt", 0) for r in gt_rows], 0.22,
               yerr=[r.get("pi_rf_on_gt_sd", 0) for r in gt_rows], capsize=3,
               label="π-weighted RF VIMP on GT", color="#ff7f0e")
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=10, ha="right")
        ax.set_ylabel("mass on ground-truth block")
        ax.set_ylim(0, 1.02)
        ax.set_title(f"{title}\ninject recovery")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / f"{slug}_inject_recovery.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["inject"] = p
    return paths


def write_dataset_latex(summary: dict, path: Path) -> None:
    gap = summary["gap"]
    names = list(gap["pi_consensus"].keys())
    cap = _tex_escape(summary.get("title", summary["slug"]))
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        rf"\caption{{{cap}: Stage-1 block shares. $Y$={_tex_escape(summary.get('label',''))}; {_tex_escape(summary.get('batch',''))}.}}",
        r"\begin{tabular}{l" + "c" * len(names) + "}",
        r"\toprule",
        "Decomposition & " + " & ".join(_tex_escape(n) for n in names) + r" \\",
        r"\midrule",
    ]
    for key, lab in [
        ("pi_vimp", "RF VIMP"),
        ("pi_auc", "block AUC"),
        ("pi_mmd", "block MMD"),
        ("pi_lomo", "LOMO drop"),
        ("pi_consensus", r"consensus $\pi$"),
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
        rf"\caption{{{cap}: Stage-2 $\pi$-guided training (5-fold CV AUC).}}",
        r"\begin{tabular}{lcccccc}",
        r"\toprule",
        r"Setting & raw RF & $\pi$-RF & pool & $X{+}Z$ & $\Delta_{\pi\mathrm{-RF}}$ & $\pi$-RF on GT \\",
        r"\midrule",
    ]
    obs = summary.get("observational_auc") or {}
    if obs:
        pirf = obs.get("domain_auc_pi_rf", obs["domain_auc_raw"])
        lines.append(
            f"observational & {obs['domain_auc_raw']:.3f} & {pirf:.3f} & "
            f"{obs.get('domain_auc_pool', float('nan')):.3f} & "
            f"{obs.get('domain_auc_stack_xz', obs['domain_auc_clever']):.3f} & "
            f"{pirf - obs['domain_auc_raw']:+.3f} & --- \\\\"
        )
    for row in summary["comparisons"]:
        pirf = row.get("domain_auc_pi_rf", row["domain_auc_raw"])
        db = pirf - row["domain_auc_raw"]
        lines.append(
            f"{_tex_escape(row['name'])} & {row['domain_auc_raw']:.3f} & "
            f"{pirf:.3f} & {row.get('domain_auc_pool', float('nan')):.3f} & "
            f"{row.get('domain_auc_stack_xz', row['domain_auc_clever']):.3f} & "
            f"{db:+.3f} & {row.get('pi_rf_on_gt', float('nan')):.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_dashboard(summaries: list[dict], out_path: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(summaries)
    fig, axes = plt.subplots(n, 2, figsize=(12.4, 3.6 * n))
    if n == 1:
        axes = np.array([axes])
    for i, summary in enumerate(summaries):
        names = list(summary["gap"]["pi_consensus"].keys())
        ax = axes[i, 0]
        x = np.arange(len(names))
        ax.bar(x, [summary["gap"]["pi_consensus"][n] for n in names], color="#111111", width=0.55)
        ax.set_xticks(x)
        ax.set_xticklabels(names, fontsize=8)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("consensus π")
        ax.set_title(summary.get("title", summary["slug"]), fontsize=10, loc="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axes[i, 1]
        rows = [r for r in summary["comparisons"] if r.get("gt") and r.get("domain_auc_raw", 0) < 0.98]
        if not rows:
            ax.axis("off")
            continue
        labs = [r["name"] for r in rows]
        x = np.arange(len(labs))
        w = 0.22
        ax.bar(x - w, [r["domain_auc_raw"] for r in rows], w, color="#9aa0a6", label="raw RF")
        ax.bar(x, [r.get("domain_auc_pi_rf", r["domain_auc_raw"]) for r in rows], w, color="#ff7f0e", label="π-RF")
        ax.bar(x + w, [r.get("domain_auc_stack_xz", r["domain_auc_clever"]) for r in rows], w, color="#2ca02c", label="X+Z")
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=12, ha="right", fontsize=8)
        ax.set_ylim(0.45, 1.02)
        ax.axhline(0.5, ls="--", lw=0.7, color="0.5")
        ax.set_ylabel("5-fold AUC")
        if i == 0:
            ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Labeled HF boards · human blocks → π → Stage-2 training", fontsize=13, y=1.01)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _relpath_for_tex(png: Path) -> str:
    """Path from docs/method/ to a results PNG."""
    docs = ROOT / "docs" / "method"
    try:
        return png.resolve().relative_to(docs.resolve()).as_posix()
    except ValueError:
        return png.as_posix()


def write_dashboard_latex(summaries: list[dict], chronoberg: dict | None, path: Path) -> None:
    """Full article-style dashboard with figures placed next to captions."""
    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=1in]{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{caption}",
        r"\usepackage{subcaption}",
        r"\usepackage{hyperref}",
        r"\graphicspath{{../../results/chronoberg_clever_cov/}{../../results/hf_block_aware/}}",
        r"\title{Block-aware $\pi$ weights on labeled Hugging Face datasets}",
        r"\author{CFPerm prototype}",
        r"\date{}",
        r"\begin{document}",
        r"\maketitle",
        r"\noindent Analyst-defined feature-blocks and Stage-1 consensus $\pi$ guide Stage-2 training: ",
        r"$\pi$-weighted RF (column replicates so $\mathtt{max\_features}$ samples $p_j\propto\pi_m/|B_m|$), ",
        r"a linear opinion pool $\sum_m\pi_m\hat e_m$, and $X{+}Z$ stacking. ",
        r"$Y$ is the dataset label; $W$ is the batch/domain used for FSDS.",
        "",
        r"\section{Overview dashboard}",
        r"\begin{figure}[h]",
        r"\centering",
        r"\includegraphics[width=\textwidth]{dashboard.png}",
        r"\caption{Each row is one labeled Hugging Face board. Left: Stage-1 consensus $\pi$ over human blocks. ",
        r"Right: 5-fold domain AUC on concentrated inject-GT (raw RF vs $\pi$-weighted RF vs $X{+}Z$).}",
        r"\label{fig:hf-dashboard}",
        r"\end{figure}",
        "",
    ]
    if chronoberg:
        lines += [
            r"\section{Chronoberg (temporal text, VAD blocks)}",
            r"Temporal batch $W=1750$ vs.\ $1950$ on \texttt{spaul25/Chronoberg}; modalities hashed text / valence / arousal / dominance.",
            r"\begin{figure}[h]",
            r"\centering",
            r"\begin{subfigure}[t]{0.48\textwidth}",
            r"\includegraphics[width=\linewidth]{chronoberg_modality_gap_shares.png}",
            r"\caption{Stage-1 gap shares.}",
            r"\end{subfigure}\hfill",
            r"\begin{subfigure}[t]{0.48\textwidth}",
            r"\includegraphics[width=\linewidth]{chronoberg_clever_vs_raw_auc.png}",
            r"\caption{$\pi$-guided Stage-2 GT board.}",
            r"\end{subfigure}",
            r"\caption{Chronoberg 1750 vs.\ 1950. Observational shift is diffuse; wins are on concentrated valence inject.}",
            r"\label{fig:chrono}",
            r"\end{figure}",
            r"\begin{figure}[h]",
            r"\centering",
            r"\begin{subfigure}[t]{0.48\textwidth}",
            r"\includegraphics[width=\linewidth]{chronoberg_inject_recovery.png}",
            r"\caption{GT mass: raw VIMP vs $\pi$ vs $\pi$-RF.}",
            r"\end{subfigure}\hfill",
            r"\begin{subfigure}[t]{0.48\textwidth}",
            r"\includegraphics[width=\linewidth]{chronoberg_sample_efficiency.png}",
            r"\caption{Sample efficiency on the observational split.}",
            r"\end{subfigure}",
            r"\caption{Chronoberg localization and $n$-scaling.}",
            r"\end{figure}",
            "",
        ]
        tex_tables = ROOT / "docs" / "method" / "Chronoberg_CleverCov_tables_only.tex"
        if tex_tables.exists():
            lines.append(r"\input{Chronoberg_CleverCov_tables_only.tex}")
            lines.append("")

    for summary in summaries:
        slug = summary["slug"]
        cap = _tex_escape(summary.get("title", slug))
        lines += [
            rf"\section{{{cap}}}",
            rf"\noindent Dataset: \texttt{{{_tex_escape(summary['dataset'])}}}. ",
            rf"Label $Y$: {_tex_escape(summary.get('label',''))}. {_tex_escape(summary.get('batch',''))}. ",
            rf"$n_0={summary['n0']}$, $n_1={summary['n1']}$.",
            r"\begin{figure}[h]",
            r"\centering",
            r"\begin{subfigure}[t]{0.48\textwidth}",
            rf"\includegraphics[width=\linewidth]{{{slug}/{slug}_gap_shares.png}}",
            r"\caption{Stage-1 block gap.}",
            r"\end{subfigure}\hfill",
            r"\begin{subfigure}[t]{0.48\textwidth}",
            rf"\includegraphics[width=\linewidth]{{{slug}/{slug}_stage2_auc.png}}",
            r"\caption{Stage-2 domain AUC.}",
            r"\end{subfigure}",
            rf"\caption{{{cap}: human blocks and $\pi$-guided training. Observational rows can be diffuse; inject-GT is the concentrated-shift check.}}",
            rf"\label{{fig:{slug}}}",
            r"\end{figure}",
            r"\begin{figure}[h]",
            r"\centering",
            rf"\includegraphics[width=0.72\textwidth]{{{slug}/{slug}_inject_recovery.png}}",
            rf"\caption{{{cap}: mass on the ground-truth block after inject. $\pi$-weighted RF should concentrate relative to raw VIMP.}}",
            r"\end{figure}",
            "",
        ]
        table_file = ROOT / "results" / "hf_block_aware" / slug / f"{slug}_tables.tex"
        if table_file.exists():
            lines.append(table_file.read_text(encoding="utf-8"))
            lines.append("")

    lines += [
        r"\section{Cross-board notes}",
        r"\begin{itemize}",
        r"\item Trees ignore monotone column scaling; $\pi$ enters sklearn RF through column multiplicity.",
        r"\item Random-subspace BAWF is a different estimator class (not shown on these boards).",
        r"\item Adaptive group-logit can saturate on mean-shift inject; the RF analogue is the headline.",
        r"\item Observational Hugging Face shifts are often diffuse across blocks; the GT inject rows isolate whether $\pi$ guides training when the shift \emph{is} concentrated.",
        r"\end{itemize}",
        r"\end{document}",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_board_readme(summaries: list[dict], out_dir: Path) -> None:
    lines = [
        "# Labeled Hugging Face π-guided boards",
        "",
        "Human feature-blocks → Stage-1 consensus π → Stage-2 (π-weighted RF, opinion pool, X+Z).",
        "Y is the dataset label; W is the batch used for FSDS.",
        "",
        "Compile the dashboard with:",
        "",
        "```bash",
        "cd docs/method && pdflatex -interaction=nonstopmode BlockAware_HF_dashboard.tex",
        "```",
        "",
        "## Datasets",
        "",
    ]
    for s in summaries:
        lines.append(f"### {s.get('title', s['slug'])}")
        lines.append("")
        lines.append(f"- Hub: `{s['dataset']}`")
        lines.append(f"- Label Y: {s.get('label')}")
        lines.append(f"- Batch W: {s.get('batch')}")
        lines.append(f"- n0={s['n0']}, n1={s['n1']}")
        pi = s["gap"]["pi_consensus"]
        lines.append("- π: " + ", ".join(f"{k} {v:.3f}" for k, v in pi.items()))
        obs = s.get("observational_auc") or {}
        if obs:
            lines.append(
                f"- Observational AUC raw {obs.get('domain_auc_raw', float('nan')):.3f} vs "
                f"π-RF {obs.get('domain_auc_pi_rf', float('nan')):.3f}"
            )
        lines.append("")
        lines.append("| setting | raw | π-RF | pool | X+Z | Δ π-RF | π-RF on GT |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for row in s["comparisons"]:
            pirf = row.get("domain_auc_pi_rf", row["domain_auc_raw"])
            lines.append(
                f"| {row['name']} | {row['domain_auc_raw']:.3f} | {pirf:.3f} | "
                f"{row.get('domain_auc_pool', float('nan')):.3f} | "
                f"{row.get('domain_auc_stack_xz', row['domain_auc_clever']):.3f} | "
                f"{pirf - row['domain_auc_raw']:+.3f} | {row.get('pi_rf_on_gt', float('nan')):.3f} |"
            )
        lines.append("")
    (out_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_all(
    *,
    n_per_batch: int = 800,
    d_text: int = 32,
    seed: int = 2026,
    n_estimators: int = 80,
    n_gt_seeds: int = 2,
    datasets: tuple[str, ...] = ("adult", "imdb_rt", "mnli"),
) -> dict:
    out_root = ROOT / "results" / "hf_block_aware"
    out_root.mkdir(parents=True, exist_ok=True)
    loaders = {
        "adult": lambda: load_adult(n_per_batch=n_per_batch, seed=seed),
        "imdb_rt": lambda: load_imdb_rt(n_per_batch=n_per_batch, d_text=d_text, seed=seed),
        "mnli": lambda: load_mnli_genre(n_per_batch=n_per_batch, d_text=max(16, d_text - 8), seed=seed),
    }
    summaries = []
    for key in datasets:
        pack = loaders[key]()
        summary = run_labeled_board(
            pack, seed=seed, n_estimators=n_estimators, n_gt_seeds=n_gt_seeds,
            out_dir=out_root / pack["slug"],
        )
        summaries.append(summary)

    dash = plot_dashboard(summaries, out_root / "dashboard.png")
    chrono_json = ROOT / "results" / "chronoberg_clever_cov" / "chronoberg_clever_cov.json"
    chrono = json.loads(chrono_json.read_text(encoding="utf-8")) if chrono_json.exists() else None
    tex_path = ROOT / "docs" / "method" / "BlockAware_HF_dashboard.tex"
    write_dashboard_latex(summaries, chrono, tex_path)
    # also keep a copy next to the PNGs
    write_dashboard_latex(summaries, chrono, out_root / "BlockAware_HF_dashboard.tex")
    write_board_readme(summaries, out_root)
    board = {
        "datasets": [s["slug"] for s in summaries],
        "dashboard_png": str(dash),
        "dashboard_tex": str(tex_path),
        "summaries": {
            s["slug"]: {
                "title": s["title"],
                "pi": s["gap"]["pi_consensus"],
                "obs_raw": s["observational_auc"]["domain_auc_raw"],
                "obs_pi_rf": s["observational_auc"].get("domain_auc_pi_rf"),
                "comparisons": [
                    {
                        "name": r["name"],
                        "raw": r["domain_auc_raw"],
                        "pi_rf": r.get("domain_auc_pi_rf"),
                        "stack": r.get("domain_auc_stack_xz"),
                        "pi_rf_on_gt": r.get("pi_rf_on_gt"),
                    }
                    for r in s["comparisons"]
                ],
            }
            for s in summaries
        },
    }
    (out_root / "board.json").write_text(json.dumps(board, indent=2), encoding="utf-8")
    return board
