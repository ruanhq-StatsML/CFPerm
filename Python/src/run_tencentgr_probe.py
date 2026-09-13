#!/usr/bin/env python3
"""Probe the dumped TencentGR subset. No training, no attribution."""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.config import DEFAULT_CFG
from tencentgr.dataset import _events_from_row


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--out-dir", default="experiments/tencentgr/probe")
    return p.parse_args()


def _action_counts(vals: np.ndarray) -> dict:
    out = {"n": int(vals.size), "exposure_0": 0, "click_1": 0, "conversion_2": 0, "other": 0}
    for k, lab in [(0, "exposure_0"), (1, "click_1"), (2, "conversion_2")]:
        out[lab] = int((vals == k).sum())
    out["other"] = int(out["n"] - out["exposure_0"] - out["click_1"] - out["conversion_2"])
    if out["n"]:
        out["click_or_conv_rate"] = (out["click_1"] + out["conversion_2"]) / out["n"]
    else:
        out["click_or_conv_rate"] = 0.0
    return out


def _len_or_zero(v) -> int:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return 0
    try:
        return int(len(v))
    except TypeError:
        return 0


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    seq_df = pd.read_parquet(cache / "seq_df.parquet")
    user_feat = pd.read_parquet(cache / "user_feat.parquet")
    item_feat = pd.read_parquet(cache / "item_feat.parquet")
    samples = pickle.load(open(cache / "samples.pkl", "rb"))
    blob = np.load(cache / "mm_embeddings.npz", allow_pickle=False)
    mm_rids = set(int(x) for x in blob["rids"].tolist())
    embs = blob["embs"]

    raw_lens = []
    n_events = 0
    all_actions = []
    last_actions = []
    first_actions = []
    pos_from_end = {k: [] for k in range(1, 101)}
    ts_min = []
    ts_max = []
    ts_span = []
    unique_items = []
    item_hits = []
    click_positions = []

    for row in seq_df.itertuples(index=False):
        events = _events_from_row(getattr(row, "seq", None))
        raw_lens.append(len(events))
        if not events:
            continue
        acts = [int(e.get("action_type", 0)) for e in events]
        items = [int(e.get("item_id", 0)) for e in events]
        ts = [int(e.get("timestamp", 0)) for e in events]
        all_actions.extend(acts)
        last_actions.append(acts[-1])
        first_actions.append(acts[0])
        unique_items.append(len(set(items)))
        item_hits.extend(items)
        n_events += len(acts)
        ts_min.append(min(ts))
        ts_max.append(max(ts))
        ts_span.append(max(ts) - min(ts))
        for k, a in enumerate(reversed(acts), start=1):
            if k <= 100:
                pos_from_end[k].append(a)
        for i, a in enumerate(acts):
            if a >= 1:
                click_positions.append(i / max(len(acts) - 1, 1))

    all_actions = np.asarray(all_actions, dtype=np.int64)
    last_actions = np.asarray(last_actions, dtype=np.int64)
    first_actions = np.asarray(first_actions, dtype=np.int64)
    raw_lens = np.asarray(raw_lens, dtype=np.int64)
    unique_items = np.asarray(unique_items, dtype=np.int64)
    ts_min_a = np.asarray(ts_min, dtype=np.int64)
    ts_max_a = np.asarray(ts_max, dtype=np.int64)
    ts_span_a = np.asarray(ts_span, dtype=np.int64)
    item_hits = np.asarray(item_hits, dtype=np.int64)

    pos_click_rate = []
    for k in range(1, 101):
        arr = np.asarray(pos_from_end[k], dtype=np.int64)
        if arr.size == 0:
            break
        pos_click_rate.append(
            {
                "k_from_end": k,
                "n": int(arr.size),
                "click_or_conv_rate": float((arr >= 1).mean()),
                "conversion_rate": float((arr == 2).mean()),
            }
        )

    sample_hist_lens = np.asarray([len(s["history_items"]) for s in samples], dtype=np.int64)
    sample_tgt = np.asarray([int(s["target_action"]) for s in samples], dtype=np.int64)
    sample_any = []
    sample_rate = []
    sample_ts = []
    sample_mm_hist = []
    sample_mm_tgt = []
    for s in samples:
        acts = [int(a) for a in s["history_actions"]] + [int(s["target_action"])]
        sample_any.append(int(max(acts) >= 1) if acts else 0)
        sample_rate.append(float(np.mean(acts)) if acts else 0.0)
        sample_ts.append(int(s.get("last_timestamp", 0)))
        hist_cov = np.mean([int(int(i) in mm_rids) for i in s["history_items"]]) if s["history_items"] else 0.0
        sample_mm_hist.append(hist_cov)
        sample_mm_tgt.append(int(int(s["target_item"]) in mm_rids))
    sample_any = np.asarray(sample_any)
    sample_rate = np.asarray(sample_rate)
    sample_ts = np.asarray(sample_ts, dtype=np.int64)
    med_ts = float(np.median(sample_ts)) if sample_ts.size else 0.0
    early = sample_ts < med_ts
    late = sample_ts >= med_ts

    uniq_seq_items = np.unique(item_hits)
    n_uniq = int(uniq_seq_items.size)
    n_mm = int(sum(1 for i in uniq_seq_items if int(i) in mm_rids))

    vc = pd.Series(item_hits).value_counts()
    pop = {
        "n_item_mentions": int(item_hits.size),
        "n_unique_items": n_uniq,
        "top1_share": float(vc.iloc[0] / item_hits.size) if len(vc) else 0.0,
        "top10_share": float(vc.iloc[:10].sum() / item_hits.size) if len(vc) else 0.0,
        "top100_share": float(vc.iloc[:100].sum() / item_hits.size) if len(vc) else 0.0,
        "median_item_freq": float(vc.median()) if len(vc) else 0.0,
        "p95_item_freq": float(vc.quantile(0.95)) if len(vc) else 0.0,
    }

    user_missing = {}
    for c in user_feat.columns:
        if c == "user_id":
            continue
        user_missing[str(c)] = {
            "null_rate": float(user_feat[c].isna().mean()),
            "list_len_mean": float(user_feat[c].map(_len_or_zero).mean()) if str(c) in {"106", "107", "108", "110"} else None,
        }

    item_null = {str(c): float(item_feat[c].isna().mean()) for c in item_feat.columns if c != "item_id"}

    norms = None
    if embs.size:
        rng = np.random.default_rng(0)
        take = rng.choice(embs.shape[0], size=min(5000, embs.shape[0]), replace=False)
        sl = embs[take]
        norms = {
            "n_mm": int(embs.shape[0]),
            "emb_dim": int(embs.shape[1]),
            "head1024_l2_mean": float(np.linalg.norm(sl[:, :1024], axis=1).mean()),
            "tail32_l2_mean": float(np.linalg.norm(sl[:, 1024:], axis=1).mean()),
            "head1024_zero_frac": float((np.linalg.norm(sl[:, :1024], axis=1) < 1e-8).mean()),
            "tail32_zero_frac": float((np.linalg.norm(sl[:, 1024:], axis=1) < 1e-8).mean()),
        }

    probe = {
        "n_seq_rows": int(len(seq_df)),
        "n_samples": int(len(samples)),
        "n_user_feat": int(len(user_feat)),
        "n_item_feat": int(len(item_feat)),
        "n_events": int(n_events),
        "seq_len": {
            "min": int(raw_lens.min()) if raw_lens.size else 0,
            "median": float(np.median(raw_lens)) if raw_lens.size else 0,
            "mean": float(raw_lens.mean()) if raw_lens.size else 0,
            "p05": float(np.quantile(raw_lens, 0.05)) if raw_lens.size else 0,
            "p95": float(np.quantile(raw_lens, 0.95)) if raw_lens.size else 0,
            "max": int(raw_lens.max()) if raw_lens.size else 0,
            "frac_len_100": float((raw_lens == 100).mean()) if raw_lens.size else 0,
        },
        "unique_items_per_user": {
            "mean": float(unique_items.mean()) if unique_items.size else 0,
            "median": float(np.median(unique_items)) if unique_items.size else 0,
        },
        "actions_all_events": _action_counts(all_actions),
        "actions_first_event": _action_counts(first_actions),
        "actions_last_event": _action_counts(last_actions),
        "sample_target_action": _action_counts(sample_tgt),
        "sample_history_len": {
            "mean": float(sample_hist_lens.mean()) if sample_hist_lens.size else 0,
            "median": float(np.median(sample_hist_lens)) if sample_hist_lens.size else 0,
        },
        "any_click_rate": float(sample_any.mean()) if sample_any.size else 0,
        "action_rate_mean": float(sample_rate.mean()) if sample_rate.size else 0,
        "click_position_from_start_mean": float(np.mean(click_positions)) if click_positions else None,
        "position_from_end_click_rate": pos_click_rate[:20],
        "timestamp": {
            "min": int(ts_min_a.min()) if ts_min_a.size else 0,
            "max": int(ts_max_a.max()) if ts_max_a.size else 0,
            "span_sec_mean": float(ts_span_a.mean()) if ts_span_a.size else 0,
            "span_day_mean": float(ts_span_a.mean() / 86400) if ts_span_a.size else 0,
            "span_day_median": float(np.median(ts_span_a) / 86400) if ts_span_a.size else 0,
        },
        "early_vs_late": {
            "median_last_ts": med_ts,
            "n_early": int(early.sum()),
            "n_late": int(late.sum()),
            "any_click_early": float(sample_any[early].mean()) if early.any() else None,
            "any_click_late": float(sample_any[late].mean()) if late.any() else None,
            "action_rate_early": float(sample_rate[early].mean()) if early.any() else None,
            "action_rate_late": float(sample_rate[late].mean()) if late.any() else None,
        },
        "mm_coverage": {
            "unique_seq_items": n_uniq,
            "unique_with_mm": n_mm,
            "frac_unique_with_mm": n_mm / max(n_uniq, 1),
            "hist_item_has_mm_mean": float(np.mean(sample_mm_hist)) if sample_mm_hist else 0,
            "target_has_mm_mean": float(np.mean(sample_mm_tgt)) if sample_mm_tgt else 0,
        },
        "item_popularity": pop,
        "user_feat_missing": user_missing,
        "item_feat_null_rate": item_null,
        "mm_norms_sample": norms,
        "note": (
            "Probe only. Last event is almost always exposure; clicks sit earlier in the window. "
            "Early vs late is a timestamp split, not a treatment."
        ),
    }
    (out_dir / "probe.json").write_text(json.dumps(probe, indent=2))
    pd.DataFrame(pos_click_rate).to_csv(out_dir / "click_rate_by_position_from_end.csv", index=False)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(json.dumps({k: probe[k] for k in ("n_seq_rows", "actions_last_event", "any_click_rate", "mm_coverage", "early_vs_late")}, indent=2))
        print(f"wrote {out_dir / 'probe.json'} (no matplotlib)")
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].hist(raw_lens, bins=30, color="#3b6d9a")
    axes[0, 0].set_title("raw sequence length")
    axes[0, 0].set_xlabel("events")
    ks = [r["k_from_end"] for r in pos_click_rate[:50]]
    rates = [r["click_or_conv_rate"] for r in pos_click_rate[:50]]
    axes[0, 1].plot(ks, rates, color="#b85c38")
    axes[0, 1].set_title("P(click or conversion) by position from end")
    axes[0, 1].set_xlabel("k (1 = last event)")
    axes[0, 1].set_ylabel("rate")
    labels = ["exposure", "click", "conversion"]
    counts = [(all_actions == i).sum() for i in range(3)]
    axes[1, 0].bar(labels, counts, color=["#888", "#3b6d9a", "#b85c38"])
    axes[1, 0].set_title("all events: action_type")
    axes[1, 0].set_yscale("log")
    axes[1, 1].bar(
        ["early", "late"],
        [probe["early_vs_late"]["any_click_early"], probe["early_vs_late"]["any_click_late"]],
        color=["#3b6d9a", "#b85c38"],
    )
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_title("any-click rate, timestamp split")
    fig.tight_layout()
    fig.savefig(out_dir / "probe_overview.png", dpi=140)
    plt.close(fig)

    if len(vc):
        fig, ax = plt.subplots(figsize=(6, 4))
        ranks = np.arange(1, min(len(vc), 5000) + 1)
        ax.loglog(ranks, vc.iloc[: len(ranks)].to_numpy(), color="#3b6d9a")
        ax.set_title("item frequency vs rank")
        ax.set_xlabel("rank")
        ax.set_ylabel("mentions")
        fig.tight_layout()
        fig.savefig(out_dir / "item_zipf.png", dpi=140)
        plt.close(fig)

    print(json.dumps(
        {
            "n_seq_rows": probe["n_seq_rows"],
            "seq_len": probe["seq_len"],
            "actions_all_events": probe["actions_all_events"],
            "actions_last_event": probe["actions_last_event"],
            "any_click_rate": probe["any_click_rate"],
            "early_vs_late": probe["early_vs_late"],
            "mm_coverage": probe["mm_coverage"],
            "item_popularity": probe["item_popularity"],
            "timestamp": probe["timestamp"],
        },
        indent=2,
    ))
    print(f"wrote {out_dir}")


if __name__ == "__main__":
    main()
