"""B1 / B2 / B3 comparison on a ChronoBerg-style stream."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .chronoberg import ChronoBergStream, SyntheticChronoBergConfig
from .metrics import attribution_consistency, forgetting_rate
from .online import AGODConfig, AGODTrainer, run_stream


BASELINES: Sequence[str] = ("B1", "B2", "B3")


@dataclass
class ExperimentResult:
    metrics: Dict[str, Dict[str, float]]
    trajectories: Dict[str, List[Dict[str, object]]]
    ranking: List[str]


def _summarize(trainer: AGODTrainer) -> Dict[str, float]:
    logs = trainer.logs
    if not logs:
        return {}
    drift_logs = [lg for lg in logs if lg.drifted_audio] or list(logs)
    stable_logs = [lg for lg in logs if not lg.drifted_audio]
    audio_alpha = float(np.mean([lg.alpha["audio"] for lg in drift_logs])) if drift_logs else 0.0
    text_alpha = float(np.mean([lg.alpha["text"] for lg in drift_logs])) if drift_logs else 0.0
    drift_recall = float(np.mean([lg.drift_recall for lg in drift_logs])) if drift_logs else 0.0
    mean_recall = float(np.mean([lg.recall for lg in logs]))
    audio_align = float(np.mean([lg.alignment["audio"] for lg in drift_logs])) if drift_logs else 0.0
    image_align_early = float(np.mean([lg.alignment["image"] for lg in stable_logs])) if stable_logs else 0.0
    image_align_late = float(np.mean([lg.alignment["image"] for lg in drift_logs])) if drift_logs else 0.0
    true_shift = {"audio": 2.0, "image": 0.0, "text": 1.0}
    attr = float(
        np.mean(
            [
                attribution_consistency(lg.gap, true_shift)
                for lg in drift_logs
            ]
        )
    ) if drift_logs else 0.0
    flops = float(np.sum([lg.distill_flops for lg in logs]))
    return {
        "drift_subgroup_recall": drift_recall,
        "mean_recall": mean_recall,
        "audio_alpha_on_drift": audio_alpha,
        "text_alpha_on_drift": text_alpha,
        "audio_alignment_on_drift": audio_align,
        "image_forgetting": forgetting_rate([image_align_early], [image_align_late]) if stable_logs and drift_logs else 0.0,
        "attribution_consistency": attr,
        "mean_loss": float(np.mean([lg.loss for lg in logs])),
        "steps": float(len(logs)),
        "distill_flops": flops,
        "skip_rate": float(np.mean([lg.skipped for lg in logs])),
    }


def _trajectory(trainer: AGODTrainer) -> List[Dict[str, object]]:
    rows = []
    for lg in trainer.logs:
        rows.append(
            {
                "year": lg.year,
                "baseline": lg.baseline,
                "alpha": lg.alpha,
                "auc": lg.auc,
                "po_risk": lg.po_risk,
                "gap": lg.gap,
                "recall": lg.recall,
                "drift_recall": lg.drift_recall,
                "alignment": lg.alignment,
                "loss": lg.loss,
                "localize": list(lg.localize),
                "distill_flops": lg.distill_flops,
                "skipped": lg.skipped,
            }
        )
    return rows


def run_baseline_comparison(
    config: Optional[SyntheticChronoBergConfig] = None,
    trainer_config: Optional[AGODConfig] = None,
    baselines: Sequence[str] = BASELINES,
    stream_factory=None,
) -> ExperimentResult:
    config = config or SyntheticChronoBergConfig()
    trainer_config = trainer_config or AGODConfig(seed=config.seed)
    metrics: Dict[str, Dict[str, float]] = {}
    trajectories: Dict[str, List[Dict[str, object]]] = {}
    for name in baselines:
        stream = stream_factory() if stream_factory is not None else ChronoBergStream(config=config)
        trainer = run_stream(stream, baseline=name, config=trainer_config)
        metrics[name] = _summarize(trainer)
        trajectories[name] = _trajectory(trainer)
    b1_flops = metrics.get("B1", {}).get("distill_flops", 0.0) or 1.0
    for row in metrics.values():
        row["rel_flops"] = float(row.get("distill_flops", 0.0) / b1_flops)
    ranking = sorted(
        baselines,
        key=lambda b: (
            metrics[b].get("drift_subgroup_recall", 0.0),
            metrics[b].get("audio_alignment_on_drift", 0.0),
            -metrics[b].get("rel_flops", 1.0),
        ),
        reverse=True,
    )
    return ExperimentResult(metrics=metrics, trajectories=trajectories, ranking=ranking)


def save_result(result: ExperimentResult, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "metrics": result.metrics,
        "ranking": result.ranking,
        "trajectories": result.trajectories,
    }
    path.write_text(json.dumps(payload, indent=2, default=float))


def plot_result(result: ExperimentResult, out_dir: Path) -> Dict[str, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}

    fig, ax = plt.subplots(figsize=(7.6, 4.3))
    names = list(result.metrics)
    x = np.arange(len(names))
    width = 0.36
    rec = [result.metrics[n]["drift_subgroup_recall"] for n in names]
    align = [result.metrics[n]["audio_alignment_on_drift"] for n in names]
    ax.bar(x - width / 2, rec, width, label="Audio Recall@5", color="#c0392b")
    ax.bar(x + width / 2, align, width, label="Audio cosine alignment", color="#3b6ea8")
    ax.set_xticks(x, names)
    ax.set_ylabel("Held-out score")
    ax.set_title("AGOD vs static and covariate-only distillation")
    ax.set_ylim(0.0, max(0.05, max(rec + align) * 1.25))
    ax.legend(frameon=False)
    fig.tight_layout()
    paths["recall"] = out_dir / "drift_recall.png"
    fig.savefig(paths["recall"], dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    for name, traj in result.trajectories.items():
        years = [row["year"] for row in traj]
        audio = [row["alpha"]["audio"] for row in traj]
        ax.plot(years, audio, marker="o", label=f"{name} α_audio")
    ax.axvline(1850, color="k", ls="--", lw=0.8, label="audio concept drift")
    ax.set_xlabel("ChronoBerg window (year)")
    ax.set_ylabel("Audio distillation weight")
    ax.set_title("Routing capacity toward the drifting modality")
    ax.legend(frameon=False)
    ax.set_ylim(0.0, 1.0)
    fig.tight_layout()
    paths["routing"] = out_dir / "audio_routing.png"
    fig.savefig(paths["routing"], dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6), sharey=False)
    b3 = result.trajectories.get("B3", [])
    years = [row["year"] for row in b3]
    for ax, key, title in zip(
        axes,
        ("auc", "po_risk", "gap"),
        ("RF-Domain AUC", "CFPerm PO-risk", "MSG gap"),
    ):
        for m in ("audio", "image", "text"):
            ax.plot(years, [row[key][m] for row in b3], marker="o", label=m)
        ax.set_title(title)
        ax.set_xlabel("year")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Modality-specific gap decomposition (AGOD / B3)", y=1.03)
    fig.tight_layout()
    paths["msg"] = out_dir / "msg_decomposition.png"
    fig.savefig(paths["msg"], dpi=140, bbox_inches="tight")
    plt.close(fig)
    return paths
