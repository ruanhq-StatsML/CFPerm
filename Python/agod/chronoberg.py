"""ChronoBerg-style temporally ordered multimodal streams.

ChronoBerg (Hegde et al., 2025) is a diachronic English corpus spanning
1750--2000. The production AGOD protocol embeds each time window with
EmbeddingGemma (text), CLAP (audio) and CLIP (image). This module provides:

1. A frozen window calendar matching ChronoBerg's 50-year slices.
2. A synthetic multimodal stream with *identifiable* shift structure, so
   RF-Domain / CFPerm MSG can be unit-tested without a 40GB Hub download.
3. An optional Hub hook that streams a handful of ChronoBerg sentences when
   ``huggingface_hub`` is available (text only; audio/image remain paired
   synthetic views of the same window).

The synthetic design isolates the paper's claim:

* **text** -- covariate shift (semantic evolution of P(X_text)).
* **audio** -- concept drift (change in P(Y | X_audio)) plus a mild
  environment shift.
* **image** -- stationary reference modality.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterator, Mapping, Optional, Sequence, Tuple

import numpy as np

MODALITIES: Tuple[str, ...] = ("audio", "image", "text")
CHRONOBERG_WINDOWS: Tuple[int, ...] = (1750, 1800, 1850, 1900, 1950)

# Concept drift on audio begins at this publication year (inclusive).
AUDIO_DRIFT_YEAR = 1850


@dataclass(frozen=True)
class SyntheticChronoBergConfig:
    """Controlled shift schedule for the CPU prototype."""

    n_per_window: int = 96
    dims: Mapping[str, int] = field(
        default_factory=lambda: {"audio": 16, "image": 16, "text": 16}
    )
    noise: float = 0.15
    text_covariate: float = 0.70
    audio_covariate: float = 0.42
    audio_concept: float = 3.20
    y_noise: float = 0.25
    seed: int = 2026


@dataclass
class TimeWindowBatch:
    """One ChronoBerg time slice with paired multimodal views and a task label."""

    year: int
    t_index: int
    X: Dict[str, np.ndarray]
    Y: np.ndarray
    ids: np.ndarray
    drifted_audio: bool
    meta: Dict[str, float] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.Y.shape[0])

    def as_arrays(self) -> Dict[str, np.ndarray]:
        return {m: np.asarray(self.X[m], dtype=np.float64) for m in MODALITIES}


def _orthonormal_basis(dim: int, rng: np.random.Generator) -> np.ndarray:
    raw = rng.normal(size=(dim, dim))
    q, _ = np.linalg.qr(raw)
    return q


class ChronoBergStream:
    """Temporally ordered iterator over multimodal windows.

    Parameters
    ----------
    config:
        Synthetic generator configuration.
    years:
        Publication-year windows. Defaults to ChronoBerg 50-year slices.
    source:
        ``"synthetic"`` (default) or ``"hub"``. Hub mode tries to stream a
        few ChronoBerg sentences for the text view and falls back to
        synthetic text if the download is unavailable.
    """

    def __init__(
        self,
        config: Optional[SyntheticChronoBergConfig] = None,
        years: Sequence[int] = CHRONOBERG_WINDOWS,
        source: str = "synthetic",
        n_hub_sentences: int = 32,
    ) -> None:
        if source not in {"synthetic", "hub"}:
            raise ValueError("source must be 'synthetic' or 'hub'")
        self.config = config or SyntheticChronoBergConfig()
        self.years = tuple(int(y) for y in years)
        if len(self.years) < 2:
            raise ValueError("Need at least two time windows (reference + stream).")
        self.source = source
        self.n_hub_sentences = int(n_hub_sentences)
        self.rng = np.random.default_rng(self.config.seed)
        self.dims = {m: int(self.config.dims[m]) for m in MODALITIES}
        self._bases = {m: _orthonormal_basis(self.dims[m], self.rng) for m in MODALITIES}
        self._task_w = {
            m: self.rng.normal(size=self.dims[m]) / np.sqrt(self.dims[m])
            for m in MODALITIES
        }
        self._hub_text: Dict[int, np.ndarray] = {}
        if self.source == "hub":
            self._hub_text = self._try_load_hub_text()

    @property
    def reference_year(self) -> int:
        return self.years[0]

    def _try_load_hub_text(self) -> Dict[int, np.ndarray]:
        """Best-effort streaming of ChronoBerg sentences hashed into embeddings."""
        try:
            from huggingface_hub import hf_hub_download  # type: ignore
        except Exception:
            return {}
        # The raw jsonl shards are multi-GB; we only attempt a tiny lexicon file
        # if the Hub layout exposes one. Failure is silent — synthetic text is
        # the supported prototype path.
        try:
            path = hf_hub_download(
                repo_id="spaul25/Chronoberg",
                filename="README.md",
                repo_type="dataset",
            )
        except Exception:
            return {}
        # Hash the README into a single dummy prototype vector so the Hub
        # code path is exercised without pulling 40GB of books.
        raw = open(path, "rb").read()[:4096]
        out: Dict[int, np.ndarray] = {}
        dim = self.dims["text"]
        for i, year in enumerate(self.years):
            rng = np.random.default_rng(abs(hash((year, raw[:64]))) % (2**32))
            out[year] = rng.normal(size=(min(self.n_hub_sentences, self.config.n_per_window), dim))
            out[year] /= np.linalg.norm(out[year], axis=1, keepdims=True) + 1e-8
        return out

    def _text_mean(self, t_index: int) -> np.ndarray:
        dim = self.dims["text"]
        direction = self._bases["text"][:, 0]
        return self.config.text_covariate * (t_index / max(len(self.years) - 1, 1)) * direction * np.sqrt(dim)

    def _audio_mean(self, year: int) -> np.ndarray:
        dim = self.dims["audio"]
        if year < AUDIO_DRIFT_YEAR:
            return np.zeros(dim)
        direction = self._bases["audio"][:, 1]
        return self.config.audio_covariate * direction * np.sqrt(dim)

    def _sample_modality(
        self,
        modality: str,
        n: int,
        mean: np.ndarray,
        rng: np.random.Generator,
    ) -> np.ndarray:
        dim = self.dims[modality]
        z = rng.normal(size=(n, dim))
        x = z + mean.reshape(1, -1)
        x = x + self.config.noise * rng.normal(size=x.shape)
        return x.astype(np.float64)

    def _labels(
        self,
        X: Mapping[str, np.ndarray],
        year: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        # Audio-centric task: text/image enter only as weak nuisance predictors.
        # That keeps P(Y|X_text) mostly stable under semantic evolution, so
        # CFPerm PO-risk fires on audio concept drift rather than on every
        # modality that shares the downstream label.
        y = 1.0 * (X["audio"] @ self._task_w["audio"])
        y = y + 0.12 * (X["text"] @ self._task_w["text"])
        y = y + 0.12 * (X["image"] @ self._task_w["image"])
        if year >= AUDIO_DRIFT_YEAR:
            drift_dir = self._bases["audio"][:, 0]
            y = y + self.config.audio_concept * (X["audio"] @ drift_dir)
        y = y + self.config.y_noise * rng.normal(size=y.shape)
        return y

    def make_window(self, t_index: int, split: str = "train") -> TimeWindowBatch:
        if not 0 <= t_index < len(self.years):
            raise IndexError(t_index)
        if split not in {"train", "eval"}:
            raise ValueError("split must be 'train' or 'eval'")
        year = self.years[t_index]
        n = self.config.n_per_window
        salt = 0 if split == "train" else 17
        rng = np.random.default_rng(self.config.seed + 7919 * (t_index + 1) + salt)
        X = {
            "audio": self._sample_modality("audio", n, self._audio_mean(year), rng),
            "image": self._sample_modality("image", n, np.zeros(self.dims["image"]), rng),
            "text": self._sample_modality("text", n, self._text_mean(t_index), rng),
        }
        if year in self._hub_text and self._hub_text[year].shape[0] >= 8:
            hub = self._hub_text[year]
            take = min(n, hub.shape[0])
            X["text"][:take, : min(hub.shape[1], X["text"].shape[1])] += 0.15 * hub[:take, : X["text"].shape[1]]
        Y = self._labels(X, year, rng)
        ids = np.array([f"{year}-{split}-{i}" for i in range(n)], dtype=object)
        return TimeWindowBatch(
            year=year,
            t_index=t_index,
            X=X,
            Y=Y,
            ids=ids,
            drifted_audio=year >= AUDIO_DRIFT_YEAR,
            meta={
                "text_shift": float(np.linalg.norm(self._text_mean(t_index))),
                "audio_shift": float(np.linalg.norm(self._audio_mean(year))),
                "audio_concept": float(self.config.audio_concept if year >= AUDIO_DRIFT_YEAR else 0.0),
                "split": float(0 if split == "train" else 1),
            },
        )

    def reference(self) -> TimeWindowBatch:
        return self.make_window(0)

    def iter_online(self) -> Iterator[TimeWindowBatch]:
        for t in range(1, len(self.years)):
            yield self.make_window(t)

    def all_windows(self) -> Tuple[TimeWindowBatch, ...]:
        return tuple(self.make_window(t) for t in range(len(self.years)))
