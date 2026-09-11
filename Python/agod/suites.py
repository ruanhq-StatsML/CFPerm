"""Additional shift suites beyond ChronoBerg.

Each suite exposes the same ChronoBergStream interface so B1–B6 can be
compared on multiple distribution-shift benchmarks:

* ``newsgroups`` — topic windows (sci → rec → talk/politics), hashed
  bag-of-words text with char-3gram ``audio`` (style) and punctuation
  ``image`` views. Concept drift = topic/label map change.
* ``amazon`` — review polarity over years; later windows inject sarcasm
  (label flip given similar lexical features) plus a catalog covariate
  shift (electronics → fashion → grocery).
"""

from __future__ import annotations

import hashlib
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .chronoberg import (
    MODALITIES,
    ChronoBergStream,
    SyntheticChronoBergConfig,
    TimeWindowBatch,
)

# Compact original snippets (not copied from 20 Newsgroups / Amazon dumps).
NEWSGROUP_WINDOWS: Dict[int, List[Tuple[str, int]]] = {
    1990: [  # sci, label 0
        ("orbit Kepler mass gravity telescope spectrum hydrogen fusion", 0),
        ("neuron synapse cortex hippocampus MRI scan blood oxygen", 0),
        ("compiler algorithm complexity polynomial matrix eigenvalue", 0),
        ("climate carbon dioxide glacier ocean current satellite data", 0),
        ("vaccine protein amino acid genome sequence CRISPR lab", 0),
        ("quantum spin entanglement photon detector laser cavity", 0),
        ("fossil sediment isotope strata volcano earthquake plate", 0),
        ("protocol packet latency bandwidth router checksum header", 0),
    ],
    1995: [  # rec, still label 0 (covariate / vocab shift)
        ("biking trail helmet cadence gear mountain weekend camp", 0),
        ("baseball inning pitcher batting average stadium crowd", 0),
        ("guitar amp chord riff studio album tour setlist", 0),
        ("recipe skillet olive garlic simmer roast weekend", 0),
        ("camera lens aperture shutter landscape hike sunset", 0),
        ("chess opening gambit endgame rating tournament clock", 0),
        ("sailing keel wind knot harbor tide weekend race", 0),
        ("yoga stretch breath studio mat weekend class", 0),
    ],
    2000: [  # talk/politics, concept drift: label 1
        ("ballot precinct turnout gerrymander senate floor vote", 1),
        ("tariff trade deficit supply chain factory wage union", 1),
        ("court precedent brief statute jury appeal verdict", 1),
        ("broadcast headline editor bias source leak interview", 1),
        ("zoning housing rent transit commute council hearing", 1),
        ("treaty border asylum visa embassy summit accord", 1),
        ("budget deficit tax credit subsidy ceiling vote", 1),
        ("privacy warrant surveillance metadata subpoena court", 1),
    ],
    2005: [
        ("filibuster cloture amendment markup hearing witness", 1),
        ("sanctions embargo tanker port inspection cargo", 1),
        ("antitrust merger monopoly platform fee store", 1),
        ("spectrum license auction carrier tower rural", 1),
        ("pension unfunded liability actuarial discount rate", 1),
        ("wildfire drought reservoir allocation farm water", 1),
        ("cyber intrusion ransomware agency briefing leak", 1),
        ("primary caucus delegate super Tuesday poll", 1),
    ],
    2010: [
        ("disinfo deepfake platform moderation civic trust", 1),
        ("chip export control foundry lithography subsidy", 1),
        ("grid battery peak load utility rate case", 1),
        ("student debt income share refinance servicer", 1),
        ("antibiotic resistance hospital stewardship trial", 1),
        ("orbital debris launch license mega constellation", 1),
        ("labor strike warehouse contract vote picket", 1),
        ("local news desert nonprofit civic reporting", 1),
    ],
}

AMAZON_WINDOWS: Dict[int, List[Tuple[str, int]]] = {
    2016: [
        ("this laptop battery lasts all day great keyboard", 1),
        ("headphones are comfortable and the mic is clear", 1),
        ("charger failed after a week do not buy", 0),
        ("ssd install was easy boot time is instant", 1),
        ("webcam quality is grainy in evening light", 0),
        ("mouse clicks feel premium and tracking is accurate", 1),
        ("fan noise is loud under light browsing", 0),
        ("monitor colors look accurate out of the box", 1),
    ],
    2018: [
        ("tablet stylus lag is tiny notes look clean", 1),
        ("bluetooth drops in the kitchen constantly", 0),
        ("usb hub powers my drives without a brick", 1),
        ("case yellowed in sunlight within a month", 0),
        ("mechanical keyboard is snappy for coding", 1),
        ("cheap cable broke at the strain relief", 0),
        ("4k stick plays local files without stutter", 1),
        ("router range does not reach the backyard", 0),
    ],
    2020: [  # fashion catalog + sarcasm concept drift
        ("yeah this jacket is 'premium' sure it is", 0),
        ("love how the dye bled onto every shirt", 0),
        ("sneakers actually support a long walk", 1),
        ("sizing is honest and the stitch is tight", 1),
        ("great return policy after the zipper died", 0),
        ("silk scarf looks expensive in daylight", 1),
        ("totally 'durable' after one rain storm", 0),
        ("belt leather smells right and holds shape", 1),
    ],
    2022: [
        ("this coat is a steal if you like holes", 0),
        ("wool sweater stays warm without itching", 1),
        ("hats everywhere except on my head", 0),
        ("denim fade looks intentional and even", 1),
        ("perfect gift if you enjoy sewing buttons back", 0),
        ("canvas tote survived the weekly market", 1),
        ("amazing lining that peeled in a week", 0),
        ("boots broke in fast and grip ice", 1),
    ],
    2024: [
        ("grocery bars taste like the box they came in", 0),
        ("oat milk foams well for morning coffee", 1),
        ("spice kit is fresh and labeled clearly", 1),
        ("yes the 'crisp' chips arrived as dust", 0),
        ("olive oil is peppery and bottled dark", 1),
        ("subscription coffee is always stale", 0),
        ("honey crystallized slowly and tastes floral", 1),
        ("protein powder mixes and is not chalky", 1),
    ],
}


def _hash_token(token: str, dim: int, salt: int) -> int:
    digest = hashlib.blake2b(f"{salt}:{token}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % dim


def hash_bag(texts: Sequence[str], dim: int, salt: int, dropout: float = 0.0, rng: np.random.Generator | None = None) -> np.ndarray:
    x = np.zeros((len(texts), dim), dtype=np.float64)
    for i, text in enumerate(texts):
        toks = [t for t in text.lower().replace("'", " ").split() if t]
        if dropout > 0 and rng is not None and toks:
            keep = [t for t in toks if rng.random() > dropout]
            toks = keep or toks[:1]
        for tok in toks:
            x[i, _hash_token(tok, dim, salt)] += 1.0
        n = max(len(toks), 1)
        x[i] /= n
    return x


def char_trigram(texts: Sequence[str], dim: int, salt: int) -> np.ndarray:
    x = np.zeros((len(texts), dim), dtype=np.float64)
    for i, text in enumerate(texts):
        s = f"  {text.lower()}  "
        for j in range(len(s) - 2):
            x[i, _hash_token(s[j : j + 3], dim, salt)] += 1.0
        x[i] /= max(len(s) - 2, 1)
    return x


def punct_image(texts: Sequence[str], dim: int, rng: np.random.Generator) -> np.ndarray:
    feats = []
    for text in texts:
        n = max(len(text), 1)
        feats.append(
            [
                text.count(" ") / n,
                text.count("'") / n,
                text.count(",") / n,
                sum(c.isupper() for c in text) / n,
                np.log1p(n) / 6.0,
            ]
        )
    base = np.asarray(feats, dtype=np.float64)
    proj = rng.normal(size=(base.shape[1], dim))
    proj /= np.linalg.norm(proj, axis=0, keepdims=True) + 1e-8
    return base @ proj


def _tile(pairs: Sequence[Tuple[str, int]], n: int, rng: np.random.Generator) -> List[Tuple[str, int]]:
    idx = rng.integers(0, len(pairs), size=n)
    out = []
    for i in idx:
        text, y = pairs[int(i)]
        jitter = rng.choice(text.split())
        out.append((text + " " + jitter, int(y)))
    return out


class TextShiftStream(ChronoBergStream):
    """ChronoBerg-compatible stream backed by hashed text windows."""

    def __init__(
        self,
        windows: Dict[int, List[Tuple[str, int]]],
        *,
        name: str,
        drift_year: int,
        config: SyntheticChronoBergConfig | None = None,
        n_per_window: int = 64,
        dim: int = 14,
        seed: int = 2026,
    ) -> None:
        cfg = config or SyntheticChronoBergConfig(
            n_per_window=n_per_window,
            dims={"audio": dim, "image": dim, "text": dim},
            seed=seed,
        )
        years = tuple(sorted(windows.keys()))
        super().__init__(config=cfg, years=years, source="synthetic")
        self.suite_name = name
        self.drift_year = int(drift_year)
        self._pairs = windows
        self._dim = dim

    def make_window(self, t_index: int, split: str = "train") -> TimeWindowBatch:
        if not 0 <= t_index < len(self.years):
            raise IndexError(t_index)
        year = self.years[t_index]
        salt = 0 if split == "train" else 17
        rng = np.random.default_rng(self.config.seed + 7919 * (t_index + 1) + salt)
        n = self.config.n_per_window
        pairs = _tile(self._pairs[year], n, rng)
        texts = [p[0] for p in pairs]
        y = np.array([p[1] for p in pairs], dtype=np.float64)
        y = 2.0 * y - 1.0
        x_text = hash_bag(texts, self._dim, salt=11)
        x_audio = char_trigram(texts, self._dim, salt=23)
        x_image = punct_image(texts, self._dim, rng)
        # Mild extra covariate on text as years move, matching ChronoBerg MSG.
        x_text = x_text + 0.15 * t_index * rng.normal(size=x_text.shape) / self._dim
        return TimeWindowBatch(
            year=year,
            t_index=t_index,
            X={"audio": x_audio, "image": x_image, "text": x_text},
            Y=y,
            ids=np.array([f"{self.suite_name}-{year}-{split}-{i}" for i in range(n)], dtype=object),
            drifted_audio=year >= self.drift_year,
            meta={"suite": float(hash(self.suite_name) % 1000), "split": float(0 if split == "train" else 1)},
        )


def make_suite(name: str, n_per_window: int = 64, seed: int = 2026, dim: int = 14) -> ChronoBergStream:
    name = name.lower()
    if name in {"chronoberg", "chrono"}:
        return ChronoBergStream(
            SyntheticChronoBergConfig(
                n_per_window=n_per_window,
                dims={"audio": dim, "image": dim, "text": dim},
                seed=seed,
            )
        )
    if name in {"newsgroups", "news", "20ng"}:
        return TextShiftStream(
            NEWSGROUP_WINDOWS,
            name="newsgroups",
            drift_year=2000,
            n_per_window=n_per_window,
            dim=dim,
            seed=seed,
        )
    if name in {"amazon", "reviews"}:
        return TextShiftStream(
            AMAZON_WINDOWS,
            name="amazon",
            drift_year=2020,
            n_per_window=n_per_window,
            dim=dim,
            seed=seed,
        )
    raise ValueError(f"Unknown suite {name!r}; expected chronoberg, newsgroups, amazon.")


SUITE_NAMES: Tuple[str, ...] = ("chronoberg", "newsgroups", "amazon")
