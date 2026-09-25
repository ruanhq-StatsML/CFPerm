# ChronoBerg AGOD snapshot

Prototype features for Attribution-Guided Online Distillation on
[`spaul25/Chronoberg`](https://huggingface.co/datasets/spaul25/Chronoberg).

- Source: era **test** splits (`1750/1800/1850/1900/1950`) chunked into 80-token passages.
- Modalities (ChronoBerg-native; no audio/image in the corpus):
  - `text`: signed hashing bag-of-words (EmbeddingGemma stand-in)
  - `valence`: period-calibrated VAD valence, feature-hashed
  - `arousal`: period-calibrated VAD arousal, feature-hashed
- Labels for PO-risk: `Y = np.arange(n)` on the stacked reference+current window.

Rebuild the window tensor with:

```bash
python3 scripts/run_chronoberg_agod.py --rebuild-features
```

Raw JSON / lexicons stay local (gitignored). The committed `chronoberg_agod_windows.npz` is the reproducible snapshot.
