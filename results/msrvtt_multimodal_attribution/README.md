# MSR-VTT multimodal FSDS attribution

Window matrix `s`: **768 video + 512 audio + 768 text + 1 label**. Batch $W$ = first vs second half of sliding windows **inside each video**.

| Method | Video | Audio | Text | Diagnostic |
|--------|------:|------:|-----:|------------|
| RF-Domain | 0.772 | 0.191 | 0.037 | AUC 0.804 |
| MMD (coord VIMP) | 0.897 | 0.101 | 0.002 | MMD 0.0000 |
| PO-risk VIMP | 0.813 | 0.174 | 0.013 | R=0.0003 |
| PO permute-group | -0.500 | -0.446 | -0.054 | -- |

Inference: video-clustered bootstrap with **B=10** (mean and variance of replicates), Friedman/Wilcoxon on per-video shares, within-video permutation of $W$ for AUC, and a group-label permutation test that the named 768/512/768 blocks are more imbalanced than random partitions of the same sizes.

n=320 windows, n_videos=16, bootstrap B=10.

## Inference readout

Across RF-Domain, coordinate-MMD, and PO-risk VIMP the **video block dominates** early-vs-late window shift (shares $\approx$ 0.77 / 0.90 / 0.81), then audio, then text.

- Per-video Friedman test of equal RF shares: $p=2.87e-07$ (n=16 videos).
- Wilcoxon signed-rank (Holm) rejects video=audio, video=text, and audio=text at $p<10^{-4}$.
- Group-label permutation (named 768/512/768 vs random partitions): RF $p=0.010$.
- Video-clustered bootstrap $B=10$: RF video$-$audio mean diff $0.635$ (SD $0.039$); the two-sided bootstrap $p$ floor with $B=10$ is $1/11\approx0.091$ (all 10 replicates had the same sign). Use Wilcoxon/Friedman as the primary tests.

## Batch 0 vs Batch 1 region board

See `msrvtt_batch_region_heatmap_board.png` (stats: `msrvtt_region_board_stats.json`). Lead heatmap: windows ordered Batch 0 then Batch 1 × 24 embedding bins (8 video / 8 audio / 8 text). Middle: cosine$(B_0,B_1)$ geometry per modality. Bottom: video × region Cohen's $d$ with Holm stars, plus signed pooled $d$ vs mean $|d|$. Video/audio regions carry a **video-heterogeneous** early-vs-late shift; text is window-invariant ($d=0$). Signed pooled $d$ cancels across videos; mean $|d|$ does not. Mean |Cohen's d|: video $0.787$, audio $0.518$, text $0.000$; 12 Holm-significant video×region cells, 0 pooled-region Holm hits (signed $d$ cancels across videos).

## Per-video AUC

| Video | n | AUC | Video | Audio | Text | Dominant |
|------:|--:|----:|------:|------:|-----:|----------|
| 0 | 20 | 0.990 | 0.952 | 0.048 | 0.000 | video |
| 1 | 20 | 0.845 | 0.681 | 0.319 | 0.000 | video |
| 2 | 20 | 1.000 | 0.934 | 0.066 | 0.000 | video |
| 3 | 20 | 1.000 | 0.959 | 0.041 | 0.000 | video |
| 4 | 20 | 0.940 | 0.904 | 0.096 | 0.000 | video |
| 5 | 20 | 0.780 | 0.746 | 0.254 | 0.000 | video |
| 6 | 20 | 0.870 | 0.895 | 0.105 | 0.000 | video |
| 7 | 20 | 0.910 | 0.705 | 0.295 | 0.000 | video |
| 8 | 20 | 1.000 | 0.968 | 0.032 | 0.000 | video |
| 9 | 20 | 0.980 | 0.838 | 0.162 | 0.000 | video |
| 10 | 20 | 1.000 | 0.743 | 0.257 | 0.000 | video |
| 11 | 20 | 1.000 | 0.967 | 0.033 | 0.000 | video |
| 12 | 20 | 0.980 | 0.842 | 0.158 | 0.000 | video |
| 13 | 20 | 0.870 | 0.484 | 0.516 | 0.000 | audio |
| 14 | 20 | 0.980 | 0.931 | 0.069 | 0.000 | video |
| 15 | 20 | 0.815 | 0.692 | 0.308 | 0.000 | video |
