# Next-batch group learning rates from π

Not online learning. Freeze π on batch A; finite-step GD on batch B with
η_m = η0 M π_m (mean-preserving). Same simplex as packing / π-RF / adaptive ridge.

| setting | uniform | π-boost | VIMP-boost | π-damp | oracle |
|---|---:|---:|---:|---:|---:|
| synthetic GT=valence | 0.927 | 0.943 | 0.928 | 0.751 | 0.950 |
| inject valence | 0.805 | 0.805 | 0.783 | 0.738 | 0.813 |

Coefficient mass on GT after 18 GD steps:

| setting | uniform | π-boost | VIMP-boost | π-damp |
|---|---:|---:|---:|---:|
| synthetic GT=valence | 0.459 | 0.666 | 0.514 | 0.149 |
| inject valence | 0.268 | 0.325 | 0.242 | 0.143 |

Path plot: `next_batch_lr_path.png`
