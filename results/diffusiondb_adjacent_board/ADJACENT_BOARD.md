# DiffusionDB adjacent-batch board

Business time cuts → adjacent batch FS → joint board. Not a new estimator; visualization of ΔȲ + late AUC over rolling pairs.

![board](adjacent_board.png)

## width = 30 min
- bins=647, adjacent_pairs=10, reported=0

| t0→t1 | n0/n1 | ΔȲ | HGB AUC | top_fsds |
|---|---|---:|---:|---|

## width = 60 min
- bins=328, adjacent_pairs=10, reported=0

| t0→t1 | n0/n1 | ΔȲ | HGB AUC | top_fsds |
|---|---|---:|---:|---|

## width = 180 min
- bins=110, adjacent_pairs=10, reported=2

| t0→t1 | n0/n1 | ΔȲ | HGB AUC | top_fsds |
|---|---|---:|---:|---|
| 1→2 | 54/47 | 0.0278 | 0.498 | with, by, face |
| 6→7 | 48/43 | 0.0286 | 0.500 | very, 4k, concept |

