# Gradual concept · continuous-time latency

No jump. Last-two hop is the wrong WHEN. Evaluable latency is excess Brier until a detector fires.

Main window n_new=40, onset batch t0=3, n_ref=480.

## Detectors on the main gradual walk

| detector | hat | delay batches | delay obs | α at hat | FA pre | area until hat | never |
| --- | --- | --- | --- | --- | --- | --- | --- |
| last-two hop 1.5× | ∞ | ∞ | ∞ | ∞ | 0 | 0.207 | yes |
| T MA vs pre-onset level | 7 | 4 | 160 | 0.346 | 0 | -0.00836 |  |
| PO MA 2× baseline | ∞ | ∞ | ∞ | ∞ | 0 | 0.207 | yes |
| Brier MA +10% vs pre | 6 | 3 | 120 | 0.269 | 0 | -0.00978 |  |
| π_PO(video)≥0.5 twice | 7 | 4 | 160 | 0.346 | 0 | -0.00836 |  |
| board leaves keep | ∞ | ∞ | ∞ | ∞ | 0 | 0.207 | yes |

## Window size = statistical vs batching latency

| n_new | batches | hop delay obs | T-MA delay obs | Brier+10% delay obs | PO 2× delay obs | hop FA |
| --- | --- | --- | --- | --- | --- | --- |
| 20 | 32 | 220 | 160 | 160 | ∞ | 0 |
| 40 | 16 | ∞ | 160 | 120 | ∞ | 0 |
| 80 | 8 | ∞ | 80 | 80 | ∞ | 0 |

LOGO wall time median 567 ms / window of 40 rows (14.2 ms per row in the window).
