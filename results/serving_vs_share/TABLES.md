# Serving error vs share error

Hop and continuous time read the same S(τ). The coupling that needs latency is serving rent vs the localization pointer.

t0=3, n_new=40. hat_serve=6, hat_share=7, hat_both=7.
lead-lag (share − serve) = 1 batches.

## Per-batch cells

| t | α | serve excess | share err | π_PO video | serve loud | share loud | cell |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 0 | -0.00789 | 0.555 | 0.555 |  |  | quiet |
| 1 | 0 | 0.0201 | 0 | 0 |  |  | quiet |
| 2 | 0 | -0.0122 | 0.241 | 0.241 |  |  | quiet |
| 3 | 0.0376 | -0.0268 | 0.578 | 0.422 |  |  | quiet |
| 4 | 0.115 | 0.0518 | 0.623 | 0.377 |  |  | quiet |
| 5 | 0.192 | 8.51e-05 | 1 | 0 |  |  | quiet |
| 6 | 0.269 | 0.0592 | 0.0945 | 0.906 | yes |  | serve_only |
| 7 | 0.346 | 0.0906 | 0 | 1 | yes | yes | both |
| 8 | 0.423 | 0.0585 | 0.369 | 0.631 | yes | yes | both |
| 9 | 0.5 | 0.175 | 0.409 | 0.591 | yes | yes | both |
| 10 | 0.577 | 0.134 | 0.269 | 0.731 | yes | yes | both |
| 11 | 0.654 | 0.138 | 0.0628 | 0.937 | yes | yes | both |
| 12 | 0.731 | 0.152 | 0.107 | 0.893 | yes | yes | both |
| 13 | 0.808 | 0.233 | 0.205 | 0.795 | yes | yes | both |
| 14 | 0.885 | 0.254 | 0.245 | 0.755 | yes | yes | both |
| 15 | 0.962 | 0.203 | 0.0943 | 0.906 | yes | yes | both |

## Policy online Brier after onset (vs always-update)

| policy | hat | area vs always |
| --- | --- | --- |
| never | ∞ | 0.207 |
| from_serve | 6 | -0.00836 |
| from_share | 7 | -0.0142 |
| from_both | 7 | -0.0142 |
| always | 0 | 0 |

## Post-onset mean errors by cell

| cell | n | mean serve excess | mean share err |
| --- | --- | --- | --- |
| quiet | 3 | 0.00835 | 0.734 |
| share_only | 0 | ∞ | ∞ |
| serve_only | 1 | 0.0592 | 0.0945 |
| both | 9 | 0.16 | 0.196 |
