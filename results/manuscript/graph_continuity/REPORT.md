# Graph continuity → OnlineRFPerm

X = GraphSAGE-mean (16-d, frozen) ⊕ hand graph features (7-d). user_id → list in the window → mean-pool. Y = pack-usable rate.

OnlineRFPerm Algorithm 1: frozen RF, $T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}$, rank $p_t$. Last-two `hop_fires` is the adjacent-window gate.

| Regime | y pre | y post | mean T | rank-p < 0.05 | hop@cut | n hop |
|---|---:|---:|---:|---:|---|---:|
| quiet local pack | 0.516 | 0.541 | 0.030 | 0 | no | 0 |
| community hop | 0.516 | 0.392 | 0.049 | 0 | no | 2 |

## Last-two around the cut (community hop)

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.193 | 0.174 | 0.901 |
| 3 | no | 0.174 | 0.242 | 1.391 |
| 4 | no | 0.242 | 0.244 | 1.007 |
| 5 | no | 0.244 | 0.134 | 0.549 |
| 6 | no | 0.134 | 0.082 | 0.616 |

Hotpot file order is not wall-clock time. Figures: `pipeline.png`, `y_rate.png`, `T_and_rank_p.png`. Rebuild: `PYTHONPATH=. python3 scripts/prototype_graph_continuity.py`.
