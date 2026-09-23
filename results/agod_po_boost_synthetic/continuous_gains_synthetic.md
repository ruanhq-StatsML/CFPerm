# PO-boost synthetic continuous gains (no Affec)

schedule_card=`M_ge_3` mods=`['m0', 'm1', 'm2', 'm3', 'm4']`

| version | FLOPs_rel | T(Acc*) | cumFLOPs@★ | Jaccard | final Acc | ΔAcc | rejects | ship |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| `equal` | 1.000 | 8 | 9.00 | 1.000 | 0.618 | +0.000 | 12 | N |
| `po_soft` | 1.000 | 10 | 11.00 | 1.000 | 0.566 | -0.052 | 12 | N |
| `po_gated` | 1.000 | 10 | 11.00 | 1.000 | 0.566 | -0.052 | 12 | N |
| `po_fuse` | 0.633 | 5 | 4.00 | 0.727 | 0.749 | +0.131 | 12 | Y |

## Reject sources

- `equal`: {'proxy': 11, 'hop_oos': 1}
- `po_soft`: {'proxy': 11, 'hop_oos': 1}
- `po_gated`: {'proxy': 11, 'hop_oos': 1}
- `po_fuse`: {'proxy': 11, 'hop_oos': 1}
