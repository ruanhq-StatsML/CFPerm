# Streaming clock note ($n_{\mathrm{per}}=1$)

Same tidy form and answer-precision $Y$; only the clock changes.

| | default prototype | this run |
|--|--|--|
| $n_{\mathrm{per}}$ | 20 | **1** |
| cut | 5 | 100 |
| $n_{\mathrm{ref}}$ | 100 | 100 |
| detection unit | batch | single turn |

Hotpot mock result here: $y$ quiet→hop = 0→1, but OnlineRFPerm **did not fire**
when each step refits on one previous point (ratio-gate unstable).

**Takeaway:** $n_{\mathrm{per}}=1$ is streaming-shaped and **will** change
behavior. Keep the default batch clock for the main prototype. True streaming
OOB = freeze probe on reference → scalar error trail → `onlinePermOOB` /
`first_k` (batch_size=1), not single-point ORF refits.

See `docs/biz/ONLINERFPERM_RESPONSES_AND_Y.md` §5.
