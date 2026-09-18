# Grad-OnlineRFPerm — Sep18 add-on (adult + nomao)

Same single-stream protocol as Sep17 MVP:
`g_t = ||∇_{θ_U} L||_2`, one OnlineRFPerm, seeds `{0..4}`,
batch `128`, `n_batches=48`, `n_burn=8`, `α=0.05`.

| dataset | mean Lead | median | P(Lead<0) | P(Lead≤0) |
|---|---:|---:|---:|---:|
| adult | −3.20 | −3 | 80% | 100% |
| nomao | −4.00 | −4 | 80% | 100% |
| **add-on (10)** | **−3.60** | — | **80%** | **100%** |

Pooled with Sep17 MVP (25+10=35 runs): mean Lead ≈ **−3.34**, P(earlier) ≈ **74%**.

Note: `MagicTelescope` was probed but MSE-OnlineRFPerm never rejects under
this horizon (Grad always does) — kept as Grad-only-detect diagnostic, not
in the lead table.
