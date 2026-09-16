# Method Landing Roadmap (Inference · Alignment · CS ¥)

CN full text: `docs/biz/METHOD_LANDING_ROADMAP.md`

## One-liner

OnlineRFPerm opens the regime window; axis-split + po_risk0 route actions; **week attribution** makes the ¥ auditable by week (gap≈0); **action payback** turns each arm into an investment with a recovery clock (all same-day in the demo).

## Method innovations (short)

1. **OnlineRFPerm** — online regime probe (not PSI / static AUROC).
2. **Axis split** — style \(P(X)\) vs preference/hallucination \(P(Y|X)\) vs RAG gap; separate tickets.
3. **po_risk0** — instance priority under fire; anneal to uniform when quiet.
4. **Gated adapt** — √PO / action arms only on fire; cost enters payback denominator.
5. **Value law** — only act-vs-ignore volumes × unit prices; no AUROC¥.

## Payback (demo)

| action | payback | net ROI |
|--------|---------|---------|
| audit_topk | 0.009 d | 112.8× |
| retrieval_refresh | 0.017 d | 57.9× |
| model_rollback | 0.033 d | 29.4× |

Fastest payback ≠ largest net ¥: open audit first to stop bleeding, rollback to lock net.

## Week attribution (demo)

| week | gross ¥ | share |
|------|---------|-------|
| 2026-09-14 | ¥7,189 | 43.2% |
| 2026-09-21 | ¥9,467 | 56.8% |

Attribution gap vs ledger: **¥0**.

## Phases

- **Now:** HF protos + CS ledger + week attr + payback + unit-econ band.
- **Next:** ignored-arm opportunity ¥; payback×price joint sensitivity; alignment merge-gate on real candidates.
- **Later:** live serving clock + online FDR; multi-surface split ledgers; auto-disable arms that miss payback.
