# Board → reason-code generation

Auto-turn sample-chunk adjacent board readings into **routing / claim-control**
codes for review paste. Not convictions; not an HGB ship gate.

## Run

```bash
# after board smoke (also auto-runs inside the board script)
PYTHONPATH=. python3 scripts/generate_board_reason_codes.py \
  --summary results/sample_chunk_adjacent_board/summary.json \
  --out results/sample_chunk_adjacent_board/reason_codes
```

Artifacts: `reason_codes.json`, `REASON_CODES.md`, `paste_for_agent.txt`.

## Logic (rule → code)

| trigger | code | severity |
|---|---|---|
| always | `RC_BOARD_NOT_SHIP` | block_claim |
| always | `RC_SHORTLIST_NOT_DRIVER` | block_claim |
| AUC \< 0.65 | `RC_WEAK_TRANSFER` | info |
| DiffDB + weak | `RC_DIFFDB_TOKEN_NO_TRAVEL` | watch |
| Tencent ops + volume tops + AUC≥0.9 | `RC_INTENSITY_BASELINE` | info |
| content tops credit/share | `RC_CONTENT_CREDIT_CANDIDATE` | investigate |
| ops−content gap ≥ 0.1 | `RC_OPS_CONTENT_GAP` | watch |
| AUC≥0.85 + Jaccard\<0.35 | `RC_SHIFTING_DRIVERS` | watch |
| AUC≥0.85 + Jaccard≥0.5 | `RC_PERSISTENT_ASSOCIATION` | info |
| \|HGB−LogReg\|≥0.15 | `RC_PROBE_DISAGREE` | watch |
| FSDS∩cmean \< 0.15 | `RC_RANK_MEAN_SPLIT` | info |
| Waymo + strong | `RC_PLANTED_SANITY_OK` | info |

Every code sets `allows_ship_model=false` and `allows_driver_claim=false`.

## Mapping to review tips

Volume / credit feature names map through the same copy dictionary as
`export_review_agent_card.TIP_BUCKETS` (末跳、份额、曝光强度…). Content
candidates are meant to land on existing tip buckets; intensity stays ops-only.

## What this is for

- Paste `paste_for_agent.txt` under a review card as structured 案由修饰
- Filter tickets by `code` / `severity`
- Keep DiffDB / ops / content / ship language consistent across packs

## Related

- [`Sample_Chunk_Board_Logic.md`](./Sample_Chunk_Board_Logic.md)
- [`Scope_Review_Agent_Card.md`](./Scope_Review_Agent_Card.md)
