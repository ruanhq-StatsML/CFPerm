# RFPerm hallucination prototype

RFPerm **regime fire** = consecutive OOS jump in P(Y|X) for hallucination labels.
Instance **po_risk0** = break of previous probe → candidate under shift.
Not a fact-checker. Not causal.

- cut: 4
- gate: 1.35
- rate before/after: 0.463 / 0.477
- fire_rate: 0.167
- first_fire_t: 4
- hop@cut ranking: `{'n_t1': 100, 'halluc_rate': 0.43, 'auroc_po_risk0': 0.6401468788249696, 'precision_at_10': 0.4}`
- instance probe AUROC: 0.631

Run: `python3 scripts/agod/rfperm_hallucination_proto.py`