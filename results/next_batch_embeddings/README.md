# Next-batch embeddings

Wu 2018 instance discrimination (memory bank), Yu 2020 SDC prototypes, and
Saha 2021 GPM gated by the existing TSS pair `(c_m, δ_m)`.

Not a new TSS loss. The probe loop in `typed_shift_stepsize.py` is unchanged.

## What ran

- InstDisc: linear projector on each FSDS block, nonparametric softmax vs a
  momentum bank; previous-batch keys are extra negatives for the next batch.
  Identification: off-diagonal cosine (collision) rises and own-slot accuracy
  falls on a video location shift.
- SDC: compensate old class means with the current-batch vector field, then
  NCM on the arriving batch. Identification: SDC NCM beats stale NCM under
  covariate-only location shift.
- GPM: SVD bases of `X^{(m)}` are the CGS of that linear head. Next gradient
  is projected iff `c_m` is loud and `δ_m` is quiet (same `τ` as TSS).

## Run

```
python3 scripts/run_next_batch_embeddings.py
python3 scripts/run_next_batch_embeddings.py --quick
```

Method note: `docs/method/Next_batch_embeddings_note.tex`.
Table: `Next_batch_embeddings.tex`.
