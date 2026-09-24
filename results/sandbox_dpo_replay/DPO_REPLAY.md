# Flywheel DPO + label pollution + path replay

**Pack:** `metro_interstate` · model=`hgb`

**Reading:** DPO prefers clean path over polluted-label path (margin=510, loss=0.3399)

## Pollution

- method: `np.random.permutation(indices_subset) label shuffle`
- n_polluted indices: 125
- sampled periods: {'n_periods': 3, 'period_len': 50, 'n_indices': 125}

## Rewards

- clean: -502.4932228992061
- polluted: -1522.4945040961675
- replay shuffled path: -502.5199728992061

## DPO

- n_pairs: 2.0
- mean_loss: 0.3399308116028774
- mean_margin: 510.0140155984807
- frac_chosen_better: 1.0

## Decision path replay example

- orig: `['observe', 'forecast:hgb', 'decide:idle', 'log']`
- shuffled: `['observe', 'forecast:hgb', 'decide:idle', 'log']`

Step JSONL under `/workspace/results/sandbox_dpo_replay/metro_interstate/`.

