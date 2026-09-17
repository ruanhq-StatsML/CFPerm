# Batch size is not when to update — electricity (NSW, ordered in time)

When to update is business logic. Small n_new has no power (500/1000: frac large = 0); the real hop only shows at 2000–5000. That flag is not an update clock.

| n_new | batches | frac large |
|---:|---:|---:|
| 500 | 60 | 0.00 |
| 1000 | 30 | 0.00 |
| 2000 | 15 | 0.13 |
| 5000 | 6 | 0.33 |
