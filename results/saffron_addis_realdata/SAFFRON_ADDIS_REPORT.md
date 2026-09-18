# SAFFRON vs ADDIS on real OnlineRFPerm streams

## Objectives (elaborate)

| Procedure | What it optimizes / controls | Typical behavior |
|---|---|---|
| **SAFFRON** | Online FDR with candidate threshold `λ` (`p≤λ`); wealth spend ×`(1-λ)` | Under conservative nulls can be **under-powered → slower** (larger mean delay / t_first) |
| **ADDIS** | SAFFRON + discard `p>τ`; spend ×`(τ-λ)`; recovers power when nulls are conservative | Higher detection / can look like **higher FAR** vs SAFFRON in theory; here vs α-investing the gap is clearer |
| **α-investing** | Simple wealth rule (no discard / candidate adaptivity) | **Most aggressive** on these streams (highest alarm rate, fastest t_first) |

**Conservative ADDIS knobs:** raise `λ` toward `τ` → shrinks `(τ-λ)` → less α spent per step → **lower alarm rate**, usually **slower**.

Settings: α=0.05, wealth=α/2, SAFFRON λ=0.5, ADDIS τ=0.5, ADDIS λ∈{0.25,0.35,0.40}, burn=8, grace=4, batch=128, seeds=0,1,2.
Datasets: synthetic, electricity, bank, eeg, adult.

## Synthetic (known shift) — delay & pre-shift FAR
```
                 mean_delay  mean_far_pre  mean_AR  mean_early10  mean_t_first  detect
procedure                                                                             
addis_lam0.25          -7.0         0.361    0.250         0.433          13.0     0.0
addis_lam0.35          -7.0         0.361    0.242         0.433          13.0     0.0
addis_lam0.4           -7.0         0.361    0.242         0.433          13.0     0.0
alpha_investing        -7.0         0.389    0.300         0.433          13.0     0.0
saffron                -7.0         0.361    0.242         0.433          13.0     0.0
```
(Note: negative delay = first post-grace reject still before planted shift — cold-start / early alarms.)

## Real data — alarm rate & mean t_first (post-grace)
```
                 mean_AR  mean_early10  mean_t_first  detect
procedure                                                   
addis_lam0.25      0.173         0.267        16.917     1.0
addis_lam0.35      0.171         0.267        16.917     1.0
addis_lam0.4       0.169         0.267        16.917     1.0
alpha_investing    0.248         0.325        14.583     1.0
saffron            0.181         0.267        16.917     1.0
```

## Overall
```
      procedure  mean_alarm_rate  mean_early10  mean_far_pre  mean_t_first  mean_delay  detect_rate  lambda_
  addis_lam0.25         0.188333      0.300000      0.285556     16.133333    5.733333          0.8     0.25
  addis_lam0.35         0.185000      0.300000      0.285556     16.133333    5.733333          0.8     0.35
   addis_lam0.4         0.183333      0.300000      0.285556     16.133333    5.733333          0.8     0.40
alpha_investing         0.258333      0.346667      0.337778     14.266667    3.866667          0.8      NaN
        saffron         0.193333      0.300000      0.285556     16.133333    5.733333          0.8     0.50
```

## Takeaway
- **SAFFRON / ADDIS are slower** than α-investing (~t_first 16.1 vs 14.3 overall; real-data 16.9 vs 14.6).
- **α-investing FAR/AR too high** (~0.25 real); SAFFRON/ADDIS sit ~0.17–0.18.
- Raising ADDIS λ 0.25→0.40: real AR **0.173→0.169** (gentle conservative move via `(τ−λ)`).
- Detection rate on real packs: all ~1.0 after grace; tradeoff shows up in **AR and speed**, not miss rate.
- AR = n_alarm/n_batch baseline under updates — not classical Type-I FAR.
