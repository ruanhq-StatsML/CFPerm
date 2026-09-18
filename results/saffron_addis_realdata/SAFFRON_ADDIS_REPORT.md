# SAFFRON / ADDIS — manuscript first1/2/3 form

Counting rule (`agod/first_k_metrics.py`, OnlinePermOOB):
```python
first_k = first 1-based end index of k consecutive rejects
P25 / median / P75 = quantiles of *all* alarm times on the trail
DetRate = n_alarm_days / n_days   # NYC-taxi / metro / beijing
```

Settings: α=0.05, wealth=α/2, SAFFRON λ=0.5, ADDIS τ=0.5, ADDIS λ∈[0.25, 0.35, 0.4], burn=8, grace=4, batch=128, seeds=[0, 1, 2].

Datasets: synthetic, covertype, bank, electricity, eeg, adult, metro_interstate, beijing_pm25, nyc_taxi, stocks_AAPL, stocks_MSFT, stocks_IWM, waymo_proxy.

## Stream form — first1 / first2 / first3 + alarm P25/median/P75

| dataset | procedure | first1 | first2 | first3 | P25 | median | P75 | SUM |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| adult | addis_lam0.25 | 13 | --- | --- | 17.2 | 20.5 | 26.2 | 4 |
| adult | addis_lam0.35 | 13 | --- | --- | 17 | 22 | 28.8 | 3.7 |
| adult | addis_lam0.4 | 13 | --- | --- | 17 | 22 | 28.8 | 3.3 |
| adult | alpha_investing | 8 | 17 | --- | 13.8 | 20 | 26.2 | 7 |
| adult | saffron | 13 | --- | --- | 17 | 21 | 28 | 4 |
| bank | addis_lam0.25 | 10 | 18 | 37 | 16 | 19 | 34.5 | 5.7 |
| bank | addis_lam0.35 | 10 | 16.5 | 37 | 16 | 23.5 | 34.5 | 5.3 |
| bank | addis_lam0.4 | 10 | 16.5 | 37 | 16 | 23.5 | 34.5 | 5.3 |
| bank | alpha_investing | 9 | 11 | 37 | 14 | 21 | 32.8 | 8.7 |
| bank | saffron | 10 | 18 | 37 | 16 | 19 | 34.5 | 5.7 |
| beijing_pm25 | addis_lam0.25 | 6 | 7 | --- | 12 | 18 | 31.5 | 7 |
| beijing_pm25 | addis_lam0.35 | 6 | 7 | --- | 12 | 18 | 31.5 | 6.7 |
| beijing_pm25 | addis_lam0.4 | 6 | 7 | --- | 8.2 | 18 | 37 | 5.3 |
| beijing_pm25 | alpha_investing | 6 | 7 | 7 | 12.5 | 18 | 33.5 | 10 |
| beijing_pm25 | saffron | 6 | 7 | --- | 12 | 18 | 31.5 | 7 |
| covertype | addis_lam0.25 | 7 | 6 | --- | 10 | 16 | 21.8 | 8.7 |
| covertype | addis_lam0.35 | 7 | 6 | --- | 10 | 14 | 20 | 8.3 |
| covertype | addis_lam0.4 | 7 | 6 | --- | 10 | 14 | 20 | 8.3 |
| covertype | alpha_investing | 7 | 6 | --- | 10.5 | 16.5 | 24.2 | 9.3 |
| covertype | saffron | 7 | 6 | --- | 10 | 14 | 20 | 8.3 |
| eeg | addis_lam0.25 | 7 | 15 | 26 | 12.5 | 22 | 29 | 5.3 |
| eeg | addis_lam0.35 | 7 | 15 | 26 | 12.5 | 22 | 28 | 5.7 |
| eeg | addis_lam0.4 | 7 | 15 | 26 | 12.5 | 22 | 28 | 5.7 |
| eeg | alpha_investing | 7 | 11 | 26 | 13.2 | 23.5 | 30.2 | 9.3 |
| eeg | saffron | 7 | 15 | 26 | 13.2 | 23 | 29.8 | 6.7 |
| electricity | addis_lam0.25 | 5 | 19 | --- | 7.8 | 17.5 | 31.2 | 5.7 |
| electricity | addis_lam0.35 | 5 | 19 | --- | 7.8 | 17.5 | 31.2 | 5.7 |
| electricity | addis_lam0.4 | 5 | 19 | --- | 7.8 | 17.5 | 31.2 | 5.7 |
| electricity | alpha_investing | 5 | 19 | 33 | 12 | 23 | 31 | 7.7 |
| electricity | saffron | 5 | 19 | --- | 7.8 | 17.5 | 31.2 | 5.7 |
| metro_interstate | addis_lam0.25 | 11 | --- | --- | 13.2 | 16.5 | 19.8 | 2.3 |
| metro_interstate | addis_lam0.35 | 11 | --- | --- | 13.2 | 16.5 | 19.8 | 2 |
| metro_interstate | addis_lam0.4 | 11 | --- | --- | 13.2 | 16.5 | 19.8 | 2 |
| metro_interstate | alpha_investing | 7 | --- | --- | 12.2 | 21 | 29 | 6 |
| metro_interstate | saffron | 11 | --- | --- | 16.5 | 23 | 27 | 3 |
| nyc_taxi | addis_lam0.25 | 9 | 13 | --- | 11.8 | 23.5 | 32.2 | 7.7 |
| nyc_taxi | addis_lam0.35 | 9 | 13 | --- | 11.8 | 23.5 | 32.2 | 7.3 |
| nyc_taxi | addis_lam0.4 | 9 | 13 | --- | 11.5 | 17 | 32.5 | 6.7 |
| nyc_taxi | alpha_investing | 7 | 8 | 9 | 8.8 | 14.5 | 31.5 | 9.7 |
| nyc_taxi | saffron | 9 | 13 | --- | 11.8 | 23.5 | 32.2 | 7.7 |
| stocks_AAPL | addis_lam0.25 | 8 | --- | --- | 8.5 | 19 | 29.5 | 4 |
| stocks_AAPL | addis_lam0.35 | 8 | --- | --- | 8.5 | 19 | 29.5 | 4 |
| stocks_AAPL | addis_lam0.4 | 8 | --- | --- | 8.5 | 19 | 29.5 | 4 |
| stocks_AAPL | alpha_investing | 8 | 30 | --- | 13 | 22 | 28.5 | 7.3 |
| stocks_AAPL | saffron | 8 | --- | --- | 9 | 17 | 29 | 5.3 |
| stocks_IWM | addis_lam0.25 | 7 | --- | --- | 7.5 | 18 | 28.5 | 4 |
| stocks_IWM | addis_lam0.35 | 7 | --- | --- | 7.5 | 18 | 28.5 | 4 |
| stocks_IWM | addis_lam0.4 | 7 | --- | --- | 7.5 | 18 | 28.5 | 4 |
| stocks_IWM | alpha_investing | 7 | --- | --- | 11 | 20 | 25 | 7 |
| stocks_IWM | saffron | 7 | --- | --- | 7.5 | 18 | 28.5 | 4 |
| stocks_MSFT | addis_lam0.25 | 8 | 9 | 10 | 8.2 | 13.5 | 24.8 | 6 |
| stocks_MSFT | addis_lam0.35 | 8 | 9 | 10 | 8.2 | 13.5 | 24.8 | 6 |
| stocks_MSFT | addis_lam0.4 | 9 | 10 | --- | 8.8 | 18 | 28 | 4 |
| stocks_MSFT | alpha_investing | 5 | 9 | 10 | 8 | 18 | 27 | 9 |
| stocks_MSFT | saffron | 8 | 9 | 10 | 8.2 | 13.5 | 24.8 | 6 |
| synthetic | addis_lam0.25 | 5 | 14 | --- | 10.5 | 15 | 21 | 7.7 |
| synthetic | addis_lam0.35 | 5 | 14 | --- | 9.8 | 13.5 | 21 | 7.3 |
| synthetic | addis_lam0.4 | 5 | 14 | --- | 9.8 | 13.5 | 21 | 7.3 |
| synthetic | alpha_investing | 5 | 14 | 14 | 11.8 | 17.5 | 26.8 | 9.7 |
| synthetic | saffron | 5 | 14 | --- | 9.8 | 13.5 | 21 | 7.3 |
| waymo_proxy | addis_lam0.25 | 5 | 9 | 24 | 16.5 | 24.5 | 30.8 | 22 |
| waymo_proxy | addis_lam0.35 | 5 | 9 | 24 | 16.5 | 24.5 | 30.8 | 22 |
| waymo_proxy | addis_lam0.4 | 5 | 9 | 24 | 16.5 | 24.5 | 30.8 | 22 |
| waymo_proxy | alpha_investing | 5 | 9 | 24 | 16.5 | 24.5 | 30.8 | 22 |
| waymo_proxy | saffron | 5 | 9 | 24 | 16.5 | 24.5 | 30.8 | 22 |

## Day-level DetRate (NYC-taxi / metro / beijing)

| dataset | procedure | n_days | n_alarm_days | DetRate | first1 med | first1 P25 | first1 P75 | first1 mean±std |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| metro_interstate | saffron | 300 | 167 | 0.557 | 3 | 3 | 3 | 3.28±0.77 |
| metro_interstate | addis_lam0.25 | 300 | 167 | 0.557 | 3 | 3 | 3 | 3.28±0.77 |
| metro_interstate | addis_lam0.35 | 300 | 167 | 0.557 | 3 | 3 | 3 | 3.28±0.77 |
| metro_interstate | addis_lam0.4 | 300 | 167 | 0.557 | 3 | 3 | 3 | 3.28±0.77 |
| metro_interstate | alpha_investing | 300 | 169 | 0.563 | 3 | 3 | 3 | 3.30±0.78 |
| beijing_pm25 | saffron | 300 | 216 | 0.720 | 3 | 3 | 4 | 3.27±0.44 |
| beijing_pm25 | addis_lam0.25 | 300 | 216 | 0.720 | 3 | 3 | 4 | 3.27±0.44 |
| beijing_pm25 | addis_lam0.35 | 300 | 216 | 0.720 | 3 | 3 | 4 | 3.27±0.44 |
| beijing_pm25 | addis_lam0.4 | 300 | 216 | 0.720 | 3 | 3 | 4 | 3.27±0.44 |
| beijing_pm25 | alpha_investing | 300 | 219 | 0.730 | 3 | 3 | 3 | 3.23±0.42 |
| nyc_taxi | saffron | 215 | 172 | 0.800 | 3 | 3 | 3 | 3.69±1.97 |
| nyc_taxi | addis_lam0.25 | 215 | 172 | 0.800 | 3 | 3 | 3 | 3.69±1.97 |
| nyc_taxi | addis_lam0.35 | 215 | 178 | 0.828 | 3 | 3 | 3 | 3.90±2.25 |
| nyc_taxi | addis_lam0.4 | 215 | 164 | 0.763 | 3 | 3 | 3 | 3.38±1.42 |
| nyc_taxi | alpha_investing | 215 | 212 | 0.986 | 3 | 3 | 7.2 | 4.53±2.51 |

### NYC-taxi headline

- **saffron**: **172/215** days alarmed (DetRate=0.800)
- **addis_lam0.25**: **172/215** days alarmed (DetRate=0.800)
- **addis_lam0.35**: **178/215** days alarmed (DetRate=0.828)
- **addis_lam0.4**: **164/215** days alarmed (DetRate=0.763)
- **alpha_investing**: **212/215** days alarmed (DetRate=0.986)

## Overall (median across datasets×seeds)
```
      procedure  first1_med  first2_med  first3_med  alarm_p25  alarm_median  alarm_p75  mean_SUM  lambda_
  addis_lam0.25         8.0        13.0        17.0      11.75          18.0       28.5  6.923077     0.25
  addis_lam0.35         8.0        13.0        17.0      11.75          18.0       28.5  6.769231     0.35
   addis_lam0.4         8.0        13.0        24.0      11.25          18.0       28.5  6.435897     0.40
alpha_investing         7.0        10.0        10.0      12.25          20.0       29.0  9.435897      NaN
        saffron         7.0        13.0        17.0      11.75          18.0       28.5  7.128205     0.50
```

## Takeaway
- Real-data readout is **first1/2/3 + running alarm-time range (P25/median/P75)**, not overall AR.
- Day packs report **how many calendar days fired** (NYC-taxi DetRate).
- Raising ADDIS λ → more conservative (smaller `(τ−λ)` spend).
