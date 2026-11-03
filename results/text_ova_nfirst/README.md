# Text multiclass OVA first-n (no tabular)

Text multiclass only — tabular skipped. CE head + sample-level weights; KPI = OVA AUC / acc on first n_first of each new category.

n_first=**100**, seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, weights=**sample-level**.

| dataset | OVA-AUC U | weighted | fused | ΔAUC fused | acc@first U | fused | Δacc |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.769 | 0.705 | 0.705 | -0.064 | 0.402 | 0.407 | 0.005 |
| news20_rec | 0.739 | 0.702 | 0.702 | -0.037 | 0.419 | 0.451 | 0.032 |
| news20_mix | 0.760 | 0.725 | 0.725 | -0.035 | 0.369 | 0.355 | -0.014 |
| cnae9_text | 0.963 | 0.965 | 0.965 | 0.002 | 0.845 | 0.874 | 0.029 |
| amazon_cat | 0.844 | 0.797 | 0.797 | -0.047 | 0.595 | 0.610 | 0.014 |

mean Δ OVA-AUC fused: **-0.036**
mean Δ acc@first fused: **0.013**

## Per-category (last seed, fused vs uniform)

### news20_comp

| cat | n_pos | AUC fused | AUC U | acc fused | acc U |
| --- | --- | --- | --- | --- | --- |
| 5 | 35 | 0.728 | 0.832 | 0.486 | 0.600 |
| 6 | 35 | 0.633 | 0.737 | 0.200 | 0.400 |
| 7 | 35 | 0.729 | 0.674 | 0.343 | 0.257 |
| 8 | 35 | 0.722 | 0.798 | 0.429 | 0.400 |

### news20_rec

| cat | n_pos | AUC fused | AUC U | acc fused | acc U |
| --- | --- | --- | --- | --- | --- |
| 4 | 35 | 0.634 | 0.726 | 0.400 | 0.371 |
| 5 | 35 | 0.843 | 0.843 | 0.657 | 0.600 |
| 6 | 35 | 0.625 | 0.673 | 0.371 | 0.400 |

### news20_mix

| cat | n_pos | AUC fused | AUC U | acc fused | acc U |
| --- | --- | --- | --- | --- | --- |
| 8 | 35 | 0.710 | 0.716 | 0.400 | 0.286 |
| 9 | 35 | 0.635 | 0.739 | 0.229 | 0.314 |
| 10 | 35 | 0.819 | 0.716 | 0.514 | 0.457 |
| 11 | 35 | 0.704 | 0.798 | 0.400 | 0.429 |
| 12 | 35 | 0.780 | 0.848 | 0.457 | 0.486 |
| 13 | 35 | 0.615 | 0.685 | 0.143 | 0.257 |
| 14 | 35 | 0.770 | 0.719 | 0.457 | 0.400 |
| 15 | 35 | 0.725 | 0.745 | 0.171 | 0.371 |

### cnae9_text

| cat | n_pos | AUC fused | AUC U | acc fused | acc U |
| --- | --- | --- | --- | --- | --- |
| 5 | 35 | 0.954 | 0.977 | 0.886 | 0.886 |
| 6 | 35 | 0.931 | 0.968 | 0.857 | 0.743 |
| 7 | 35 | 0.997 | 0.978 | 1.000 | 0.800 |
| 8 | 35 | 0.927 | 0.831 | 0.629 | 0.657 |

### amazon_cat

| cat | n_pos | AUC fused | AUC U | acc fused | acc U |
| --- | --- | --- | --- | --- | --- |
| 2 | 35 | 0.838 | 0.857 | 0.714 | 0.657 |
| 3 | 35 | 0.828 | 0.858 | 0.686 | 0.629 |

```bash
python3 scripts/run_text_ova_nfirst.py --n-first 100 --multiseed
```
