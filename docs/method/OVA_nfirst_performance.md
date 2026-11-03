# OVA first-n performance (CE head)

Method still uses CE (log_loss). Empowerment locus = new category first-n=100 window; KPI = OVA AUC + acc_on_first; weights vs uniform.

n_first=**100**, seeds=[2026, 2027, 2028], loss=`CE (log_loss)`.

| dataset | OVA-AUC uniform | weighted | fused | ΔAUC fused | acc@first U | fused | Δacc |
| --- | --- | --- | --- | --- | --- | --- | --- |
| har | 0.956 | 0.874 | 0.874 | -0.082 | 0.806 | 0.673 | -0.133 |
| mice | 0.989 | 0.981 | 0.981 | -0.009 | 0.902 | 0.914 | 0.012 |
| cnae9 | 0.959 | 0.966 | 0.966 | 0.007 | 0.826 | 0.833 | 0.007 |
| jvowel | 0.964 | 0.949 | 0.949 | -0.016 | 0.843 | 0.743 | -0.100 |
| news20b | 0.746 | 0.690 | 0.690 | -0.056 | 0.352 | 0.395 | 0.043 |

mean Δ OVA-AUC (fused−uniform): **-0.031**
mean Δ acc@first (fused−uniform): **-0.034**

## Per-category detail (last seed, fused)

### har

| category | n_pos | OVA-AUC fused | OVA-AUC uniform | acc@first fused | acc@first U |
| --- | --- | --- | --- | --- | --- |
| 3 | 35 | 0.564 | 0.865 | 0.057 | 0.257 |
| 4 | 35 | 1.000 | 0.999 | 1.000 | 1.000 |
| 5 | 35 | 0.985 | 1.000 | 0.971 | 1.000 |

### mice

| category | n_pos | OVA-AUC fused | OVA-AUC uniform | acc@first fused | acc@first U |
| --- | --- | --- | --- | --- | --- |
| 4 | 35 | 0.983 | 0.984 | 0.914 | 0.857 |
| 5 | 35 | 0.989 | 1.000 | 0.943 | 1.000 |
| 6 | 35 | 0.957 | 0.979 | 0.857 | 0.971 |
| 7 | 35 | 0.996 | 0.992 | 0.943 | 0.914 |

### cnae9

| category | n_pos | OVA-AUC fused | OVA-AUC uniform | acc@first fused | acc@first U |
| --- | --- | --- | --- | --- | --- |
| 5 | 35 | 0.971 | 0.973 | 0.800 | 0.857 |
| 6 | 35 | 0.953 | 0.972 | 0.771 | 0.829 |
| 7 | 35 | 0.996 | 0.997 | 0.971 | 0.800 |
| 8 | 35 | 0.918 | 0.812 | 0.686 | 0.657 |

### jvowel

| category | n_pos | OVA-AUC fused | OVA-AUC uniform | acc@first fused | acc@first U |
| --- | --- | --- | --- | --- | --- |
| 5 | 35 | 0.999 | 1.000 | 0.971 | 0.943 |
| 6 | 35 | 0.984 | 0.991 | 0.829 | 0.914 |
| 7 | 35 | 0.894 | 0.915 | 0.714 | 0.829 |
| 8 | 35 | 0.975 | 0.924 | 0.486 | 0.743 |

### news20b

| category | n_pos | OVA-AUC fused | OVA-AUC uniform | acc@first fused | acc@first U |
| --- | --- | --- | --- | --- | --- |
| 2 | 35 | 0.744 | 0.802 | 0.514 | 0.429 |
| 3 | 35 | 0.569 | 0.659 | 0.229 | 0.143 |

```bash
python3 scripts/run_ova_nfirst_performance.py --n-first 100 --multiseed
```
