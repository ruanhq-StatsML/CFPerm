# Prediction table — one Y, many X

Manuscript use-case **4. 大模型审核**. Not attribution. Not VIMP. A supervised table:

- `y`: binary audit decision (pass=1 / fail=0)
- `x_*`: features of the reply (`x_n_toks` … `x_refuse` … `x_thank`)
- `batch`: arrival window
- CFPerm-ready files also have `T` (reviewer queue) and `Y`

HH **chosen/rejected is not Y**. Those fields only say which traffic queue the row came from (helpful vs harmless). The probe predicts the auditor’s pass/fail, not the preference pair.

## Files

| File | n | Y | T / stream |
|---|---:|---|---|
| `xy_hh_helpful_consistent.csv` | 1200 | prototype auditor (policy fixed) | helpful queue |
| `xy_hh_helpful_hop.csv` | 1200 | same X; Y flipped after cut | helpful queue |
| `xy_hh_harmless_consistent.csv` | 1200 | prototype auditor (policy fixed) | harmless queue |
| `xy_hh_harmless_hop.csv` | 1200 | same X; Y flipped after cut | harmless queue |
| `xy_hh_two_stream.csv` | 2400 | auditor Y, not chosen | `T=0` helpful, `T=1` harmless |
| `xy_beavertails.csv` | 1200 | `is_safe` | real safety labels |
| `xy_wildguard.csv` | 1200 | `response_harm_label==unharmful` | real response-harm labels |
| `xy_toxicchat.csv` | 1200 | human `toxicity==0` | real moderation-queue labels |
| `xy_real_two_stream.csv` | 2400 | real labels | `T=0` BeaverTails, `T=1` ToxicChat |

Schema of `xy_*.csv` (except the `*_two_stream.csv` files):

`y,batch,x_n_toks,x_n_chars,x_avg_word,x_qmark,x_bang,x_hedge,x_formal,x_i_count,x_newlines,x_upper,x_refuse,x_please,x_thank`

Online bootstrap + last-two overlay on every label file:

```bash
PYTHONPATH=. python3 scripts/llm_audit_online_bootstrap_prototype.py
```

Numbers: `results/manuscript/llm_audit_online_bootstrap/`. Frozen-ref \(\Delta_t=s_t-\mu_{\mathrm{ref}}\), fire on \(\mathrm{CI}_{\mathrm{lo}}>0\). Last-two `hop_fires` is the other gate.

Rebuild (parquet cache under `data/hf_cache/audit/`, not committed):

```bash
python3 scripts/build_llm_audit_xy.py
```

## How CFPerm reads it

```r
tab <- read.csv("results/manuscript/llm_audit/xy_hh_two_stream.csv")
X <- as.matrix(tab[, grep("^x_", names(tab))])
Y <- tab$Y
T <- tab$T
# T = reviewer queue (helpful vs harmless). Y = pass/fail. Not HH chosen.
```
