# PO-risk variance — how we characterize it, how we evaluate it

## 1. Unit risk (what varies)

With labels $Y$, unit PO-risk on example $i$ is

$$
\mathrm{PO}_i \;=\; 1 - \hat p(Y_i \mid x_i)
$$

(optional twin: $(Y_i-\hat p_i)^2$). This is a **true-label** risk proxy: high when the detector is wrong or uncertain on the true class.

**Magnitude** (first axis) is the batch mean:

$$
\mu \;=\; \frac1n\sum_{i=1}^n \mathrm{PO}_i.
$$

## 2. Two variance objects (what we mean by “variance of PO-risk”)

We do **not** treat “variance” as a single vague OOD score. Two related estimators:

### A. Sub-environment variance (primary for Mag×Var)

Draw $K$ bootstrap windows $b=1..K$ of size $w=\max(12,\lfloor n/5\rfloor)$ without replacement. Window mean:

$$
r_b \;=\; \frac1{|b|}\sum_{i\in b}\mathrm{PO}_i.
$$

**PO-risk variance (ops headline):**

$$
s \;=\; \mathrm{Std}(r_1,\ldots,r_K),\qquad
v \;=\; \mathrm{Var}(r_1,\ldots,r_K).
$$

Report $\Delta s = s_{\mathrm{LIVE}}-s_{\mathrm{REF}}$ (and $\Delta v$).

**Why this form (justify):**  
Uniform benign shift lifts almost every $r_b$ together → $\mu\uparrow$, $s$ stays quiet.  
Concentrated malignant subset $S$ makes some windows hit $S$ (high $r_b$) and some miss (near REF) → **heteroscedastic** $\{r_b\}$ → $s\uparrow$.  
So $s$ is the variance of *environment-level* risk, which matches the generative model of sparse rule-breaking.

### B. Unit-level variance (confirmatory)

$$
v_{\mathrm{unit}} \;=\; \mathrm{Var}_i(\mathrm{PO}_i).
$$

Useful as a secondary check ($\Delta v_{\mathrm{unit}}$); primary gate still uses **sub-env** $s$, because mixture attacks are defined across environments, not only as a fatter unit histogram.

Derived summaries of the same $\{r_b\}$ (CVaR, max−median) may be logged; they are **not** the characterization headline.

## 3. How we evaluate the variance characterization

Evaluation target is **shift resolution**, not “is variance pretty”:

1. Build paired LIVE with known truth: **benign** (uniform synonym) vs **malignant** (attack subset).  
2. Both are required to sit in a raised-$\mu$ band (so Magnitude alone is uninformative).  
3. Gate on $(\Delta\mu,\,\Delta s)$ (and $\Delta\mathrm{Error}$ when Y exists).  
4. Score with a **2×2 contingency** (true × pred ∈ {benign, malignant}).

**Justification of this eval:** if $s$ really tracks heterogeneity of PO across sub-envs, malignant cells must show systematically larger $\Delta s$ / $J$-gap than benign cells, and the 2×2 must separate them. On five boards (HC3×3 + news20×2) pooled contingency was `[[5,0],[0,5]]`.

Supporting metrics (same panel): `d_s_benign`, `d_s_malignant`, `J_gap`, `malignant_recall`, vs `magnitude_only_acc`.

## 4. One-line justification

> Variance of PO-risk := **dispersion of bootstrap sub-environment mean PO** ($s=\mathrm{Std}(r_b)$).  
> It is justified because sparse malignant mass induces heteroscedastic $\{r_b\}$ while uniform benign mass does not — and that claim is evaluated by a benign/malignant **2×2**, not by Peak aesthetics.
