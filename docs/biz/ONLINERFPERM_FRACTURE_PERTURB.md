# Fracture: manuscript first1 / first2 (detection delay)

> Smooth quiet → ---/---；中间 fracture → first1=0（delay 0），first2=0（头一回连续两个 trail unit 就火）.

Logic (`docs/biz/FIRST_K_LOGIC.py`):

```python
def first_k_consecutive(det, k):
    det = np.asarray(det, dtype=bool)
    for i in range(len(det) - k + 1):
        if det[i:i+k].all():
            return int(i)  # 0 = first trail unit after cut
    return None
```

| dataset | perturbation | y quiet→frac | smooth first1/first2 | fracture first1/first2 |
|---------|--------------|--------------|----------------------|------------------------|
| halueval | `invent_fracture` | 0.00→1.00 | ---/--- | 0/0 |
| halueval | `label_flip` | 0.00→1.00 | ---/--- | 0/0 |
| halueval | `answer_corrupt` | 0.00→1.00 | ---/--- | 0/0 |
| squad | `invent_fracture` | 0.00→1.00 | ---/--- | 0/0 |
| squad | `label_flip` | 0.00→1.00 | ---/--- | 0/0 |
| squad | `answer_corrupt` | 0.00→1.00 | ---/--- | 0/0 |
| hotpotqa | `invent_fracture` | 0.00→1.00 | ---/--- | 0/0 |
| hotpotqa | `label_flip` | 0.00→1.00 | ---/--- | 0/0 |
| hotpotqa | `answer_corrupt` | 0.00→1.00 | ---/--- | 0/0 |
| truthfulqa | `invent_fracture` | 0.00→1.00 | ---/--- | 0/0 |
| truthfulqa | `label_flip` | 0.00→1.00 | ---/--- | 0/0 |
| truthfulqa | `answer_corrupt` | 0.00→1.00 | ---/--- | 0/0 |
| boolq | `invent_fracture` | 0.00→1.00 | ---/--- | 0/0 |
| boolq | `label_flip` | 0.00→1.00 | ---/--- | 0/0 |
| boolq | `answer_corrupt` | 0.00→1.00 | ---/--- | 0/0 |
| nq_open | `invent_fracture` | 0.00→1.00 | ---/--- | 0/0 |
| nq_open | `label_flip` | 0.00→1.00 | ---/--- | 0/0 |
| nq_open | `answer_corrupt` | 0.00→1.00 | ---/--- | 0/0 |

LaTeX: `docs/biz/ONLINERFPERM_SMOOTH_VS_FIRE.tex`
