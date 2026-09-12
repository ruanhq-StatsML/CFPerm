# AGOD Amazon — modality MSG → per-modality LR

Not inference latency. Online adaptation: attribute text vs image distribution-shift contribution, rescale modality projection LRs, hard-gate inactive modality grads.

Ref category: `Tools & Home Improvement`

- B1: flops=1.000, loss=0.692, α_text=0.500, lr×_text/image=1.00/1.00
- B2: flops=1.000, loss=0.692, α_text=0.481, lr×_text/image=0.97/1.03
- B3: flops=0.500, loss=0.692, α_text=0.299, lr×_text/image=0.66/1.34
