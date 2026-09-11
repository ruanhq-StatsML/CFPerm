# MSR-VTT multimodal attribution

- n=5000, fit_n=2000, layout=`video[0:768]|audio[768:1280]|text[1280:2048]|W=labelsmsr`
- RF AUC=0.526; PO-risk=0.000015; MMD²=0.001604

## Modality shares

| method | video | audio | text |
|---|---:|---:|---:|
| RF | 0.453 | 0.159 | 0.388 |
| PO-LOGO | 0.411 | 0.589 | 0.000 |
| MMD-LOGO | 0.336 | 0.493 | 0.172 |
