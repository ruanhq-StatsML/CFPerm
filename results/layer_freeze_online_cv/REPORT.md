# 从哪一层开始 update

看板直接读 PO-risk。Incoming batch is **T=1**, n_ref=10000, n_new ≥ 5000（不做 online-bootstrap）。

- [electricity (NSW, ordered in time)](electricity/board.html) — start updating from **model_1**
- [covertype (geographic order, class 2 vs rest)](covertype/board.html) — start updating from **model_3**
