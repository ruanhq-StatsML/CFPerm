"""LaTeX performance tables from board summaries."""
from __future__ import annotations

from agod.performance_tex import real_tables, synth_tables, write_all_tex


def _synth():
    hop = {str(t): 1.0 + 0.1 * t for t in range(1, 7)}
    hop["3"] = 5.0
    hop["4"] = 2.0
    methods = {
        name: {"mse_mean": mse, "fire_mean": fire, "by_t": dict(hop)}
        for name, mse, fire in (
            ("uniform_pair", 1.5, 0.0),
            ("rfperm", 1.4, 0.17),
            ("resid", 1.3, 0.17),
            ("oracle", 1.35, 1.0),
        )
    }
    return {
        "rf": {
            "similar": {
                "methods": {
                    k: {"mse_mean": 0.54, "fire_mean": 0.0, "by_t": hop}
                    for k in methods
                }
            },
            "concept": {"methods": methods},
        }
    }


def test_synth_tables_have_fire_and_hop_path():
    tex = "\n".join(synth_tables(_synth(), 1.5, ["rf"]))
    assert r"\begin{tabular}" in tex
    assert "similar" in tex
    assert "fire" in tex
    assert r"$\leftarrow$ cut" in tex
    assert r"\mathbf{1.300}" in tex or r"\mathbf{1.3}" in tex


def test_real_tables_rmse_and_acc():
    rows = []
    for method, mse, acc, fire in (
        ("uniform_pair", 100.0, 0.80, 0.0),
        ("rfperm", 81.0, 0.82, 0.05),
        ("resid", 121.0, 0.79, 0.10),
    ):
        rows.append(
            {
                "dataset": "nyc_taxi",
                "clock": "time",
                "task": "mse",
                "learner": "rf",
                "method": method,
                "online_score": mse,
                "fire_rate": fire,
                "history": [],
            }
        )
        rows.append(
            {
                "dataset": "occupancy",
                "clock": "time",
                "task": "acc",
                "learner": "rf",
                "method": method,
                "online_score": acc,
                "fire_rate": fire,
                "history": [
                    {
                        "t": 5,
                        "score": acc,
                        "fired": method == "rfperm",
                    }
                ],
            }
        )
    tex = "\n".join(real_tables(rows, 1.5, 200))
    assert "nyc\\_taxi" in tex
    assert "RMSE" in tex
    assert "Acc" in tex
    assert r"\mathbf{9.000}" in tex or r"\mathbf{9.00}" in tex
    assert r"\mathbf{0.820}" in tex


def test_write_all_tex(tmp_path):
    path = write_all_tex(
        {"summary": _synth(), "gate": 1.5, "learners": ["rf"]},
        {"rows": [], "gate": 1.5, "n_per": 200},
        tmp_path / "all.tex",
    )
    text = path.read_text(encoding="utf-8")
    assert "tab:po-refit-mse-rf" in text
