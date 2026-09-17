#!/usr/bin/env python3
"""Layer-freeze online CV: k+1 models, T=1 streaming PO-risk."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from dl_model_registry import (  # noqa: E402
    AnyMLP,
    DLModelRegistry,
    TabularDataset,
    apply_train_top_i,
    construct_dataloader,
    spawn_layer_models,
)
from layer_freeze_cv import (  # noqa: E402
    attach_layer_dicts,
    run_deviation_gate,
    run_layer_freeze_cv,
    stack_layer_metric_dicts,
)
from streaming_po_risk import (  # noqa: E402
    ACTION_FREEZE,
    ACTION_KEEP,
    ACTION_WATCH,
    MIN_STREAM_N,
    REF_N,
    TabularPORisk,
    annotate_moving_average,
    annotate_po_mse_contrast,
    batch_mse,
    large_deviation,
    ma_window,
    moving_average,
    pack_ref_new,
    po_mse_action,
    po_risk,
    ref_split_baseline,
    streaming_po_and_mse,
    streaming_po_risk,
)


class RegistryTests(unittest.TestCase):
    def test_tabular_dataset_returns_row_and_label(self):
        ds = TabularDataset([[1.0, 2.0], [3.0, 4.0]], [0, 1])
        x, y = ds[1]
        self.assertEqual(tuple(x.shape), (2,))
        self.assertEqual(int(y.item()), 1)

    def test_any_mlp_k_plus_one_models_train_top_i(self):
        m = AnyMLP(4, 1, (8, 4), dropout=0.0)
        self.assertEqual(m.n_layer_groups, 3)
        clones = spawn_layer_models(m)
        self.assertEqual(len(clones), 4)
        self.assertEqual([c._registry_name for c in clones], ["model_0", "model_1", "model_2", "model_3"])
        n_train = []
        for c in clones:
            n_train.append(sum(p.requires_grad for p in c.parameters()))
        self.assertEqual(n_train[0], 0)
        self.assertGreater(n_train[1], 0)
        self.assertGreater(n_train[-1], n_train[1])
        apply_train_top_i(clones[1], 1)
        head_ok = all(p.requires_grad for p in clones[1].head.parameters())
        bottom_frozen = all(not p.requires_grad for p in clones[1].blocks["fc0"].parameters())
        self.assertTrue(head_ok)
        self.assertTrue(bottom_frozen)

    def test_forward_logits_shape(self):
        m = AnyMLP(5, 1, (6,), dropout=0.0)
        y = m(torch.randn(7, 5))
        self.assertEqual(tuple(y.shape), (7, 1))


class StreamingPOTests(unittest.TestCase):
    def test_ref_and_stream_are_large_enough_to_read_po_risk(self):
        self.assertEqual(REF_N, 10_000)
        self.assertGreaterEqual(MIN_STREAM_N, 5_000)

    def test_new_batch_is_t1(self):
        X0 = np.zeros((10, 2))
        X1 = np.ones((6, 2))
        _, _, T = pack_ref_new(X0, np.zeros(10), X1, np.ones(6))
        self.assertTrue(np.all(T[:10] == 0))
        self.assertTrue(np.all(T[10:] == 1))

    def test_po_risk_fires_when_y_hops_with_t(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(400, 4))
        T = np.array([0] * 200 + [1] * 200)
        Y = (X[:, 0] + 1.8 * T > 0).astype(float)
        hopped = po_risk(X, Y, T)
        Y_null = (X[:, 0] > 0).astype(float)
        quiet = po_risk(X, Y_null, T)
        self.assertGreater(hopped, quiet)

    def test_streaming_wrapper_matches_pack(self):
        rng = np.random.default_rng(1)
        X0, Y0 = rng.normal(size=(80, 3)), rng.integers(0, 2, size=80).astype(float)
        X1, Y1 = rng.normal(size=(80, 3)) + 0.5, rng.integers(0, 2, size=80).astype(float)
        a = streaming_po_risk(X0, Y0, X1, Y1, seed=1)
        X, Y, T = pack_ref_new(X0, Y0, X1, Y1)
        b = po_risk(X, Y, T, seed=1)
        self.assertAlmostEqual(a, b, places=12)

    def test_outcome_and_propensity_are_separate_models(self):
        rng = np.random.default_rng(4)
        X = rng.normal(size=(240, 4))
        T = np.array([0] * 120 + [1] * 120)
        Y = (X[:, 0] + 1.5 * T > 0).astype(float)
        est = TabularPORisk(seed=4)
        out = est.risk(X, Y, T)
        self.assertIsNotNone(est.outcome)
        self.assertIsNotNone(est.propensity)
        self.assertGreater(out["po_risk"], 0.0)
        self.assertEqual(len(out["mu"]), 240)
        self.assertEqual(len(out["e"]), 240)

    def test_large_deviation_reads_stream_against_ref_split(self):
        rng = np.random.default_rng(5)
        X0 = rng.normal(size=(200, 3))
        Y0 = (X0[:, 0] > 0).astype(float)
        base = ref_split_baseline(X0, Y0, seed=5)
        X1 = rng.normal(size=(200, 3))
        Y1 = (1.0 - (X1[:, 0] > 0).astype(float))
        hopped = streaming_po_risk(X0, Y0, X1, Y1, seed=5)
        self.assertTrue(large_deviation(hopped, base))
        quiet = streaming_po_risk(X0[:100], Y0[:100], X0[100:], Y0[100:], seed=5)
        self.assertFalse(large_deviation(base, base))
        self.assertTrue(large_deviation(2.1 * (base + 1e-12), base + 1e-12))
        self.assertGreaterEqual(hopped, quiet)

    def test_moving_average_is_smoother_and_preserves_level(self):
        rng = np.random.default_rng(0)
        y = 0.01 + 0.2 * rng.normal(size=80)
        ma = moving_average(y, 10)
        self.assertEqual(len(ma), len(y))
        self.assertLess(float(ma.std()), float(y.std()))
        self.assertAlmostEqual(ma[0], y[0])
        self.assertAlmostEqual(ma[9], float(y[:10].mean()))
        self.assertAlmostEqual(ma[10], float(y[1:11].mean()))
        self.assertEqual(len(moving_average([], 5)), 0)

    def test_ma_stable_means_all_layer_backprop(self):
        self.assertEqual(ma_window(20), 50)
        baseline = 1e-6
        quiet = [{"t": t, "n_new": 20, "po_stream": baseline * 0.4} for t in range(40)]
        info = annotate_moving_average(quiet, baseline, n_new=20)
        self.assertTrue(info["ma_stable"])
        self.assertTrue(info["all_layer_backprop"])
        self.assertLess(info["ma_max"], 2.0 * baseline)
        hopped = [{"t": t, "n_new": 20, "po_stream": baseline * 3.0} for t in range(40)]
        info2 = annotate_moving_average(hopped, baseline, n_new=20)
        self.assertFalse(info2["ma_stable"])
        self.assertFalse(info2["all_layer_backprop"])


class FreezeCvSmokeTests(unittest.TestCase):
    def test_model_0_weights_do_not_move(self):
        rng = np.random.default_rng(2)
        X = rng.normal(size=(80, 3))
        Y = (X[:, 0] > 0).astype(float)
        registry = DLModelRegistry(epochs=1, patience=2, amp=False, lr=1e-2)
        adapter = registry.make_any_mlp(input_dim=3, output_dim=1, hidden_dims=(6,))
        model = adapter.build()
        apply_train_top_i(model, 0)
        before = {k: v.clone() for k, v in model.state_dict().items()}
        loader = construct_dataloader(X, Y, batch_size=32, shuffle=False)
        registry._train_loop(model, loader, epochs=2, restore_best=False)
        for k, v in model.state_dict().items():
            self.assertTrue(torch.allclose(before[k], v))

    def test_cv_gate_then_rows(self):
        rng = np.random.default_rng(3)
        n_ref, n_new = 220, 180
        X0 = rng.normal(size=(n_ref, 4))
        Y0 = (X0[:, 0] > 0).astype(float)
        X1 = rng.normal(size=(n_new, 4))
        Y1 = (X1[:, 0] + 1.2 > 0).astype(float)
        X = np.vstack([X0, X1])
        Y = np.concatenate([Y0, Y1])
        registry = DLModelRegistry(epochs=2, patience=2, amp=False, lr=5e-3, verbose=False)
        out = run_layer_freeze_cv(
            X,
            Y,
            n_ref=n_ref,
            batch_size_stream=90,
            hidden_dims=(8, 4),
            online_epochs=1,
            max_batches=2,
            registry=registry,
            n_ref_eval=80,
        )
        self.assertEqual(out["k"], 3)
        self.assertEqual(len(out["layer_names"]), 4)
        self.assertEqual(out["n_batches"], 2)
        self.assertEqual(len(out["rows"]), 2)
        self.assertEqual(set(out["MSE_Dict"]), {"layer0", "layer1", "layer2", "layer3"})
        self.assertEqual(len(out["PO_Dict"]["layer0"]), 2)
        self.assertTrue(all(r["n_new"] == 90 for r in out["rows"]))
        self.assertIn("large_deviation", out["rows"][0])
        self.assertIn("po_base", out["rows"][0])
        for r in out["rows"]:
            self.assertIn(r["action"], {ACTION_KEEP, ACTION_WATCH, ACTION_FREEZE})
            self.assertIn("mse_stream", r)
            if r["action"] == ACTION_FREEZE:
                self.assertEqual(len(r["layers"]), 4)
                self.assertIn(r["i_star"], {0, 1, 2, 3})
            else:
                self.assertTrue(r["all_trainable"])
                self.assertIsNone(r["freeze_from"])

    def test_deviation_gate_does_not_claim_when_to_update(self):
        rng = np.random.default_rng(6)
        X = rng.normal(size=(400, 3))
        Y = (X[:, 0] > 0).astype(float)
        out = run_deviation_gate(X, Y, n_ref=200, batch_size_stream=50, max_batches=4)
        self.assertEqual(out["n_batches"], 4)
        self.assertEqual(out["n_new"], 50)
        self.assertIn("frac_large", out)
        self.assertIn("all_layer_backprop", out)
        self.assertIn("ma_window", out)
        self.assertIn("po_ma", out["rows"][0])

    def test_n_new_20_keeps_n_ref_and_skips_bootstrap(self):
        rng = np.random.default_rng(8)
        X = rng.normal(size=(260, 4))
        Y = (X[:, 0] > 0).astype(float)
        out = run_deviation_gate(X, Y, n_ref=200, batch_size_stream=20, max_batches=3)
        self.assertEqual(out["n_ref"], 200)
        self.assertEqual(out["n_new"], 20)
        self.assertEqual(out["n_batches"], 3)
        self.assertTrue(all(r["n_new"] == 20 for r in out["rows"]))
        self.assertNotIn("bootstrap", out)
        self.assertTrue(all("po_ma" in r for r in out["rows"]))

    def test_po_and_mse_dicts_align_to_stream(self):
        rows = [
            {"t": 0, "n_new": 20, "po_stream": 1e-6, "layers": []},
            {
                "t": 1,
                "n_new": 20,
                "po_stream": 3e-6,
                "layers": [
                    {"i": 0, "po_fit": 1.0, "mse": 0.40},
                    {"i": 1, "po_fit": 0.8, "mse": 0.30},
                    {"i": 2, "po_fit": 0.5, "mse": 0.22},
                ],
            },
        ]
        stacked = stack_layer_metric_dicts(rows, k=2)
        self.assertEqual(list(stacked["MSE_Dict"]), ["layer0", "layer1", "layer2"])
        self.assertEqual(list(stacked["PO_Dict"]), ["layer0", "layer1", "layer2"])
        np.testing.assert_allclose(stacked["MSE_Dict"]["layer1"], [np.nan, 0.30], equal_nan=True)
        np.testing.assert_allclose(stacked["PO_Dict"]["layer2"], [np.nan, 0.5], equal_nan=True)
        wrapped = attach_layer_dicts(
            {"k": 2, "po_base": 1e-6, "n_new": 20, "rows": rows, "n_batches": 2}
        )
        self.assertIn("all_layer_backprop", wrapped)
        self.assertEqual(len(wrapped["MSE_Dict"]["layer0"]), 2)

    def test_batch_mse_and_paired_po(self):
        y = np.array([0.0, 1.0, 1.0, 0.0])
        self.assertAlmostEqual(batch_mse(y, y), 0.0)
        self.assertGreater(batch_mse(y, 1.0 - y), 0.0)
        rng = np.random.default_rng(9)
        X0 = rng.normal(size=(60, 3))
        Y0 = (X0[:, 0] > 0).astype(float)
        X1 = rng.normal(size=(60, 3)) + 0.8
        Y1 = (X1[:, 0] > 0).astype(float)
        po, mse = streaming_po_and_mse(X0, Y0, X1, Y1, mu_fn=None, seed=9)
        self.assertGreater(po, 0.0)
        self.assertGreater(mse, 0.0)
        self.assertLess(mse, 1.0)

    def test_cv_fills_mse_on_large_hops(self):
        rng = np.random.default_rng(11)
        n_ref, n_new = 220, 180
        X0 = rng.normal(size=(n_ref, 4))
        Y0 = (X0[:, 0] > 0).astype(float)
        X1 = rng.normal(size=(n_new, 4))
        Y1 = (1.0 - (X1[:, 0] > 0).astype(float))
        X = np.vstack([X0, X1])
        Y = np.concatenate([Y0, Y1])
        registry = DLModelRegistry(epochs=2, patience=2, amp=False, lr=5e-3, verbose=False)
        out = run_layer_freeze_cv(
            X,
            Y,
            n_ref=n_ref,
            batch_size_stream=90,
            hidden_dims=(8, 4),
            online_epochs=1,
            max_batches=2,
            registry=registry,
            n_ref_eval=80,
        )
        self.assertEqual(set(out["MSE_Dict"]), {"layer0", "layer1", "layer2", "layer3"})
        self.assertEqual(len(out["PO_Dict"]["layer1"]), out["n_batches"])
        self.assertIn("board_action", out)
        self.assertTrue(all("mse_stream" in r and "action" in r for r in out["rows"]))
        for r in out["rows"]:
            if r["action"] == ACTION_FREEZE:
                self.assertEqual(len(r["layers"]), 4)
                self.assertTrue(all("mse" in x for x in r["layers"]))
                self.assertIn("i_star_mse", r)
                self.assertTrue(np.isfinite(out["MSE_Dict"][f"layer{r['i_star']}"][r["t"]]))
            else:
                self.assertTrue(r["all_trainable"])


    def test_po_mse_contrast_is_the_board_readout(self):
        self.assertEqual(po_mse_action(False, False), ACTION_KEEP)
        self.assertEqual(po_mse_action(True, False), ACTION_WATCH)
        self.assertEqual(po_mse_action(True, True), ACTION_FREEZE)
        self.assertEqual(po_mse_action(False, True), ACTION_KEEP)
        rows = [
            {"t": 0, "n_new": 20, "po_stream": 1e-6, "mse_stream": 0.10},
            {"t": 1, "n_new": 20, "po_stream": 3e-6, "mse_stream": 0.11},
            {"t": 2, "n_new": 20, "po_stream": 4e-6, "mse_stream": 0.50},
        ]
        info = annotate_po_mse_contrast(rows, po_base=1e-6, mse_base=0.10, n_new=20)
        self.assertEqual(rows[0]["action"], ACTION_KEEP)
        self.assertEqual(rows[1]["action"], ACTION_WATCH)
        self.assertEqual(rows[2]["action"], ACTION_FREEZE)
        self.assertEqual(info["board_action"], ACTION_FREEZE)
        self.assertEqual(info["n_watch"], 1)
        self.assertEqual(info["n_freeze"], 1)


if __name__ == "__main__":
    unittest.main()
