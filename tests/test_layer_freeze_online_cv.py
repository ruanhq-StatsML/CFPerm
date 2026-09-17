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
from layer_freeze_cv import run_layer_freeze_cv  # noqa: E402
from streaming_po_risk import pack_ref_new, po_risk, streaming_po_risk  # noqa: E402


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
        a = streaming_po_risk(X0, Y0, X1, Y1)
        X, Y, T = pack_ref_new(X0, Y0, X1, Y1)
        b = po_risk(X, Y, T)
        self.assertAlmostEqual(a, b, places=12)


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

    def test_cv_returns_k_plus_one_and_i_star(self):
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
        self.assertEqual(len(out["rows"]), 2 * 4)
        self.assertIn(out["rows"][0]["i_star"], {0, 1, 2, 3})
        self.assertTrue(all(r["n_new"] == 90 for r in out["rows"]))


if __name__ == "__main__":
    unittest.main()
