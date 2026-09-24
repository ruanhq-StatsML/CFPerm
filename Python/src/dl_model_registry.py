"""Tabular AnyMLP registry: k layers → model_0 … model_k.

model_i trains only the top i layer groups (i=0 is fully frozen).
Online updates use AdamW / SGD + CosineAnnealingLR.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any, Callable, Optional, Sequence

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset

BuildFn = Callable[[], nn.Module]
FitFn = Callable[..., nn.Module]
PredictFn = Callable[..., np.ndarray]
SaveFn = Callable[..., None]
LoadFn = Callable[[str], nn.Module]
OptimizerFactory = Callable[[nn.Module], torch.optim.Optimizer]
SchedulerFactory = Callable[[torch.optim.Optimizer], Any]


class TabularDataset(Dataset):
    def __init__(self, X, Y, as_sequence: bool = False):
        self.X = torch.as_tensor(np.asarray(X), dtype=torch.float32)
        y = np.asarray(Y)
        if y.dtype == np.float64 or y.dtype == np.float32 or np.issubdtype(y.dtype, np.floating):
            self.Y = torch.as_tensor(y, dtype=torch.float32)
        else:
            self.Y = torch.as_tensor(y)
        if self.Y.ndim > 1 and self.Y.size(-1) == 1:
            self.Y = self.Y.view(-1)
        if as_sequence and self.X.ndim == 2:
            self.X = self.X.unsqueeze(-1)

    def __len__(self):
        return int(self.X.size(0))

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]


def construct_dataloader(
    X,
    Y,
    batch_size: int = 128,
    shuffle: bool = True,
    num_workers: int = 0,
    collate_fn=None,
    as_sequence: bool = False,
):
    ds = TabularDataset(X, Y, as_sequence=as_sequence)
    return DataLoader(
        ds,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        collate_fn=collate_fn,
        pin_memory=torch.cuda.is_available(),
    )


@dataclass(frozen=True)
class DLModelAdapter:
    name: str
    build: BuildFn
    fit: FitFn
    predict: PredictFn
    save: Optional[SaveFn] = None
    load: Optional[LoadFn] = None


class SwiGLUDropIn(nn.Module):
    """SiLU stand-in: true SwiGLU needs paired projections, not a pointwise act."""

    def forward(self, x):
        return torch.nn.functional.silu(x)


ACTIVATIONS = {
    "relu": nn.ReLU,
    "gelu": nn.GELU,
    "silu": nn.SiLU,
    "swiglu": SwiGLUDropIn,
    "tanh": nn.Tanh,
    "leaky_relu": nn.LeakyReLU,
    "identity": nn.Identity,
}


def _get_activation(name):
    if name is None:
        return nn.Identity
    key = str(name).lower()
    if key not in ACTIVATIONS:
        raise KeyError(f"Unsupported activation '{name}'. Available: {sorted(ACTIVATIONS)}")
    return ACTIVATIONS[key]


def _flatten_tabular(X: torch.Tensor) -> torch.Tensor:
    if X.dim() == 3:
        return X.squeeze(1) if X.size(1) == 1 else X.mean(dim=1)
    return X


class AnyMLP(nn.Module):
    def __init__(
        self,
        input_dim,
        output_dim,
        hidden_dims,
        dropout=0.1,
        activation="gelu",
        output_activation=None,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.output_dim = int(output_dim)
        self.hidden_dims = tuple(int(h) for h in hidden_dims)
        self.blocks = nn.ModuleDict()
        dims = [self.input_dim] + list(self.hidden_dims)
        act_cls = _get_activation(activation)
        for i in range(len(dims) - 1):
            self.blocks[f"fc{i}"] = nn.Linear(dims[i], dims[i + 1])
            self.blocks[f"act{i}"] = act_cls()
            self.blocks[f"drop{i}"] = nn.Dropout(float(dropout))
        self.head = nn.Linear(dims[-1], self.output_dim)
        out_cls = _get_activation(output_activation)
        self.out_act = out_cls()

    @property
    def n_hidden(self) -> int:
        return len(self.hidden_dims)

    @property
    def n_layer_groups(self) -> int:
        """k: hidden blocks + head. k+1 models are model_0 … model_k."""
        return self.n_hidden + 1

    def layer_param_groups(self) -> list[list[nn.Parameter]]:
        groups: list[list[nn.Parameter]] = []
        for i in range(self.n_hidden):
            groups.append(list(self.blocks[f"fc{i}"].parameters()))
        groups.append(list(self.head.parameters()))
        return groups

    def forward(self, X):
        X = _flatten_tabular(X)
        for i in range(self.n_hidden):
            X = self.blocks[f"drop{i}"](self.blocks[f"act{i}"](self.blocks[f"fc{i}"](X)))
        return self.out_act(self.head(X))


class AnyTransformer(nn.Module):
    """Feature-as-token encoder for tabular rows. Optional; the board uses AnyMLP."""

    def __init__(
        self,
        input_dim,
        output_dim,
        d_model=64,
        nhead=4,
        num_layers=2,
        dim_ff=128,
        dropout=0.1,
        activation="gelu",
        pooling="mean",
    ):
        super().__init__()
        act = activation.lower()
        if act in {"swiglu", "silu"}:
            act = "gelu"
        self.pooling = pooling
        self.input_proj = nn.Linear(1, d_model)
        self.col_embed = nn.Parameter(torch.zeros(1, int(input_dim), d_model))
        nn.init.trunc_normal_(self.col_embed, std=0.02)
        if pooling == "cls":
            self.cls_token = nn.Parameter(torch.zeros(1, 1, d_model))
            nn.init.trunc_normal_(self.cls_token, std=0.02)
        self.layers = nn.ModuleList(
            [
                nn.TransformerEncoderLayer(
                    d_model=d_model,
                    nhead=nhead,
                    dim_feedforward=dim_ff,
                    dropout=dropout,
                    batch_first=True,
                    activation=act,
                )
                for _ in range(num_layers)
            ]
        )
        self.norm = nn.LayerNorm(d_model)
        self.head = nn.Linear(d_model, output_dim)

    def forward(self, X):
        if X.dim() == 2:
            X = X.unsqueeze(-1)
        elif X.dim() == 3 and X.size(1) == 1 and X.size(-1) != 1:
            X = X.transpose(1, 2)
        X = self.input_proj(X) + self.col_embed[:, : X.size(1), :]
        if self.pooling == "cls":
            cls_tk = self.cls_token.expand(X.size(0), -1, -1)
            X = torch.cat([cls_tk, X], dim=1)
        for layer in self.layers:
            X = layer(X)
        X = self.norm(X)
        if self.pooling == "cls":
            X = X[:, 0]
        elif self.pooling == "last":
            X = X[:, -1]
        else:
            X = X.mean(dim=1)
        return self.head(X)


class FineTuneWrapper(nn.Module):
    def __init__(self, backbone, head, forward_fn=None):
        super().__init__()
        self.backbone = backbone
        self.head = head
        self._forward_fn = forward_fn

    def forward(self, X):
        if self._forward_fn is not None:
            return self._forward_fn(self.backbone, self.head, X)
        Z = (
            self.backbone(X)
            if not isinstance(self.backbone, nn.ModuleDict)
            else self._default_multi(self.backbone, X)
        )
        if Z.dim() > 2:
            Z = Z.mean(dim=1)
        return self.head(Z)

    @staticmethod
    def _default_multi(backbones: nn.ModuleDict, X):
        feats = [backbones[k](v) for k, v in X.items()]
        return torch.cat([f.mean(dim=1) if f.dim() > 2 else f for f in feats], dim=-1)


def _freeze_backbone(backbone, freeze: bool, unfreeze_last_n: int = 0) -> None:
    params = list(backbone.parameters())
    for p in params:
        p.requires_grad = not freeze
    if freeze and unfreeze_last_n > 0:
        n = min(unfreeze_last_n, len(params))
        for p in params[-n:]:
            p.requires_grad = True


def apply_train_top_i(model: AnyMLP, i: int) -> AnyMLP:
    """Train the top i layer groups; i=0 freezes everything."""
    groups = model.layer_param_groups()
    n = len(groups)
    i = int(i)
    for j, ps in enumerate(groups):
        trainable = i > 0 and j >= n - i
        for p in ps:
            p.requires_grad = trainable
    return model


def apply_train_stem(model: AnyMLP) -> AnyMLP:
    """Train only the bottom layer group. Covariate-local stem adapt, not concept freeze."""
    groups = model.layer_param_groups()
    for j, ps in enumerate(groups):
        trainable = j == 0
        for p in ps:
            p.requires_grad = trainable
    return model


def spawn_layer_models(pretrained: AnyMLP) -> list[AnyMLP]:
    """k layers → k+1 clones: model_0 … model_k."""
    k = pretrained.n_layer_groups
    out = []
    for i in range(k + 1):
        m = copy.deepcopy(pretrained)
        apply_train_top_i(m, i)
        m._freeze_i = i  # type: ignore[attr-defined]
        m._registry_name = f"model_{i}"  # type: ignore[attr-defined]
        out.append(m)
    return out


class DLModelRegistry:
    def __init__(
        self,
        device=None,
        lr=1e-3,
        weight_decay=1e-5,
        epochs=8,
        grad_clip=1.0,
        amp=True,
        patience=3,
        optimizer="adamw",
        scheduler="cosine",
        optimizer_factory=None,
        scheduler_factory=None,
        loss_fn=None,
        verbose=False,
    ):
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device)
        self.lr = float(lr)
        self.weight_decay = float(weight_decay)
        self.epochs = int(epochs)
        self.grad_clip = float(grad_clip)
        self.patience = int(patience)
        self.verbose = bool(verbose)
        self.amp = bool(amp)
        self.optimizer = str(optimizer).lower()
        self.scheduler = str(scheduler).lower()
        self._optimizer_factory = optimizer_factory
        self._scheduler_factory = scheduler_factory
        self._loss_fn = loss_fn

    def _default_optimizer_factory(self):
        if self._optimizer_factory is not None:
            return self._optimizer_factory

        def factory(model):
            params = [p for p in model.parameters() if p.requires_grad]
            if not params:
                return None
            if self.optimizer == "sgd":
                return torch.optim.SGD(params, lr=self.lr, momentum=0.9, weight_decay=self.weight_decay)
            return torch.optim.AdamW(params, lr=self.lr, weight_decay=self.weight_decay)

        return factory

    def _default_scheduler_factory(self, n_steps: int | None = None):
        if self._scheduler_factory is not None:
            return self._scheduler_factory
        kind = self.scheduler
        t_max = max(1, int(n_steps or self.epochs))

        def factory(optimizer):
            if optimizer is None:
                return None
            if kind == "cosine":
                return torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=t_max)
            if kind == "step":
                return torch.optim.lr_scheduler.StepLR(optimizer, step_size=max(1, t_max // 3), gamma=0.1)
            if kind == "plateau":
                return torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=2)
            if kind == "onecycle":
                return torch.optim.lr_scheduler.OneCycleLR(optimizer, max_lr=self.lr, total_steps=t_max)
            if kind == "linear":
                return torch.optim.lr_scheduler.LinearLR(
                    optimizer, start_factor=1.0, end_factor=0.01, total_iters=t_max
                )
            if kind == "constant":
                return torch.optim.lr_scheduler.ConstantLR(optimizer, factor=1.0)
            raise ValueError(f"Unknown scheduler: {self.scheduler}")

        return factory

    def _loss(self, task: str):
        if self._loss_fn is not None:
            return self._loss_fn
        if task == "regression":
            return nn.MSELoss()
        if task == "multiclass":
            return nn.CrossEntropyLoss()
        return nn.BCEWithLogitsLoss()

    def _prepare_yb(self, Yb: torch.Tensor, task: str, output_dim: int) -> torch.Tensor:
        if task == "multiclass":
            return Yb.long().view(-1)
        if task == "regression":
            return Yb.float().view(-1, output_dim) if output_dim > 1 else Yb.float().view(-1, 1)
        return Yb.float().view(-1, 1)

    def _train_loop(
        self,
        model,
        train_loader,
        val_loader=None,
        optimizer_factory=None,
        scheduler_factory=None,
        loss_fn=None,
        task="classification",
        epochs=None,
        restore_best=True,
    ):
        model = model.to(self.device)
        epochs = int(epochs or self.epochs)
        trainable = [p for p in model.parameters() if p.requires_grad]
        if not trainable:
            return model
        opt_factory = optimizer_factory or self._default_optimizer_factory()
        optimizer = opt_factory(model)
        sch_factory = scheduler_factory or self._default_scheduler_factory(epochs)
        scheduler = sch_factory(optimizer)
        _loss_fn = loss_fn or self._loss(task)
        use_amp = self.amp and self.device.type == "cuda"
        scaler = torch.amp.GradScaler("cuda", enabled=use_amp)
        output_dim = int(getattr(model, "output_dim", 1) or 1)
        best_state, best_loss, bad = None, float("inf"), 0
        for epoch in range(epochs):
            model.train()
            train_total = 0.0
            n = 0
            for Xb, Yb in train_loader:
                Xb = Xb.to(self.device)
                Yb = self._prepare_yb(Yb.to(self.device), task, output_dim)
                optimizer.zero_grad(set_to_none=True)
                with torch.amp.autocast("cuda", enabled=use_amp):
                    pred = model(Xb)
                    if task != "multiclass" and pred.dim() == 1:
                        pred = pred.view(-1, 1)
                    loss = _loss_fn(pred, Yb)
                scaler.scale(loss).backward()
                if self.grad_clip > 0:
                    scaler.unscale_(optimizer)
                    nn.utils.clip_grad_norm_(trainable, self.grad_clip)
                scaler.step(optimizer)
                scaler.update()
                train_total += float(loss.item()) * Xb.size(0)
                n += Xb.size(0)
            train_loss = train_total / max(1, n)
            if scheduler is not None:
                if isinstance(scheduler, torch.optim.lr_scheduler.ReduceLROnPlateau):
                    scheduler.step(train_loss)
                else:
                    scheduler.step()
            if val_loader is None:
                val_loss = train_loss
            else:
                model.eval()
                val_total = 0.0
                n = 0
                with torch.no_grad():
                    for Xb, Yb in val_loader:
                        Xb = Xb.to(self.device)
                        Yb = self._prepare_yb(Yb.to(self.device), task, output_dim)
                        with torch.amp.autocast("cuda", enabled=use_amp):
                            pred = model(Xb)
                            if task != "multiclass" and pred.dim() == 1:
                                pred = pred.view(-1, 1)
                            loss = _loss_fn(pred, Yb)
                        val_total += float(loss.item()) * Xb.size(0)
                        n += Xb.size(0)
                val_loss = val_total / max(1, n)
            if self.verbose:
                print(f"[{epoch + 1}/{epochs}] train={train_loss:.4f} val={val_loss:.4f}")
            if val_loss < best_loss - 1e-8:
                best_loss = val_loss
                bad = 0
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
                if bad >= self.patience:
                    break
        if restore_best and best_state is not None:
            model.load_state_dict(best_state)
        return model

    def predict_numpy(self, model, X_new, task="classification") -> np.ndarray:
        model.eval()
        with torch.no_grad():
            Xb = X_new if torch.is_tensor(X_new) else torch.as_tensor(np.asarray(X_new), dtype=torch.float32)
            Xb = Xb.to(self.device)
            out = model(Xb)
            if task == "classification":
                if out.size(-1) == 1:
                    out = torch.sigmoid(out.view(-1))
                else:
                    out = torch.softmax(out, dim=-1)
            return out.detach().cpu().numpy()

    def make_any_mlp(
        self,
        *,
        input_dim: int,
        output_dim: int = 1,
        hidden_dims: Sequence[int] = (64, 32),
        dropout: float = 0.1,
        activation: str = "gelu",
        task: str = "classification",
    ) -> DLModelAdapter:
        def build():
            return AnyMLP(input_dim, output_dim, hidden_dims, dropout, activation)

        def fit(model, train_loader, val_loader=None, config=None):
            cfg = config or {}
            return self._train_loop(
                model,
                train_loader,
                val_loader,
                optimizer_factory=cfg.get("optimizer_factory"),
                scheduler_factory=cfg.get("scheduler_factory"),
                loss_fn=cfg.get("loss_fn"),
                task=task,
                epochs=cfg.get("epochs"),
            )

        def predict(fitted, X_new):
            return self.predict_numpy(fitted, X_new, task=task)

        def save(fitted, path):
            torch.save(fitted.state_dict(), path)

        def load(path):
            m = build().to(self.device)
            m.load_state_dict(torch.load(path, map_location=self.device))
            return m

        return DLModelAdapter(name="any_mlp", build=build, fit=fit, predict=predict, save=save, load=load)

    def make_any_transformer_encoder(
        self,
        *,
        input_dim: int,
        output_dim: int = 1,
        d_model: int = 64,
        nhead: int = 4,
        num_layers: int = 2,
        dim_ff: int = 128,
        dropout: float = 0.1,
        pooling: str = "mean",
        task: str = "classification",
    ) -> DLModelAdapter:
        def build():
            return AnyTransformer(
                input_dim,
                output_dim,
                d_model=d_model,
                nhead=nhead,
                num_layers=num_layers,
                dim_ff=dim_ff,
                dropout=dropout,
                pooling=pooling,
            )

        def fit(model, train_loader, val_loader=None, config=None):
            cfg = config or {}
            return self._train_loop(
                model,
                train_loader,
                val_loader,
                optimizer_factory=cfg.get("optimizer_factory"),
                scheduler_factory=cfg.get("scheduler_factory"),
                loss_fn=cfg.get("loss_fn"),
                task=task,
                epochs=cfg.get("epochs"),
            )

        def predict(fitted, X_new):
            return self.predict_numpy(fitted, X_new, task=task)

        def save(fitted, path):
            torch.save(fitted.state_dict(), path)

        def load(path):
            m = build().to(self.device)
            m.load_state_dict(torch.load(path, map_location=self.device))
            return m

        return DLModelAdapter(
            name="any_transformer_encoder",
            build=build,
            fit=fit,
            predict=predict,
            save=save,
            load=load,
        )

    def make_any_finetune(
        self,
        *,
        backbone: nn.Module,
        head: nn.Module,
        freeze_backbone: bool = True,
        unfreeze_last_n: int = 0,
        forward_fn: Optional[Callable] = None,
        task: str = "classification",
        output_dim: int = 1,
    ) -> DLModelAdapter:
        def build():
            _freeze_backbone(backbone, freeze_backbone, unfreeze_last_n)
            for p in head.parameters():
                p.requires_grad = True
            m = FineTuneWrapper(backbone, head, forward_fn)
            m.output_dim = output_dim  # type: ignore[attr-defined]
            return m

        def fit(model, train_loader, val_loader=None, config=None):
            cfg = config or {}
            return self._train_loop(
                model,
                train_loader,
                val_loader,
                optimizer_factory=cfg.get("optimizer_factory"),
                scheduler_factory=cfg.get("scheduler_factory"),
                loss_fn=cfg.get("loss_fn"),
                task=task,
                epochs=cfg.get("epochs"),
            )

        def predict(fitted, X_new):
            return self.predict_numpy(fitted, X_new, task=task)

        def save(fitted, path):
            torch.save(fitted.state_dict(), path)

        def load(path):
            m = build().to(self.device)
            m.load_state_dict(torch.load(path, map_location=self.device))
            return m

        return DLModelAdapter(name="any_finetune", build=build, fit=fit, predict=predict, save=save, load=load)
