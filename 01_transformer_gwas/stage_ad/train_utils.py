"""Shared model, loader, checkpoint, and evaluation helpers."""

from __future__ import annotations

from argparse import Namespace
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch.utils.data import DataLoader, WeightedRandomSampler

from .constants import (
    CELL_FEATURES,
    CELL_TYPES,
    CLASS_NAMES,
    GENE_FEATURES,
    MECHANISMS,
    TISSUE_FEATURES,
    TISSUES,
    VARIANT_FEATURES,
)
from .dataset import ScalerState, SexRegV3Dataset
from .metrics import binary_metrics, multiclass_metrics, ranking_metrics
from .model import SexRegADTransformerV3


DEFAULT_MODEL_ARGS = {
    "d_model": 160,
    "n_heads": 8,
    "n_layers": 5,
    "dropout": 0.1,
    "no_sex_token": False,
}


def load_checkpoint(
    path: str | Path, device: str | torch.device
) -> dict[str, Any]:
    """Load a repository checkpoint with safe-mode fallback for older PyTorch."""

    try:
        return torch.load(path, map_location=device, weights_only=True)
    except TypeError:
        return torch.load(path, map_location=device)


def make_model(args: Namespace, *, use_sex_token: bool = True) -> SexRegADTransformerV3:
    """Construct a model from argparse/checkpoint-compatible attributes."""

    return SexRegADTransformerV3(
        n_variant_features=len(VARIANT_FEATURES),
        n_gene_features=len(GENE_FEATURES),
        n_tissue_features=len(TISSUE_FEATURES),
        n_cell_features=len(CELL_FEATURES),
        n_gene_candidates=6,
        n_tissues=len(TISSUES),
        n_cells=len(CELL_TYPES),
        n_classes=len(CLASS_NAMES),
        n_mechanisms=len(MECHANISMS),
        d_model=args.d_model,
        n_heads=args.n_heads,
        n_layers=args.n_layers,
        dropout=args.dropout,
        use_sex_token=use_sex_token and not getattr(args, "no_sex_token", False),
    )


def model_args_from_checkpoint(checkpoint: Mapping[str, Any]) -> Namespace:
    """Recover model construction arguments, including legacy checkpoints."""

    values = dict(DEFAULT_MODEL_ARGS)
    values.update(
        {key.replace("-", "_"): value for key, value in checkpoint.get("args", {}).items()}
    )
    return Namespace(**values)


def make_loaders(
    train_path: str,
    valid_path: str,
    test_path: str,
    batch_size: int,
    max_train: int | None = None,
    seed: int = 1,
    ablate_groups: list[str] | None = None,
) -> tuple[
    SexRegV3Dataset,
    SexRegV3Dataset,
    SexRegV3Dataset,
    DataLoader,
    DataLoader,
    DataLoader,
    np.ndarray,
    ScalerState | None,
]:
    """Create splits while fitting preprocessing on the training set only."""

    train_ds = SexRegV3Dataset(
        train_path,
        max_rows=max_train,
        seed=seed,
        ablate_groups=ablate_groups,
    )
    scaler = train_ds.scaler_state
    valid_ds = SexRegV3Dataset(
        valid_path,
        scaler=scaler,
        ablate_groups=ablate_groups,
    )
    test_ds = SexRegV3Dataset(
        test_path,
        scaler=scaler,
        ablate_groups=ablate_groups,
    )

    counts = np.bincount(train_ds.class_y, minlength=len(CLASS_NAMES))
    class_weights = 1 / np.maximum(counts, 1)
    generator = torch.Generator().manual_seed(seed)
    sampler = WeightedRandomSampler(
        class_weights[train_ds.class_y],
        num_samples=len(train_ds),
        replacement=True,
        generator=generator,
    )
    loader_kwargs = {"batch_size": batch_size, "num_workers": 0}
    train_loader = DataLoader(train_ds, sampler=sampler, **loader_kwargs)
    valid_loader = DataLoader(valid_ds, shuffle=False, **loader_kwargs)
    test_loader = DataLoader(test_ds, shuffle=False, **loader_kwargs)
    return (
        train_ds,
        valid_ds,
        test_ds,
        train_loader,
        valid_loader,
        test_loader,
        counts,
        scaler,
    )


def move_batch(
    batch: Mapping[str, torch.Tensor], device: str | torch.device
) -> dict[str, torch.Tensor]:
    return {key: value.to(device) for key, value in batch.items()}


def evaluate_model(
    model: SexRegADTransformerV3,
    loader: DataLoader,
    device: str | torch.device,
) -> dict[str, Any]:
    model.eval()
    binary_labels: list[np.ndarray] = []
    binary_scores: list[np.ndarray] = []
    class_labels: list[np.ndarray] = []
    class_probabilities: list[np.ndarray] = []
    gene_labels: list[np.ndarray] = []
    gene_scores: list[np.ndarray] = []
    tissue_labels: list[np.ndarray] = []
    tissue_scores: list[np.ndarray] = []
    cell_labels: list[np.ndarray] = []
    cell_scores: list[np.ndarray] = []

    with torch.no_grad():
        for batch in loader:
            moved = move_batch(batch, device)
            output = model(
                moved["variant_x"],
                moved["gene_x"],
                moved["tissue_x"],
                moved["cell_x"],
                moved["sex_context"],
            )
            binary_labels.append(moved["binary_y"].cpu().numpy())
            binary_scores.append(torch.sigmoid(output["ad_logit"]).cpu().numpy())
            class_labels.append(moved["class_y"].cpu().numpy())
            class_probabilities.append(
                torch.softmax(output["class_logits"], dim=1).cpu().numpy()
            )
            gene_labels.append(moved["target_gene_y"].cpu().numpy())
            gene_scores.append(output["gene_logits"].cpu().numpy())
            tissue_labels.append(moved["tissue_y"].cpu().numpy())
            tissue_scores.append(output["tissue_logits"].cpu().numpy())
            cell_labels.append(moved["cell_y"].cpu().numpy())
            cell_scores.append(output["cell_logits"].cpu().numpy())

    metrics: dict[str, Any] = {}
    metrics.update(
        binary_metrics(np.concatenate(binary_labels), np.concatenate(binary_scores))
    )
    metrics.update(
        multiclass_metrics(
            np.concatenate(class_labels), np.vstack(class_probabilities)
        )
    )
    metrics.update(
        ranking_metrics(
            np.concatenate(gene_labels),
            np.vstack(gene_scores),
            np.concatenate(tissue_labels),
            np.vstack(tissue_scores),
            np.concatenate(cell_labels),
            np.vstack(cell_scores),
        )
    )
    return metrics
