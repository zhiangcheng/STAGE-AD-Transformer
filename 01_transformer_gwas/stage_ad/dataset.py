"""Dataset and train-only feature scaling for SexReg-AD."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

from .constants import (
    CELL_FEATURES,
    CELL_TYPES,
    GENE_FEATURES,
    TISSUE_FEATURES,
    TISSUES,
    VARIANT_FEATURES,
)


ScalerState = dict[str, dict[str, list[Any]]]


class SexRegV3Dataset(Dataset):
    """Load a model-ready Parquet table.

    If ``scaler`` is omitted, scaling statistics are fitted on this dataset.
    Training code should save ``scaler_state`` and pass it to every validation,
    test, prediction, and ablation dataset.
    """

    def __init__(
        self,
        parquet_path: str | Path,
        *,
        standardize: bool = True,
        scaler: Mapping[str, Mapping[str, Any]] | None = None,
        ablate_groups: Sequence[str] | None = None,
        max_rows: int | None = None,
        seed: int = 1,
    ) -> None:
        self.path = Path(parquet_path)
        self.df = pd.read_parquet(self.path)
        if max_rows is not None and len(self.df) > max_rows:
            self.df = self.df.sample(max_rows, random_state=seed).reset_index(drop=True)
        if self.df.empty:
            raise ValueError(f"Dataset is empty: {self.path}")

        self.ablate_groups = set(ablate_groups or [])
        self.variant_cols = list(VARIANT_FEATURES)
        self.gene_cols = [
            [f"gene{k}_{feature}" for feature in GENE_FEATURES] for k in range(6)
        ]
        self.tissue_cols = [
            [f"tissue{t}_{feature}" for feature in TISSUE_FEATURES]
            for t in range(len(TISSUES))
        ]
        self.cell_cols = [
            [f"cell{c}_{feature}" for feature in CELL_FEATURES]
            for c in range(len(CELL_TYPES))
        ]
        self._validate_columns()

        self.variant_x = self._matrix(self.variant_cols)
        self.gene_x = self._stack(self.gene_cols)
        self.tissue_x = self._stack(self.tissue_cols)
        self.cell_x = self._stack(self.cell_cols)

        self._scaler: dict[str, dict[str, np.ndarray]] | None = None
        if standardize:
            arrays = self._feature_arrays()
            self._scaler = (
                self._fit_scaler(arrays)
                if scaler is None
                else self._load_scaler(scaler, arrays)
            )
            self.variant_x = self._scale("variant", self.variant_x)
            self.gene_x = self._scale("gene", self.gene_x)
            self.tissue_x = self._scale("tissue", self.tissue_x)
            self.cell_x = self._scale("cell", self.cell_x)

        self._ablate()
        n_rows = len(self.df)
        self.sex_context = self.df.get(
            "sex_context_id", pd.Series(np.full(n_rows, 2), index=self.df.index)
        ).astype(int).to_numpy()
        self.class_y = self.df["sex_label_id"].astype(int).to_numpy()
        self.binary_y = self.df["is_positive"].astype(float).to_numpy()
        self.target_gene_y = self.df["target_gene_label"].astype(int).to_numpy()
        self.tissue_y = self.df["top_tissue_label"].astype(int).to_numpy()
        self.cell_y = self.df["top_cell_label"].astype(int).to_numpy()
        self.mech_y = self.df["mechanism_label"].astype(int).to_numpy()

    def _validate_columns(self) -> None:
        feature_columns = (
            self.variant_cols
            + [column for group in self.gene_cols for column in group]
            + [column for group in self.tissue_cols for column in group]
            + [column for group in self.cell_cols for column in group]
        )
        label_columns = [
            "sex_label_id",
            "is_positive",
            "target_gene_label",
            "top_tissue_label",
            "top_cell_label",
            "mechanism_label",
        ]
        missing = sorted(set(feature_columns + label_columns) - set(self.df.columns))
        if missing:
            preview = ", ".join(missing[:12])
            suffix = " ..." if len(missing) > 12 else ""
            raise ValueError(
                f"{self.path} is missing {len(missing)} required columns: "
                f"{preview}{suffix}. See docs/DATA_SCHEMA.md."
            )

    def _matrix(self, columns: Sequence[str]) -> np.ndarray:
        return (
            self.df[list(columns)]
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0)
            .to_numpy(dtype=np.float32)
        )

    def _stack(self, groups: Sequence[Sequence[str]]) -> np.ndarray:
        return np.stack([self._matrix(columns) for columns in groups], axis=1)

    def _feature_arrays(self) -> dict[str, np.ndarray]:
        return {
            "variant": self.variant_x,
            "gene": self.gene_x,
            "tissue": self.tissue_x,
            "cell": self.cell_x,
        }

    @staticmethod
    def _fit_scaler(
        arrays: Mapping[str, np.ndarray],
    ) -> dict[str, dict[str, np.ndarray]]:
        scaler: dict[str, dict[str, np.ndarray]] = {}
        for name, values in arrays.items():
            mean = values.mean(axis=0, keepdims=True, dtype=np.float64).astype(
                np.float32
            )
            std = values.std(axis=0, keepdims=True, dtype=np.float64).astype(
                np.float32
            )
            std[std == 0] = 1.0
            scaler[name] = {"mean": mean, "std": std}
        return scaler

    @staticmethod
    def _load_scaler(
        scaler: Mapping[str, Mapping[str, Any]],
        arrays: Mapping[str, np.ndarray],
    ) -> dict[str, dict[str, np.ndarray]]:
        loaded: dict[str, dict[str, np.ndarray]] = {}
        for name, values in arrays.items():
            if name not in scaler or "mean" not in scaler[name] or "std" not in scaler[name]:
                raise ValueError(f"Scaler is missing statistics for feature group {name!r}")
            mean = np.asarray(scaler[name]["mean"], dtype=np.float32)
            std = np.asarray(scaler[name]["std"], dtype=np.float32)
            expected = (1, *values.shape[1:])
            if mean.shape != expected or std.shape != expected:
                raise ValueError(
                    f"Scaler shape mismatch for {name}: expected {expected}, "
                    f"got mean={mean.shape}, std={std.shape}"
                )
            std = std.copy()
            std[std == 0] = 1.0
            loaded[name] = {"mean": mean, "std": std}
        return loaded

    def _scale(self, name: str, values: np.ndarray) -> np.ndarray:
        if self._scaler is None:
            return values.astype(np.float32)
        stats = self._scaler[name]
        return ((values - stats["mean"]) / stats["std"]).astype(np.float32)

    @property
    def scaler_state(self) -> ScalerState | None:
        """Return checkpoint-safe training preprocessing statistics."""

        if self._scaler is None:
            return None
        return {
            name: {
                "mean": stats["mean"].tolist(),
                "std": stats["std"].tolist(),
            }
            for name, stats in self._scaler.items()
        }

    def _zero_variant_cols(self, keys: Sequence[str]) -> None:
        indices = [
            i
            for i, column in enumerate(self.variant_cols)
            if any(key in column for key in keys)
        ]
        if indices:
            self.variant_x[:, indices] = 0

    def _ablate(self) -> None:
        groups = self.ablate_groups
        if "no_gwas" in groups:
            self._zero_variant_cols(["Z_", "P_", "Beta", "neglog10", "sex_delta"])
        if "no_pip" in groups:
            self._zero_variant_cols(["PIP"])
        if "no_x" in groups:
            self._zero_variant_cols(["x_"])
        if "no_sequence" in groups:
            self._zero_variant_cols(
                ["deepsea", "enformer", "sei", "cadd", "conservation"]
            )
        if "no_qtl" in groups:
            self.gene_x[:, :, [1, 2, 3, 4, 5]] = 0
            self.tissue_x[:, :, [0, 1, 2]] = 0
        if "no_epigenomics" in groups:
            self.tissue_x[:, :, [3, 4]] = 0
            self.cell_x[:, :, [2]] = 0
        if "no_single_cell" in groups:
            self.cell_x[:] = 0

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        return {
            "variant_x": torch.tensor(self.variant_x[index], dtype=torch.float32),
            "gene_x": torch.tensor(self.gene_x[index], dtype=torch.float32),
            "tissue_x": torch.tensor(self.tissue_x[index], dtype=torch.float32),
            "cell_x": torch.tensor(self.cell_x[index], dtype=torch.float32),
            "sex_context": torch.tensor(self.sex_context[index], dtype=torch.long),
            "class_y": torch.tensor(self.class_y[index], dtype=torch.long),
            "binary_y": torch.tensor(self.binary_y[index], dtype=torch.float32),
            "target_gene_y": torch.tensor(
                self.target_gene_y[index], dtype=torch.long
            ),
            "tissue_y": torch.tensor(self.tissue_y[index], dtype=torch.long),
            "cell_y": torch.tensor(self.cell_y[index], dtype=torch.long),
            "mech_y": torch.tensor(self.mech_y[index], dtype=torch.long),
        }
