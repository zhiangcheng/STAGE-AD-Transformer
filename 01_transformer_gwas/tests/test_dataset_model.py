from __future__ import annotations

from argparse import Namespace

import numpy as np
import pandas as pd
import torch

from sexreg_ad.constants import (
    CELL_FEATURES,
    CELL_TYPES,
    CLASS_NAMES,
    GENE_FEATURES,
    MECHANISMS,
    TISSUE_FEATURES,
    TISSUES,
    VARIANT_FEATURES,
)
from sexreg_ad.dataset import SexRegV3Dataset
from sexreg_ad.train_utils import make_model


def model_ready_frame(n_rows: int, offset: float = 0.0) -> pd.DataFrame:
    rows = np.arange(n_rows, dtype=np.float32) + offset
    data: dict[str, object] = {
        "sex_label_id": np.arange(n_rows) % len(CLASS_NAMES),
        "is_positive": np.arange(n_rows) % 2,
        "target_gene_label": np.arange(n_rows) % 6,
        "top_tissue_label": np.arange(n_rows) % len(TISSUES),
        "top_cell_label": np.arange(n_rows) % len(CELL_TYPES),
        "mechanism_label": np.arange(n_rows) % len(MECHANISMS),
        "sex_context_id": np.arange(n_rows) % 3,
    }
    for index, feature in enumerate(VARIANT_FEATURES):
        data[feature] = rows + index
    for candidate in range(6):
        for index, feature in enumerate(GENE_FEATURES):
            data[f"gene{candidate}_{feature}"] = rows + candidate + index
    for tissue in range(len(TISSUES)):
        for index, feature in enumerate(TISSUE_FEATURES):
            data[f"tissue{tissue}_{feature}"] = rows + tissue + index
    for cell in range(len(CELL_TYPES)):
        for index, feature in enumerate(CELL_FEATURES):
            data[f"cell{cell}_{feature}"] = rows + cell + index
    return pd.DataFrame(data)


def test_scaler_is_fitted_on_training_data_only(tmp_path):
    train_path = tmp_path / "train.parquet"
    test_path = tmp_path / "test.parquet"
    model_ready_frame(8).to_parquet(train_path, index=False)
    model_ready_frame(4, offset=100).to_parquet(test_path, index=False)

    train = SexRegV3Dataset(train_path)
    test = SexRegV3Dataset(test_path, scaler=train.scaler_state)

    assert np.allclose(train.variant_x.mean(axis=0), 0, atol=1e-6)
    assert not np.allclose(test.variant_x.mean(axis=0), 0, atol=1e-3)
    assert test.scaler_state == train.scaler_state


def test_model_forward_shapes(tmp_path):
    data_path = tmp_path / "data.parquet"
    model_ready_frame(4).to_parquet(data_path, index=False)
    dataset = SexRegV3Dataset(data_path)
    batch = {
        key: torch.stack([dataset[0][key], dataset[1][key]])
        for key in dataset[0]
    }
    args = Namespace(
        d_model=32,
        n_heads=4,
        n_layers=1,
        dropout=0.0,
        no_sex_token=False,
    )
    model = make_model(args)
    output = model(
        batch["variant_x"],
        batch["gene_x"],
        batch["tissue_x"],
        batch["cell_x"],
        batch["sex_context"],
    )

    assert output["ad_logit"].shape == (2,)
    assert output["class_logits"].shape == (2, len(CLASS_NAMES))
    assert output["gene_logits"].shape == (2, 6)
    assert output["tissue_logits"].shape == (2, len(TISSUES))
    assert output["cell_logits"].shape == (2, len(CELL_TYPES))
    assert output["mechanism_logits"].shape == (2, len(MECHANISMS))
    assert output["variant_reconstruction"].shape == (2, len(VARIANT_FEATURES))
