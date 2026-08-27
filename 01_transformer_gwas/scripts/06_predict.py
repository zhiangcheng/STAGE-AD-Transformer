#!/usr/bin/env python3
"""Step 06: run checkpointed SexReg-AD inference."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.constants import CELL_TYPES, CLASS_NAMES, MECHANISMS, TISSUES
from sexreg_ad.dataset import SexRegV3Dataset
from sexreg_ad.train_utils import (
    load_checkpoint,
    make_model,
    model_args_from_checkpoint,
    move_batch,
)
from sexreg_ad.utils import ensure_dir


METADATA_COLUMNS = [
    "variant_id",
    "rsid",
    "CHR",
    "BP",
    "nearest_gene",
    "sex_label",
    "sex_label_id",
    "is_positive",
    "P_total",
    "P_female",
    "P_male",
    "P_interaction",
    "Beta_female",
    "Beta_male",
    "chr_type",
    "xci_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate SexReg-AD predictions.")
    parser.add_argument("--model", required=True, help="Fine-tuned checkpoint.")
    parser.add_argument("--input", required=True, help="Model-ready Parquet file.")
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--batch-size", type=int, default=2048)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    prediction_dir = Path(args.outdir) / "predictions"
    ensure_dir(prediction_dir)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    checkpoint = load_checkpoint(args.model, device)
    model_args = model_args_from_checkpoint(checkpoint)
    model = make_model(model_args).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    scaler = checkpoint.get("scaler")
    if scaler is None:
        warnings.warn(
            "Legacy checkpoint has no training scaler; fitting preprocessing on "
            "the prediction input. Retrain with this repository for strict reuse.",
            stacklevel=2,
        )
    dataset = SexRegV3Dataset(args.input, scaler=scaler)
    missing_metadata = sorted(set(METADATA_COLUMNS) - set(dataset.df.columns))
    gene_name_columns = [f"gene{k}_name" for k in range(6)]
    missing_gene_names = sorted(set(gene_name_columns) - set(dataset.df.columns))
    if missing_metadata or missing_gene_names:
        missing = ", ".join(missing_metadata + missing_gene_names)
        raise ValueError(f"Prediction input is missing metadata columns: {missing}")

    loader = torch.utils.data.DataLoader(
        dataset, batch_size=args.batch_size, shuffle=False, num_workers=0
    )
    metadata = dataset.df[METADATA_COLUMNS].copy()
    prediction_batches = []

    with torch.no_grad():
        for batch in tqdm(loader, desc="predict-v3"):
            moved = move_batch(batch, device)
            output = model(
                moved["variant_x"],
                moved["gene_x"],
                moved["tissue_x"],
                moved["cell_x"],
                moved["sex_context"],
            )
            class_prob = torch.softmax(output["class_logits"], dim=1).cpu().numpy()
            ad_prob = torch.sigmoid(output["ad_logit"]).cpu().numpy()
            gene_prob = torch.softmax(output["gene_logits"], dim=1).cpu().numpy()
            tissue_prob = torch.softmax(output["tissue_logits"], dim=1).cpu().numpy()
            cell_prob = torch.softmax(output["cell_logits"], dim=1).cpu().numpy()
            mechanism_prob = (
                torch.softmax(output["mechanism_logits"], dim=1).cpu().numpy()
            )

            predicted_class_id = class_prob.argmax(axis=1)
            prediction = pd.DataFrame(
                {
                    "sexreg_ad_score": ad_prob,
                    "predicted_class_id": predicted_class_id,
                    "predicted_class": [
                        CLASS_NAMES[index] for index in predicted_class_id
                    ],
                    "predicted_gene_index": gene_prob.argmax(axis=1),
                    "predicted_tissue": [
                        TISSUES[index] for index in tissue_prob.argmax(axis=1)
                    ],
                    "predicted_cell_type": [
                        CELL_TYPES[index] for index in cell_prob.argmax(axis=1)
                    ],
                    "predicted_mechanism": [
                        MECHANISMS[index]
                        for index in mechanism_prob.argmax(axis=1)
                    ],
                }
            )
            for index, class_name in enumerate(CLASS_NAMES):
                prediction[f"prob_class_{class_name}"] = class_prob[:, index]
            prediction_batches.append(prediction)

    predictions = pd.concat(prediction_batches, ignore_index=True)
    output = pd.concat([metadata.reset_index(drop=True), predictions], axis=1)
    gene_names = dataset.df[gene_name_columns].to_numpy()
    gene_indices = output["predicted_gene_index"].to_numpy(dtype=int)
    output[gene_name_columns] = gene_names
    output["predicted_target_gene"] = gene_names[
        np.arange(len(gene_names)), gene_indices
    ]
    output["female_context_score"] = (
        output["prob_class_female_biased"]
        + output["prob_class_sex_interaction"]
    )
    output["male_context_score"] = (
        output["prob_class_male_biased"] + output["prob_class_sex_interaction"]
    )
    output["delta_sex_score"] = (
        output["female_context_score"] - output["male_context_score"]
    )

    csv_path = prediction_dir / "predictions.csv"
    parquet_path = prediction_dir / "predictions.parquet"
    output.to_csv(csv_path, index=False)
    output.to_parquet(parquet_path, index=False)
    print(output.head().to_string(index=False))
    print(f"Wrote {len(output):,} predictions to {prediction_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
