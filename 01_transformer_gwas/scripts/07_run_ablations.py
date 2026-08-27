#!/usr/bin/env python3
"""Step 07: evaluate fixed-model, inference-time feature ablations."""

from __future__ import annotations

import argparse
import sys
import warnings
from argparse import Namespace
from pathlib import Path

import pandas as pd
import torch

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.dataset import SexRegV3Dataset
from sexreg_ad.train_utils import (
    evaluate_model,
    load_checkpoint,
    make_model,
    model_args_from_checkpoint,
)
from sexreg_ad.utils import ensure_dir


ABLATIONS = [
    "full",
    "no_gwas",
    "no_pip",
    "no_qtl",
    "no_epigenomics",
    "no_single_cell",
    "no_x",
    "no_sequence",
    "no_sex_token",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run inference-time SexReg-AD ablations."
    )
    parser.add_argument("--model", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--outdir", default="results/metrics")
    parser.add_argument("--batch-size", type=int, default=2048)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_dir(args.outdir)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    checkpoint = load_checkpoint(args.model, device)
    base_args = model_args_from_checkpoint(checkpoint)
    scaler = checkpoint.get("scaler")
    if scaler is None:
        warnings.warn(
            "Legacy checkpoint has no training scaler; fitting preprocessing on "
            "the test input. Retrain with this repository for strict reuse.",
            stacklevel=2,
        )

    rows = []
    for ablation in ABLATIONS:
        print(f"Ablation: {ablation}")
        model_values = vars(base_args).copy()
        model_values["no_sex_token"] = ablation == "no_sex_token"
        model = make_model(Namespace(**model_values)).to(device)
        model.load_state_dict(checkpoint["model_state_dict"])

        feature_ablation = (
            [] if ablation in {"full", "no_sex_token"} else [ablation]
        )
        dataset = SexRegV3Dataset(
            args.test,
            scaler=scaler,
            ablate_groups=feature_ablation,
        )
        loader = torch.utils.data.DataLoader(
            dataset, batch_size=args.batch_size, shuffle=False, num_workers=0
        )
        metrics = evaluate_model(model, loader, device)
        metrics["ablation"] = ablation
        rows.append(metrics)

    output = pd.DataFrame(rows)
    output.to_csv(Path(args.outdir) / "ablation_metrics.csv", index=False)
    print(output.sort_values("AUPRC", ascending=False).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
