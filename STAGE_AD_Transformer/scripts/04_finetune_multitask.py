#!/usr/bin/env python3
"""Step 04: supervised multi-task fine-tuning and held-out evaluation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.constants import CLASS_NAMES
from sexreg_ad.losses import MultiTaskLoss
from sexreg_ad.train_utils import (
    evaluate_model,
    load_checkpoint,
    make_loaders,
    make_model,
    move_batch,
)
from sexreg_ad.utils import ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fine-tune SexReg-AD.")
    parser.add_argument("--train", required=True)
    parser.add_argument("--valid", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--pretrained", default="")
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--max-train", type=int, default=500_000)
    parser.add_argument("--d-model", type=int, default=160)
    parser.add_argument("--n-heads", type=int, default=8)
    parser.add_argument("--n-layers", type=int, default=5)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=8e-5)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument("--no-sex-token", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    checkpoint_dir = Path(args.outdir) / "checkpoints"
    metrics_dir = Path(args.outdir) / "metrics"
    ensure_dir(checkpoint_dir)
    ensure_dir(metrics_dir)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    (
        _,
        _,
        _,
        train_loader,
        valid_loader,
        test_loader,
        counts,
        scaler,
    ) = make_loaders(
        args.train,
        args.valid,
        args.test,
        args.batch_size,
        args.max_train,
        args.seed,
    )

    model = make_model(args).to(device)
    if args.pretrained:
        checkpoint = load_checkpoint(args.pretrained, device)
        model.load_state_dict(checkpoint["model_state_dict"], strict=False)
        print(f"Loaded pretrained weights from {args.pretrained}")

    class_weights = torch.tensor(
        1 / np.maximum(counts, 1), dtype=torch.float32, device=device
    )
    class_weights = class_weights / class_weights.mean()
    criterion = MultiTaskLoss(class_weights=class_weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=max(1, args.epochs)
    )

    history = []
    best_score = float("-inf")
    checkpoint_path = checkpoint_dir / "sexreg_ad_finetuned.pt"
    for epoch in range(1, args.epochs + 1):
        model.train()
        training_losses = []
        for batch in tqdm(train_loader, desc=f"finetune {epoch}"):
            moved = move_batch(batch, device)
            output = model(
                moved["variant_x"],
                moved["gene_x"],
                moved["tissue_x"],
                moved["cell_x"],
                moved["sex_context"],
            )
            loss, parts = criterion(output, moved)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            training_losses.append(parts["loss_total"])

        scheduler.step()
        validation = evaluate_model(model, valid_loader, device)
        row = {
            "epoch": epoch,
            "train_loss": float(np.mean(training_losses)),
            **validation,
        }
        history.append(row)
        print(row)

        score = validation["AUPRC"]
        if not np.isfinite(score):
            score = -row["train_loss"]
        if epoch == 1 or score > best_score:
            best_score = score
            torch.save(
                {
                    "format_version": 2,
                    "model_state_dict": model.state_dict(),
                    "args": vars(args),
                    "scaler": scaler,
                },
                checkpoint_path,
            )

    pd.DataFrame(history).to_csv(
        metrics_dir / "finetuning_history.csv", index=False
    )
    best_checkpoint = load_checkpoint(checkpoint_path, device)
    model.load_state_dict(best_checkpoint["model_state_dict"])
    test_metrics = evaluate_model(model, test_loader, device)
    test_metrics["model"] = "SexReg-AD Transformer v3"
    pd.DataFrame([test_metrics]).to_csv(
        metrics_dir / "transformer_test_metrics.csv", index=False
    )
    print(test_metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
