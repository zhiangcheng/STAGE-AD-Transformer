#!/usr/bin/env python3
"""Step 03: masked regulatory pretraining."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.losses import MultiTaskLoss
from sexreg_ad.train_utils import make_loaders, make_model, move_batch
from sexreg_ad.utils import ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pretrain SexReg-AD.")
    parser.add_argument("--train", required=True)
    parser.add_argument("--valid", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--max-train", type=int, default=400_000)
    parser.add_argument("--d-model", type=int, default=160)
    parser.add_argument("--n-heads", type=int, default=8)
    parser.add_argument("--n-layers", type=int, default=5)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--mask-rate", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument("--no-sex-token", action="store_true")
    return parser.parse_args()


def masked_loss(
    model: torch.nn.Module,
    batch: dict[str, torch.Tensor],
    criterion: MultiTaskLoss,
    device: str,
    mask_rate: float,
) -> tuple[torch.Tensor, dict[str, float]]:
    moved = move_batch(batch, device)
    target = moved["variant_x"].clone()
    mask = (torch.rand_like(target) < mask_rate).to(target.dtype)
    masked = target * (1 - mask)
    output = model(
        masked,
        moved["gene_x"],
        moved["tissue_x"],
        moved["cell_x"],
        moved["sex_context"],
    )
    return criterion(output, moved, recon_target=target)


def validation_loss(
    model: torch.nn.Module,
    loader: torch.utils.data.DataLoader,
    criterion: MultiTaskLoss,
    device: str,
    mask_rate: float,
) -> float:
    model.eval()
    losses = []
    with torch.no_grad():
        for batch in loader:
            _, parts = masked_loss(model, batch, criterion, device, mask_rate)
            losses.append(parts["loss_total"])
    return float(np.mean(losses))


def main() -> int:
    args = parse_args()
    if not 0 < args.mask_rate < 1:
        raise ValueError("--mask-rate must be between 0 and 1")
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
        _,
        _,
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
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    criterion = MultiTaskLoss(
        w_class=0.1,
        w_binary=0.1,
        w_gene=0.5,
        w_tissue=0.8,
        w_cell=0.8,
        w_mechanism=0.8,
        w_recon=1.0,
    )

    history = []
    best_validation_loss = float("inf")
    checkpoint_path = checkpoint_dir / "sexreg_ad_pretrained.pt"
    for epoch in range(1, args.epochs + 1):
        model.train()
        training_losses = []
        for batch in tqdm(train_loader, desc=f"pretrain {epoch}"):
            loss, parts = masked_loss(
                model, batch, criterion, device, args.mask_rate
            )
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            training_losses.append(parts["loss_total"])

        valid_loss = validation_loss(
            model, valid_loader, criterion, device, args.mask_rate
        )
        row = {
            "epoch": epoch,
            "train_pretrain_loss": float(np.mean(training_losses)),
            "valid_pretrain_loss": valid_loss,
        }
        history.append(row)
        print(row)
        if valid_loss < best_validation_loss:
            best_validation_loss = valid_loss
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
        metrics_dir / "pretraining_history.csv", index=False
    )
    print(f"Saved best checkpoint to {checkpoint_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
