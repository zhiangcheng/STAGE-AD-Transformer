"""Multi-task loss used for pretraining and fine-tuning."""

from __future__ import annotations

from collections.abc import Mapping

import torch
import torch.nn as nn


class MultiTaskLoss(nn.Module):
    def __init__(
        self,
        class_weights: torch.Tensor | None = None,
        w_class: float = 1.0,
        w_binary: float = 0.7,
        w_gene: float = 0.8,
        w_tissue: float = 0.6,
        w_cell: float = 0.6,
        w_mechanism: float = 0.5,
        w_recon: float = 0.0,
    ) -> None:
        super().__init__()
        self.w_class = w_class
        self.w_binary = w_binary
        self.w_gene = w_gene
        self.w_tissue = w_tissue
        self.w_cell = w_cell
        self.w_mechanism = w_mechanism
        self.w_recon = w_recon
        self.ce_class = nn.CrossEntropyLoss(weight=class_weights)
        self.ce = nn.CrossEntropyLoss()
        self.bce = nn.BCEWithLogitsLoss()
        self.mse = nn.MSELoss()

    def forward(
        self,
        output: Mapping[str, torch.Tensor],
        batch: Mapping[str, torch.Tensor],
        recon_target: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        class_loss = self.ce_class(output["class_logits"], batch["class_y"])
        binary_loss = self.bce(output["ad_logit"], batch["binary_y"])
        gene_loss = self.ce(output["gene_logits"], batch["target_gene_y"])
        tissue_loss = self.ce(output["tissue_logits"], batch["tissue_y"])
        cell_loss = self.ce(output["cell_logits"], batch["cell_y"])
        mechanism_loss = self.ce(output["mechanism_logits"], batch["mech_y"])
        total = (
            self.w_class * class_loss
            + self.w_binary * binary_loss
            + self.w_gene * gene_loss
            + self.w_tissue * tissue_loss
            + self.w_cell * cell_loss
            + self.w_mechanism * mechanism_loss
        )

        reconstruction_loss = torch.tensor(0.0, device=total.device)
        if recon_target is not None and self.w_recon > 0:
            reconstruction_loss = self.mse(
                output["variant_reconstruction"], recon_target
            )
            total = total + self.w_recon * reconstruction_loss

        parts = {
            "loss_total": float(total.detach().cpu()),
            "loss_class": float(class_loss.detach().cpu()),
            "loss_binary": float(binary_loss.detach().cpu()),
            "loss_gene": float(gene_loss.detach().cpu()),
            "loss_tissue": float(tissue_loss.detach().cpu()),
            "loss_cell": float(cell_loss.detach().cpu()),
            "loss_mechanism": float(mechanism_loss.detach().cpu()),
            "loss_reconstruction": float(reconstruction_loss.detach().cpu()),
        }
        return total, parts
