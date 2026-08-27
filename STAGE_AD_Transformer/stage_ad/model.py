"""Hierarchical sex-aware multi-task Transformer."""

from __future__ import annotations

import torch
import torch.nn as nn


class ScalarFeatureTokenizer(nn.Module):
    """Represent every scalar feature as its own learned token."""

    def __init__(self, n_features: int, d_model: int) -> None:
        super().__init__()
        self.weight = nn.Parameter(torch.randn(n_features, d_model) * 0.02)
        self.bias = nn.Parameter(torch.zeros(n_features, d_model))

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        return (
            features.unsqueeze(-1) * self.weight.unsqueeze(0)
            + self.bias.unsqueeze(0)
        )


class SexRegADTransformerV3(nn.Module):
    """Fuse variant, gene, tissue, cell, and sex-context tokens."""

    def __init__(
        self,
        n_variant_features: int,
        n_gene_features: int,
        n_tissue_features: int,
        n_cell_features: int,
        n_gene_candidates: int = 6,
        n_tissues: int = 9,
        n_cells: int = 7,
        n_classes: int = 5,
        n_mechanisms: int = 8,
        d_model: int = 160,
        n_heads: int = 8,
        n_layers: int = 5,
        dropout: float = 0.1,
        use_sex_token: bool = True,
    ) -> None:
        super().__init__()
        if d_model % n_heads != 0:
            raise ValueError("d_model must be divisible by n_heads")
        if min(
            n_variant_features,
            n_gene_features,
            n_tissue_features,
            n_cell_features,
        ) <= 0:
            raise ValueError("Every feature group must contain at least one feature")

        self.use_sex_token = use_sex_token
        self.n_variant_features = n_variant_features
        self.n_gene_candidates = n_gene_candidates
        self.n_tissues = n_tissues
        self.n_cells = n_cells

        self.variant_tokenizer = ScalarFeatureTokenizer(n_variant_features, d_model)
        self.gene_encoder = nn.Sequential(
            nn.Linear(n_gene_features, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
        )
        self.tissue_encoder = nn.Sequential(
            nn.Linear(n_tissue_features, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
        )
        self.cell_encoder = nn.Sequential(
            nn.Linear(n_cell_features, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
        )
        self.sex_embedding = nn.Embedding(3, d_model)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, d_model))
        self.type_embedding = nn.Embedding(6, d_model)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)
        self.norm = nn.LayerNorm(d_model)

        self.gene_cross = nn.MultiheadAttention(
            d_model, n_heads, dropout=dropout, batch_first=True
        )
        self.tissue_cross = nn.MultiheadAttention(
            d_model, n_heads, dropout=dropout, batch_first=True
        )
        self.cell_cross = nn.MultiheadAttention(
            d_model, n_heads, dropout=dropout, batch_first=True
        )
        self.cross_norm = nn.LayerNorm(d_model)

        self.class_head = self._classification_head(d_model, n_classes, dropout)
        self.ad_head = self._classification_head(d_model, 1, dropout)
        self.mechanism_head = self._classification_head(
            d_model, n_mechanisms, dropout
        )
        self.gene_query = nn.Linear(d_model, d_model)
        self.tissue_query = nn.Linear(d_model, d_model)
        self.cell_query = nn.Linear(d_model, d_model)
        self.variant_reconstruction = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.GELU(),
            nn.Linear(d_model, n_variant_features),
        )

    @staticmethod
    def _classification_head(
        d_model: int, n_outputs: int, dropout: float
    ) -> nn.Sequential:
        return nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, n_outputs),
        )

    def _add_type(self, tokens: torch.Tensor, type_id: int) -> torch.Tensor:
        type_ids = torch.full(
            (tokens.shape[1],),
            type_id,
            dtype=torch.long,
            device=tokens.device,
        )
        return tokens + self.type_embedding(type_ids).unsqueeze(0)

    def forward(
        self,
        variant_x: torch.Tensor,
        gene_x: torch.Tensor,
        tissue_x: torch.Tensor,
        cell_x: torch.Tensor,
        sex_context: torch.Tensor,
        return_attention: bool = False,
    ) -> dict[str, object]:
        batch_size = variant_x.shape[0]
        variant_tokens = self.variant_tokenizer(variant_x)
        gene_tokens = self.gene_encoder(gene_x)
        tissue_tokens = self.tissue_encoder(tissue_x)
        cell_tokens = self.cell_encoder(cell_x)
        sex_token = self.sex_embedding(sex_context).unsqueeze(1)
        if not self.use_sex_token:
            sex_token = torch.zeros_like(sex_token)
        cls_token = self.cls_token.expand(batch_size, -1, -1)

        tokens = torch.cat(
            [
                self._add_type(cls_token, 0),
                self._add_type(sex_token, 1),
                self._add_type(variant_tokens, 2),
                self._add_type(gene_tokens, 3),
                self._add_type(tissue_tokens, 4),
                self._add_type(cell_tokens, 5),
            ],
            dim=1,
        )
        encoded = self.encoder(tokens)

        cls_encoded = encoded[:, 0:1, :]
        offset = 2
        variant_encoded = encoded[
            :, offset : offset + self.n_variant_features, :
        ]
        offset += self.n_variant_features
        gene_encoded = encoded[:, offset : offset + self.n_gene_candidates, :]
        offset += self.n_gene_candidates
        tissue_encoded = encoded[:, offset : offset + self.n_tissues, :]
        offset += self.n_tissues
        cell_encoded = encoded[:, offset : offset + self.n_cells, :]

        gene_context, gene_attention = self.gene_cross(
            cls_encoded,
            gene_encoded,
            gene_encoded,
            need_weights=return_attention,
        )
        tissue_context, tissue_attention = self.tissue_cross(
            cls_encoded,
            tissue_encoded,
            tissue_encoded,
            need_weights=return_attention,
        )
        cell_context, cell_attention = self.cell_cross(
            cls_encoded,
            cell_encoded,
            cell_encoded,
            need_weights=return_attention,
        )
        representation = self.norm(
            self.cross_norm(
                cls_encoded + gene_context + tissue_context + cell_context
            ).squeeze(1)
        )

        output: dict[str, object] = {
            "class_logits": self.class_head(representation),
            "ad_logit": self.ad_head(representation).squeeze(-1),
            "gene_logits": torch.einsum(
                "bd,bkd->bk", self.gene_query(representation), gene_encoded
            ),
            "tissue_logits": torch.einsum(
                "bd,btd->bt", self.tissue_query(representation), tissue_encoded
            ),
            "cell_logits": torch.einsum(
                "bd,bcd->bc", self.cell_query(representation), cell_encoded
            ),
            "mechanism_logits": self.mechanism_head(representation),
            "variant_reconstruction": self.variant_reconstruction(
                variant_encoded.mean(dim=1)
            ),
            "embedding": representation,
        }
        if return_attention:
            output["attention"] = {
                "gene": gene_attention,
                "tissue": tissue_attention,
                "cell": cell_attention,
            }
        return output
